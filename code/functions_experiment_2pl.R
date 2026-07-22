Gate2PL_marginal <- nimbleCode({
        # --- Gate hyperparameters & priors
        alpha0 ~ dnorm(0, sd = 5)
        
        sigma_alpha ~ T(dt(0, sigma = 1, df = 3), 0, )  # half-t(3,1)
        sigma_beta  ~ T(dt(0, sigma = 1, df = 3), 0, )
        
        # Person gate effects (centered to sum to zero)
        for (n in 1:N) {
                alpha_raw[n] ~ dnorm(0, sd = sigma_alpha)
        }
        alpha_bar <- mean(alpha_raw[1:N])
        for (n in 1:N) {
                alpha_gate[n] <- alpha_raw[n] - alpha_bar
        }
        
        # Item gate effects (centered to sum to zero)
        for (i in 1:I) {
                beta_raw[i] ~ dnorm(0, sd = sigma_beta)
        }
        beta_bar <- mean(beta_raw[1:I])
        for (i in 1:I) {
                beta_gate[i] <- beta_raw[i] - beta_bar
        }
        
        # --- Persons for the two Rasch trees
        for (n in 1:N) {
                theta_dir[n]  ~ dnorm(0, sd = 1)  # Tree A person
                kappa_ext[n]  ~ dnorm(0, sd = 1)  # Tree B person
        }
        
        # --- Tree A item difficulties (centered) + positive gap
        for (i in 1:I) {
                bA1_raw[i] ~ dnorm(0, sd = 5)
                gapA[i]    ~ dexp(1.0)               # >0
        }
        bA1_bar <- mean(bA1_raw[1:I])
        for (i in 1:I) {
                bA1[i] <- bA1_raw[i] - bA1_bar
                bA2[i] <- bA1[i] + gapA[i]
        }
        
        ## --- Tree B first node: center bB1 ---
        for (i in 1:I) {
                bB1_raw[i] ~ dnorm(0, sd = 5)
        }
        bB1_bar <- mean(bB1_raw[1:I])
        for (i in 1:I) {
                bB1[i] <- bB1_raw[i] - bB1_bar
        }
        
        ## --- Tree B second node: inside (bA1, bA2) ---
        for (i in 1:I) {
                tau[i]  ~ dbeta(3, 3) 
                bB2[i]  <- bA1[i] + tau[i] * gapA[i]
        }
        
        
        # --- Likelihood with marginalized gate
        for (n in 1:N) for (i in 1:I) {
                eta_gate[n,i] <- alpha0 + alpha_gate[n] + beta_gate[i]
                p_gate[n,i]   <- ilogit(eta_gate[n,i])
                
                # dMix3catTree mixes Tree A and Tree B outcome probabilities with weight p_gate
                Y[n,i] ~ dMix3catTree(theta_dir[n], kappa_ext[n],
                                      bA1[i], bA2[i], bB1[i], bB2[i],
                                      p_gate[n,i])
        }
})

## ---------- Log-lik at posterior-mean params (marginalized 2PL gate) ----------
ll_at_posterior_mean <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        
        TH   <- colMeans(get_vec(M, "theta_dir", N))      # length N
        KAP  <- colMeans(get_vec(M, "kappa_ext", N))      # length N
        bA1  <- colMeans(get_vec(M, "bA1", I))            # length I
        bB1  <- colMeans(get_vec(M, "bB1", I))            # length I
        bA2 <- colMeans(get_vec(M, "bA2", I))
        bB2 <- colMeans(get_vec(M, "bB2", I))
        a0   <- mean(M[, "alpha0"])                       # scalar
        aN   <- colMeans(get_vec(M, "alpha_gate", N))     # length N (centered)
        bI   <- colMeans(get_vec(M, "beta_gate", I))      # length I (centered)
        
        eps <- 1e-12
        out <- numeric(N * I); k <- 1L
        
        for (i in 1:I) for (n in 1:N) {
                # Gate weight
                pg <- plogis(a0 + aN[n] + bI[i])
                
                # Tree A (theta-only)
                pA1 <- plogis(TH[n] - bA1[i]); pA2 <- plogis(TH[n] - bA2[i])
                pA_yes <- pA1 * pA2; pA_no <- 1 - pA1; pA_prh <- 1 - pA_no - pA_yes
                
                # Tree B (kappa + theta)
                pB1 <- plogis(KAP[n] - bB1[i]); pB2 <- plogis(TH[n] - bB2[i])
                pB_prh <- 1 - pB1; pB_yes <- pB1 * pB2; pB_no <- 1 - pB_prh - pB_yes
                
                # Mixture
                p_no  <- pg * pA_no  + (1 - pg) * pB_no
                p_prh <- pg * pA_prh + (1 - pg) * pB_prh
                p_yes <- pg * pA_yes + (1 - pg) * pB_yes
                
                y  <- Y[n, i]
                pr <- c(p_no, p_prh, p_yes)[y]
                pr <- min(max(pr, eps), 1 - eps)
                out[k] <- log(pr); k <- k + 1L
        }
        out
}

## ---------- Per-draw log-lik matrix: S x (N*I) (marginalized 2PL gate) ----------
llmat_gate2_2pl <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        
        TH   <- get_vec(M, "theta_dir", N)        # S x N
        KAP  <- get_vec(M, "kappa_ext", N)        # S x N
        bA1  <- get_vec(M, "bA1", I)              # S x I
        bA2  <- get_vec(M, "bA2", I)
        bB1  <- get_vec(M, "bB1", I)              # S x I
        bB2  <- get_vec(M, "bB2", I)
        a0   <- M[, "alpha0"]                     # S
        aN   <- get_vec(M, "alpha_gate", N)       # S x N (centered)
        bI   <- get_vec(M, "beta_gate", I)        # S x I (centered)
        
        S <- nrow(TH)
        out <- matrix(NA_real_, nrow = S, ncol = N * I)
        colnames(out) <- paste0("Y[", rep(1:N, each = I), ",", rep(1:I, times = N), "]")
        
        eps <- 1e-12; k <- 1L
        for (i in 1:I) for (n in 1:N) {
                pg <- plogis(a0 + aN[, n] + bI[, i])         # S
                
                # Tree A
                pA1 <- plogis(TH[, n] - bA1[, i]); pA2 <- plogis(TH[, n] - bA2[, i])
                pA_yes <- pA1 * pA2; pA_no <- 1 - pA1; pA_prh <- 1 - pA_no - pA_yes
                
                # Tree B
                pB1 <- plogis(KAP[, n] - bB1[, i]); pB2 <- plogis(TH[, n] - bB2[, i])
                pB_prh <- 1 - pB1; pB_yes <- pB1 * pB2; pB_no <- 1 - pB_prh - pB_yes
                
                # Mixture
                p_no  <- pg * pA_no  + (1 - pg) * pB_no
                p_prh <- pg * pA_prh + (1 - pg) * pB_prh
                p_yes <- pg * pA_yes + (1 - pg) * pB_yes
                
                y  <- Y[n, i]
                pr <- cbind(p_no, p_prh, p_yes)[, y]
                pr <- pmin(pmax(pr, eps), 1 - eps)
                out[, k] <- log(pr); k <- k + 1L
        }
        out
}

## ---------- Build & fit ----------
build_and_fit_gate2_2pl <- function(code, Y,
                                    niter = 4000, nburn = 2000, thin = 10, nchains = 4,
                                    inits = NULL) {
        N <- nrow(Y); I <- ncol(Y)
        constants <- list(N = N, I = I)
        data <- list(Y = Y)
        
        if (is.null(inits)) {
                inits <- replicate(nchains, simplify = FALSE, list(
                        alpha0 = rnorm(1, 0, 1),
                        sigma_alpha = abs(rt(1, df = 3)) + 0.2,
                        sigma_beta  = abs(rt(1, df = 3)) + 0.2,
                        alpha_raw = rnorm(N, 0, 0.5),
                        beta_raw  = rnorm(I, 0, 0.5),
                        theta_dir = rnorm(N, 0, 0.7),
                        kappa_ext = rnorm(N, 0, 0.7),
                        bA1_raw = rnorm(I, 0, 1.0),
                        gapA    = pmax(rexp(I, 1.0), 1e-4),   # strictly > 0
                        bB1_raw = rnorm(I, 0, 1.0),
                        tau     = pmin(pmax(rbeta(I, 3, 3), 1e-3), 1 - 1e-3)
                        # DO NOT set bA1, bA2, bB1, bB2 (deterministic)
                ))
        }
        
        
        Rmodel <- nimbleModel(code, constants = constants, data = data, inits = inits[[1]], check = FALSE)
        Cmodel <- compileNimble(Rmodel)
        
        # Monitor derived nodes (centered) and deterministics (bA1, bA2, bB1, bB2)
        mons <- c("alpha0", "alpha_gate", "beta_gate", "sigma_alpha", "sigma_beta",
                  "theta_dir", "kappa_ext",
                  "bA1", "bA2", "bB1", "bB2", "tau")
        
        conf <- configureMCMC(Rmodel, monitors = mons, thin = thin, enableWAIC = FALSE)
        MCMC <- buildMCMC(conf)
        CMCMC <- compileNimble(MCMC, project = Rmodel)
        
        samp <- runMCMC(CMCMC, niter = niter, nburnin = nburn, thin = thin, nchains = nchains,
                        samplesAsCodaMCMC = TRUE, setSeed = TRUE)
        
        list(samples = samp, N = N, I = I)
}

## ---------- End-to-end wrapper ----------
fit_gate2_2pl <- function(Y,
                          niter = 4000, nburn = 2000, thin = 10, nchains = 4,
                          auto_mm = TRUE, k_threshold = 0.7) {
        
        fit <- build_and_fit_gate2_2pl(Gate2PL_marginal, Y, niter, nburn, thin, nchains)
        M <- as_mat(fit$samples)
        
        chain_id <- if (inherits(fit$samples, "mcmc.list")) {
                S_per <- vapply(fit$samples, nrow, 0L); rep(seq_along(S_per), S_per)
        } else rep(1L, nrow(M))
        
        ll   <- llmat_gate2_2pl(M, Y)
        ic   <- info_criteria(ll, chain_id = chain_id, auto_mm = auto_mm,
                              k_threshold = k_threshold, mm_post_draws = as_mat(fit$samples))
        llhat <- ll_at_posterior_mean(M, Y)
        dic   <- dic_from_loglik_and_mean(ll, llhat)
        
        list(model = "Gate2PL_marginal", N = fit$N, I = fit$I, samples = fit$samples,
             IC = append(dic, ic),
             param_summaries = summarize_params(fit$samples, fit$N, fit$I, model = "G2"))
}
