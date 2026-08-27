## ---------------- Gate 6: Between-person mixture (person switch) ----------------
Gate6_personSwitch <- nimbleCode({
        # prevalence of Tree A at the person level
        alpha_person ~ dnorm(0, sd = 2)
        
        # persons
        for (n in 1:N) {
                theta_dir[n] ~ dnorm(0, sd=1)
                kappa_ext[n] ~ dnorm(0, sd=1)
                SP[n] ~ dbern(ilogit(alpha_person))   # T[n]=1 => Tree A on all items
        }
        
        # items: free thresholds, with anchor for bB1
        for (i in 1:I) {
                bA1_raw[i] ~ dnorm(0, sd = 5)
                gapA[i]    ~ dexp(1.0)       # > 0
        }
        bA1_bar <- mean(bA1_raw[1:I])
        for (i in 1:I) {
                bA1[i] <- bA1_raw[i] - bA1_bar
                bA2[i] <- bA1[i] + gapA[i]   # inherits centering from bA1
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
        
        # likelihood
        for (n in 1:N) for (i in 1:I) {
                p_gate[n,i] <- SP[n]   # hard switch per person
                Y[n,i] ~ dMix3catTree(theta_dir[n], kappa_ext[n],
                                      bA1[i], bA2[i], bB1[i], bB2[i],
                                      p_gate[n,i])
        }
})

## ---------- Log-lik at posterior-mean params (Gate 6: person switch) ----------
ll_at_posterior_mean_G6 <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        TH   <- colMeans(get_vec(M, "theta_dir", N))  # length N
        KAP  <- colMeans(get_vec(M, "kappa_ext", N))  # length N
        bA1  <- colMeans(get_vec(M, "bA1", I))        # length I
        bA2  <- colMeans(get_vec(M, "bA2", I))        # length I (deterministic = bA1+gapA in model)
        bB1  <- colMeans(get_vec(M, "bB1", I))        # length I
        bB2  <- colMeans(get_vec(M, "bB2", I))        # length I
        pg   <- colMeans(get_vec(M, "SP", N))         # soft gate: E[SP_n|y] in [0,1]
        
        eps <- 1e-12
        out <- numeric(N * I); k <- 1L
        for (i in 1:I) for (n in 1:N) {
                ## Tree A probs
                pA1 <- plogis(TH[n]  - bA1[i])
                pA2 <- plogis(TH[n]  - bA2[i])
                pA_yes <- pA1 * pA2
                pA_no  <- 1 - pA1
                pA_prh <- 1 - pA_no - pA_yes
                
                ## Tree B probs
                pB1 <- plogis(KAP[n] - bB1[i])   # Extreme
                pB2 <- plogis(TH[n]  - bB2[i])   # Yes | Extreme
                pB_prh <- 1 - pB1
                pB_yes <- pB1 * pB2
                pB_no  <- 1 - pB_prh - pB_yes
                
                ## Mixture at posterior-mean gate
                p_no  <- pg[n] * pA_no  + (1 - pg[n]) * pB_no
                p_prh <- pg[n] * pA_prh + (1 - pg[n]) * pB_prh
                p_yes <- pg[n] * pA_yes + (1 - pg[n]) * pB_yes
                
                y  <- Y[n, i]                         # 1..3 = {no, perhaps, yes}
                pr <- c(p_no, p_prh, p_yes)[y]
                pr <- min(max(pr, eps), 1 - eps)
                out[k] <- log(pr); k <- k + 1L
        }
        out
}

## ---------- Response-level log-lik matrix S x (N*I) for WAIC/LOO ----------
llmat_gate6_personSwitch <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        TH  <- get_vec(M, "theta_dir", N)   # S x N
        KAP <- get_vec(M, "kappa_ext", N)   # S x N
        bA1 <- get_vec(M, "bA1", I)         # S x I
        bA2 <- get_vec(M, "bA2", I)         # S x I  (deterministic in model but saved via get_vec)
        bB1 <- get_vec(M, "bB1", I)         # S x I
        bB2 <- get_vec(M, "bB2", I)         # S x I
        SP  <- get_vec(M, "SP", N)          # S x N (0/1 by draw)
        
        S <- nrow(TH)
        out <- matrix(NA_real_, nrow = S, ncol = N * I)
        colnames(out) <- paste0("Y[", rep(1:N, each = I), ",", rep(1:I, times = N), "]")
        
        eps <- 1e-12
        k <- 1L
        for (i in 1:I) {
                for (n in 1:N) {
                        # Tree A (vectors length S)
                        pA1 <- plogis(TH[, n]  - bA1[, i])
                        pA2 <- plogis(TH[, n]  - bA2[, i])
                        pA_yes <- pA1 * pA2
                        pA_no  <- 1 - pA1
                        pA_prh <- 1 - pA_no - pA_yes
                        
                        # Tree B
                        pB1 <- plogis(KAP[, n] - bB1[, i])
                        pB2 <- plogis(TH[, n]  - bB2[, i])
                        pB_prh <- 1 - pB1
                        pB_yes <- pB1 * pB2
                        pB_no  <- 1 - pB_prh - pB_yes
                        
                        # Mixture with per-draw hard gate SP (0/1), treated as weight
                        w <- SP[, n]                       # 0 or 1 per draw
                        p_no  <- w * pA_no  + (1 - w) * pB_no
                        p_prh <- w * pA_prh + (1 - w) * pB_prh
                        p_yes <- w * pA_yes + (1 - w) * pB_yes
                        
                        y  <- Y[n, i]
                        pr <- cbind(p_no, p_prh, p_yes)[, y]
                        pr <- pmin(pmax(pr, eps), 1 - eps)
                        out[, k] <- log(pr)
                        k <- k + 1L
                }
        }
        out
}

## ---------- End-to-end wrapper (fixes LOO r_eff via chain_id) ----------
fit_Gate6_personSwitch <- function(Y, niter = 4000, nburn = 2000, thin = 10, nchains = 4,
                                   auto_mm = TRUE, k_threshold = 0.7) {
        fit <- build_and_fit_gate6_personSwitch(Gate6_personSwitch, Y, niter, nburn, thin, nchains)
        
        M <- as_mat(fit$samples)
        
        # chain_id for r_eff
        chain_id <- if (inherits(fit$samples, "mcmc.list")) {
                S_per_chain <- vapply(fit$samples, nrow, 0L)
                rep(seq_along(S_per_chain), times = S_per_chain)
        } else if (is.matrix(fit$samples)) rep(1L, nrow(fit$samples)) else rep(1L, nrow(M))
        
        # Response-level log-lik
        ll <- llmat_gate6_personSwitch(M, Y)
        
        # WAIC + LOO (with r_eff; auto moment-matching optional)
        ic <- info_criteria(ll, chain_id = chain_id, auto_mm = auto_mm, 
                            k_threshold = k_threshold,mm_post_draws = as_mat(fit$samples))
        
        # DIC via posterior-mean parameters (soft gate)
        llhat <- ll_at_posterior_mean_G6(M, Y)
        dic   <- dic_from_loglik_and_mean(ll, llhat)
        
        list(
                model = "Gate6_personSwitch",
                N = fit$N, I = fit$I,
                samples = fit$samples,
                IC = append(dic, ic),
                param_summaries = summarize_params(fit$samples, fit$N, fit$I, model = "G6")
        )
}

## ---------- Build/fit ----------
build_and_fit_gate6_personSwitch <- function(code, Y, niter = 4000, nburn = 2000,
                                             thin = 10, nchains = 4, inits = NULL) {
        N <- nrow(Y); I <- ncol(Y)
        constants <- list(N = N, I = I)
        data <- list(Y = Y)
        
        if (is.null(inits)) {
                inits <- replicate(nchains, simplify = FALSE, list(
                        alpha_person = rnorm(1, 0, 1),
                        theta_dir = rnorm(N, 0, 0.5),
                        kappa_ext = rnorm(N, 0, 0.5),
                        SP = rbinom(N, 1, 0.5),
                        bA1_raw = rnorm(I, 0, 1),
                        gapA = rexp(I, 1),              # positive
                        # bA2 is deterministic (bA1 + gapA)
                        bB1_raw = rnorm(I, 0, 1),
                        #bB2 = bA1+gapA
                        tau = rbeta(I,3,3)
                        #mu_d2 = rnorm(1, 0.5, 0.5),
                        #sigma_d2 = abs(rt(1, df = 3)) * 0.5 + 0.05,  # positive, roughly half-t
                        #d2 = rnorm(I, 0.5, 0.5)        # can be negative
                        # bB2 is deterministic (bA2 - d2)
                ))
        }
        
        Rmodel <- nimbleModel(code, constants = constants, data = data,
                              inits = inits[[1]], check = FALSE)
        Cmodel <- compileNimble(Rmodel)
        
        # Monitors: include SP and all thresholds used by ll functions
        mons <- c("alpha_person", "theta_dir", "kappa_ext", "SP",
                  "bA1", "bA2", "bB1", "bB2", "tau")
        
        conf <- configureMCMC(Rmodel, monitors = mons, thin = thin, enableWAIC = FALSE)
        MCMC <- buildMCMC(conf)
        CMCMC <- compileNimble(MCMC, project = Rmodel)
        
        samp <- runMCMC(CMCMC, niter = niter, nburnin = nburn, thin = thin,
                        nchains = nchains, samplesAsCodaMCMC = TRUE, setSeed = TRUE)
        
        list(samples = samp, N = N, I = I)
}
