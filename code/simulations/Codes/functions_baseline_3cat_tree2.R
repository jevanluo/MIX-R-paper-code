BaselineB_code <- nimbleCode({
        for (n in 1:N) {
                theta_dir[n] ~ dnorm(0, sd = 1)
                kappa_ext[n] ~ dnorm(0, sd = 1)
        }
        for (i in 1:I) {
                bB1_raw[i] ~ dnorm(0, sd = 5)
        }
        bB1_bar <- mean(bB1_raw[1:I])
        for (i in 1:I) {
                bB1[i] <- bB1_raw[i] - bB1_bar
                bB2[i] ~ dnorm(0, sd = 5) 
        }
        for (n in 1:N) for (i in 1:I) {
                pB1[n,i]  <- ilogit(kappa_ext[n] - bB1[i])     # Extreme
                pB2[n,i]  <- ilogit(theta_dir[n] - bB2[i])     # Yes|Extreme
                p_perhaps[n,i] <- 1 - pB1[n,i]
                p_yes[n,i]     <- pB1[n,i] * pB2[n,i]
                p_no[n,i]      <- 1 - p_perhaps[n,i] - p_yes[n,i]
                probs[n,i,1] <- p_no[n,i]
                probs[n,i,2] <- p_perhaps[n,i]
                probs[n,i,3] <- p_yes[n,i]
                Y[n,i] ~ dcat(prob = probs[n,i,1:3])
        }
})

## ---------- utilities ----------
ilogit <- function(x) 1/(1+exp(-x))
log_mean_exp <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }

## ---------- Log-lik at posterior-mean params (for DIC if you want it) ----------
# M: your fit object with stored draws (use your existing get_vec helper)
# Y: N x I matrix with categories coded 1..3 (for {no, perhaps, yes})
ll_at_posterior_mean_B <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        TH  <- colMeans(get_vec(M, "theta_dir", N))   # length N
        KAP <- colMeans(get_vec(M, "kappa_ext", N))   # length N
        bB1 <- colMeans(get_vec(M, "bB1", I))         # length I
        bB2 <- colMeans(get_vec(M, "bB2", I))         # length I
        
        eps <- 1e-12
        out <- numeric(N * I); k <- 1L
        for (i in 1:I) for (n in 1:N) {
                pB1 <- plogis(KAP[n] - bB1[i])             # Extreme
                pB2 <- plogis(TH[n]  - bB2[i])             # Yes | Extreme
                p_perhaps <- 1 - pB1
                p_yes     <- pB1 * pB2
                p_no      <- 1 - p_perhaps - p_yes
                y  <- Y[n, i]                               # 1..3 = {no, perhaps, yes}
                pr <- c(p_no, p_perhaps, p_yes)[y]
                pr <- min(max(pr, eps), 1 - eps)
                out[k] <- log(pr); k <- k + 1L
        }
        out
}

llmat_baselineB <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        # Draw matrices (S x dim)
        TH  <- get_vec(M, "theta_dir", N)   # S x N
        KAP <- get_vec(M, "kappa_ext", N)   # S x N
        bB1 <- get_vec(M, "bB1", I)         # S x I
        bB2 <- get_vec(M, "bB2", I)         # S x I
        
        S <- nrow(TH)
        out <- matrix(NA_real_, nrow = S, ncol = N * I)
        colnames(out) <- paste0("Y[", rep(1:N, each = I), ",", rep(1:I, times = N), "]")
        
        eps <- 1e-12
        k <- 1L
        for (i in 1:I) {
                for (n in 1:N) {
                        pB1 <- plogis(KAP[, n] - bB1[, i])           # S
                        pB2 <- plogis(TH[, n]  - bB2[, i])           # S
                        p_perhaps <- 1 - pB1
                        p_yes     <- pB1 * pB2
                        p_no      <- 1 - p_perhaps - p_yes
                        y  <- Y[n, i]                                 # scalar 1..3
                        pr <- cbind(p_no, p_perhaps, p_yes)[, y]      # length S
                        pr <- pmin(pmax(pr, eps), 1 - eps)
                        out[, k] <- log(pr)
                        k <- k + 1L
                }
        }
        out
}

fit_BaselineB <- function(Y, niter = 4000, nburn = 2000, thin = 10, nchains = 4,
                          auto_mm = TRUE, k_threshold = 0.7) {
        # 1) Fit
        fit <- build_and_fit_baselineB(BaselineB_code, Y, niter, nburn, thin, nchains)
        
        # 2) Draws matrix
        M <- as_mat(fit$samples)
        
        # 3) chain_id for r_eff (fixes loo warning)
        chain_id <- NULL
        if (inherits(fit$samples, "mcmc.list")) {
                S_per_chain <- vapply(fit$samples, nrow, 0L)
                chain_id <- rep(seq_along(S_per_chain), times = S_per_chain)
        } else if (is.matrix(fit$samples)) {
                chain_id <- rep(1L, nrow(fit$samples))
        } else {
                chain_id <- rep(1L, nrow(M))
        }
        
        # 4) Response-level log-lik matrix
        ll <- llmat_baselineB(M, Y)
        
        # 5) WAIC + LOO with r_eff (and optional moment matching)
        ic <- info_criteria(ll, chain_id = chain_id, auto_mm = auto_mm, k_threshold = k_threshold)
        
        # 6) DIC via log-lik at posterior-mean params
        llhat <- ll_at_posterior_mean_B(M, Y)
        dic   <- dic_from_loglik_and_mean(ll, llhat)
        
        # 7) Return
        list(
                model = "BaselineB",
                N = fit$N, I = fit$I,
                samples = fit$samples,
                IC = append(dic, ic),  # pD, DIC, lppd, p_waic, elpd_waic, WAIC, LOO_*, etc.
                param_summaries = summarize_params(fit$samples, fit$N, fit$I, "B")
        )
}

## ---------- Build/fit ----------
build_and_fit_baselineB <- function(code, Y, niter=4000, nburn=2000, thin=10, nchains=4, inits=NULL) {
        N <- nrow(Y); I <- ncol(Y)
        constants <- list(N=N, I=I)
        data <- list(Y = Y)
        if (is.null(inits)) {
                # broad defaults; nimble ignores extras not used by a model
                inits <- replicate(nchains, simplify=FALSE,
                                   list(theta_dir = rnorm(N),kappa_ext = rnorm(N),
                                        bB1 = rnorm(I), bB2 = rnorm(I))
                )
        }
        Rmodel <- nimbleModel(code, constants=constants, data=data, inits=inits[[1]], check=FALSE)
        Cmodel <- compileNimble(Rmodel)
        conf   <- configureMCMC(Rmodel, monitors =  c("theta_dir","kappa_ext","bB1","bB2"),
                                thin=thin, enableWAIC = TRUE)
        MCMC   <- buildMCMC(conf); CMCMC <- compileNimble(MCMC, project=Rmodel)
        samp   <- runMCMC(CMCMC, niter=niter, nburnin=nburn, thin=thin, nchains=nchains,
                          samplesAsCodaMCMC=TRUE, setSeed=TRUE)
        list(samples = samp, N=N, I=I)
}