BaselineA_code <- nimbleCode({
        for (n in 1:N) theta_dir[n] ~ dnorm(0, sd = 1)
        for (i in 1:I) {
                bA1_raw[i] ~ dnorm(0, sd = 5)
                gapA[i]    ~ dexp(1.0)       # > 0
        }
        bA1_bar <- mean(bA1_raw[1:I])
        for (i in 1:I) {
                bA1[i] <- bA1_raw[i] - bA1_bar
                bA2[i] <- bA1[i] + gapA[i]   # inherits centering from bA1
        }
        for (n in 1:N) for (i in 1:I) {
                pA1[n,i] <- ilogit(theta_dir[n] - bA1[i])
                pA2[n,i] <- ilogit(theta_dir[n] - bA2[i])
                p_yes[n,i] <- pA1[n,i] * pA2[n,i]
                p_no[n,i]  <- 1 - pA1[n,i]
                p_perhaps[n,i] <- 1 - p_no[n,i] - p_yes[n,i]
                probs[n,i,1] <- p_no[n,i]
                probs[n,i,2] <- p_perhaps[n,i]
                probs[n,i,3] <- p_yes[n,i]
                Y[n,i] ~ dcat(prob = probs[n,i,1:3])
        }
})



## ---------- Log-lik at posterior-mean params (for DIC if you want it) ----------
# M: your fit object with stored draws (use your existing get_vec helper)
# Y: N x I matrix with categories coded 1..3 (for {no, perhaps, yes})
ll_at_posterior_mean_A <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y)
        TH  <- colMeans(get_vec(M, "theta_dir", N))
        bA1 <- colMeans(get_vec(M, "bA1", I))
        bA2 <- colMeans(get_vec(M, "bA2", I))
        eps <- 1e-12
        
        out <- numeric(N * I); k <- 1L
        for (i in 1:I) for (n in 1:N) {
                pA1 <- ilogit(TH[n] - bA1[i])
                pA2 <- ilogit(TH[n] - bA2[i])
                p_no     <- 1 - pA1
                p_yes    <- pA1 * pA2
                p_perhaps<- 1 - p_no - p_yes
                y  <- Y[n, i]
                pr <- c(p_no, p_perhaps, p_yes)[y]
                pr <- min(max(pr, eps), 1 - eps)
                out[k] <- log(pr); k <- k + 1L
        }
        out
}

## ---------- Pointwise log-likelihood matrices ----------
# returns loglik matrix: S x (N*I)
llmat_baselineA <- function(M, Y) {
        N <- nrow(Y); I <- ncol(Y); S <- nrow(M)
        TH <- get_vec(M, "theta_dir", N)         # SxN
        bA1 <- get_vec(M, "bA1", I)              # SxI
        bA2 <- get_vec(M, "bA2", I)              # SxI
        if (is.null(TH) || is.null(bA1) || is.null(bA2)) stop("Missing parameters in samples (A).")
        eps <- 1e-12
        out <- matrix(NA_real_, S, N*I)
        k <- 1L
        for (i in 1:I) {
                pA1 <- plogis(TH - bA1[, i])           # SxN
                pA2 <- plogis(TH - bA2[, i])
                p_no <- 1 - pA1
                p_yes <- pA1 * pA2
                p_maybe <- 1 - p_no - p_yes
                for (n in 1:N) {
                        y <- Y[n,i]
                        pr <- if (y==1L) p_no[,n] else if (y==2L) p_maybe[,n] else p_yes[,n]
                        pr <- pmin(pmax(pr, eps), 1-eps)
                        out[, k] <- log(pr)
                        k <- k + 1L
                }
        }
        out
}
## ---------- End-to-end wrappers ----------
fit_BaselineA <- function(Y, niter=10000, nburn=2000, thin=10, nchains=4) {
        fit <- build_and_fit_baselineA(BaselineA_code, Y, niter, nburn, thin, nchains)
        M   <- as_mat(fit$samples)
        ll  <- llmat_baselineA(M, Y)
        ic  <- info_criteria(ll)
        llhat <- ll_at_posterior_mean_A(M, Y)
        dic <- dic_from_ll(ll, llhat)
        list(
                model = "BaselineA",
                N = fit$N, I = fit$I,
                samples = fit$samples,
                IC = c(dic, ic),                    # pD, DIC, WAIC, LOO, etc.
                param_summaries = summarize_params(fit$samples, fit$N, fit$I, "A")
        )
}

## ---------- Build/fit ----------
build_and_fit_baselineA <- function(code, Y, niter=4000, nburn=2000, thin=10, nchains=4, inits=NULL) {
        N <- nrow(Y); I <- ncol(Y)
        constants <- list(N=N, I=I)
        data <- list(Y = Y)
        if (is.null(inits)) {
                # broad defaults; nimble ignores extras not used by a model
                inits <- replicate(nchains, simplify=FALSE,
                                   list(theta_dir = rnorm(N),
                                        bA1 = rnorm(I), bA2 = rnorm(I))
                )
        }
        Rmodel <- nimbleModel(code, constants=constants, data=data, inits=inits[[1]], check=FALSE)
        Cmodel <- compileNimble(Rmodel)
        conf   <- configureMCMC(Rmodel, monitors =  c("theta_dir","bA1","bA2"),
                                thin=thin, enableWAIC = TRUE)
        MCMC   <- buildMCMC(conf); CMCMC <- compileNimble(MCMC, project=Rmodel)
        samp   <- runMCMC(CMCMC, niter=niter, nburnin=nburn, thin=thin, nchains=nchains,
                          samplesAsCodaMCMC=TRUE, setSeed=TRUE)
        list(samples = samp, N=N, I=I)
}