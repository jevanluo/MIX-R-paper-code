# alpha0_hat: scalar
# alpha_gate_hat: length N (centered, sum≈0)
# beta_gate_hat : length I (centered, sum≈0)
prevalence_from_gate_point <- function(alpha0_hat, alpha_gate_hat, beta_gate_hat) {
        N <- length(alpha_gate_hat); I <- length(beta_gate_hat)
        eta  <- outer(alpha_gate_hat, beta_gate_hat, "+") + alpha0_hat  # N x I
        p    <- plogis(eta)
        list(
                person  = as.numeric(rowMeans(p)),   # length N
                item    = as.numeric(colMeans(p)),   # length I
                overall = mean(p)                    # scalar
        )
}


# Logistic and derivatives
.sigmoid  <- function(x) 1/(1 + exp(-x))
.sigmoid2 <- function(x) { s <- .sigmoid(x); s*(1 - s)*(1 - 2*s) }

# Requires your helper get_vec(M, name, dim) that returns an S x dim matrix.
# M can be an mcmc.list or a matrix of draws.
prevalence_from_gate_draws <- function(M, N, I) {
        if (inherits(M, "mcmc.list")) {
                M <- do.call(rbind, lapply(M, as.matrix))
        } else {
                M <- as.matrix(M)
        }
        
        alpha0  <- M[, "alpha0"]                         # S
        aN      <- get_vec(M, "alpha_gate", N)           # S x N
        bI      <- get_vec(M, "beta_gate",  I)           # S x I
        S       <- nrow(M)
        
        stopifnot(ncol(aN) == N, ncol(bI) == I)
        
        P_person <- matrix(NA_real_, S, N)
        P_item   <- matrix(NA_real_, S, I)
        P_over   <- numeric(S)
        
        for (s in 1:S) {
                eta <- outer(aN[s, ], bI[s, ], "+") + alpha0[s]  # N x I
                p   <- plogis(eta)
                P_person[s, ] <- rowMeans(p)
                P_item[s, ]   <- colMeans(p)
                P_over[s]     <- mean(p)
        }
        
        qsum <- function(x) c(mean = mean(x),
                              q025 = as.numeric(quantile(x, .025)),
                              q500 = as.numeric(quantile(x, .50)),
                              q975 = as.numeric(quantile(x, .975)))
        
        person_df <- as.data.frame(t(apply(P_person, 2, qsum)))
        person_df$id <- 1:N
        item_df   <- as.data.frame(t(apply(P_item,   2, qsum)))
        item_df$item <- 1:I
        overall   <- qsum(P_over)
        
        list(
                person  = person_df[, c("id","mean","q025","q500","q975")],
                item    = item_df[,   c("item","mean","q025","q500","q975")],
                overall = overall
        )
}

# Requires your helper get_vec(M, name, dim) that returns an S x dim matrix.
# M can be an mcmc.list or a matrix of draws.
prevalence_from_gate_draws_plug_in <- function(M, N, I) {
        # unify draws
        if (inherits(M, "mcmc.list")) {
                M <- do.call(rbind, lapply(M, as.matrix))
        } else {
                M <- as.matrix(M)
        }
        
        # pull draws
        alpha0 <- M[, "alpha0"]                    # S
        aN     <- get_vec(M, "alpha_gate", N)      # S x N
        bI     <- get_vec(M, "beta_gate",  I)      # S x I
        S      <- nrow(M)
        stopifnot(ncol(aN) == N, ncol(bI) == I)
        
        # -------- probability-scale functional (for uncertainty) --------
        P_person <- matrix(NA_real_, S, N)
        P_item   <- matrix(NA_real_, S, I)
        P_over   <- numeric(S)
        
        for (s in seq_len(S)) {
                eta <- outer(aN[s, ], bI[s, ], "+") + alpha0[s]   # N x I
                p   <- plogis(eta)
                P_person[s, ] <- rowMeans(p)
                P_item[s, ]   <- colMeans(p)
                P_over[s]     <- mean(p)
        }
        
        # -------- plug-in on the logit scale (bias-reduced point estimate) --------
        alpha0_bar <- mean(alpha0)
        a_bar      <- colMeans(aN)   # length N
        b_bar      <- colMeans(bI)   # length I
        
        eta_bar      <- outer(a_bar, b_bar, "+") + alpha0_bar  # N x I
        p_person_plg <- rowMeans(plogis(eta_bar))              # length N
        p_item_plg   <- colMeans(plogis(eta_bar))              # length I
        p_over_plg   <- mean(plogis(eta_bar))                  # scalar
        
        # -------- posterior summaries for the probability-scale functional --------
        qsum <- function(x) c(
                mean_draw = mean(x),
                q025 = as.numeric(stats::quantile(x, .025)),
                q500 = as.numeric(stats::quantile(x, .50)),
                q975 = as.numeric(stats::quantile(x, .975))
        )
        
        person_q <- as.data.frame(t(apply(P_person, 2, qsum)))  # N x 4
        item_q   <- as.data.frame(t(apply(P_item,   2, qsum)))  # I x 4
        overall_q <- qsum(P_over)
        
        # -------- assemble data frames --------
        person_df <- transform(
                person_q,
                id   = seq_len(N),
                mean = p_person_plg          # <- plug-in point estimate (bias-reduced)
        )[ , c("id","mean","mean_draw","q025","q500","q975")]
        
        item_df <- transform(
                item_q,
                item = seq_len(I),
                mean = p_item_plg            # <- plug-in point estimate (bias-reduced)
        )[ , c("item","mean","mean_draw","q025","q500","q975")]
        
        overall <- c(
                mean      = p_over_plg,      # <- plug-in point estimate (bias-reduced)
                overall_q # mean_draw, q025, q500, q975
        )
        
        list(
                person  = person_df,
                item    = item_df,
                overall = overall
        )
}



# ---- Delta-method prevalence (overall, per-person, per-item) ----
prevalence_from_gate_draws_delta <- function(M, N, I, include_cross = FALSE) {
        if (inherits(M, "mcmc.list")) M <- do.call(rbind, lapply(M, as.matrix)) else M <- as.matrix(M)
        S <- nrow(M)
        
        a0   <- as.numeric(M[, "alpha0"])
        A    <- get_vec(M, "alpha_gate", N)   # S x N (persons)
        B    <- get_vec(M, "beta_gate",  I)   # S x I (items)
        
        # ----- (1) Posterior-draw functional (baseline) -----
        P_person <- P_item <- matrix(NA_real_, S, max(N, I))
        P_over   <- numeric(S)
        for (s in seq_len(S)) {
                eta <- outer(A[s, ], B[s, ], "+") + a0[s]           # N x I
                p   <- .sigmoid(eta)
                P_person[s, 1:N] <- rowMeans(p)
                P_item[s, 1:I]   <- colMeans(p)
                P_over[s]        <- mean(p)
        }
        qsum <- function(x) c(mean_draw = mean(x), q025 = unname(quantile(x,.025)),
                              q500 = unname(quantile(x,.5)), q975 = unname(quantile(x,.975)))
        person_draw <- as.data.frame(t(apply(P_person[, 1:N, drop=FALSE], 2, qsum))); person_draw$id <- 1:N
        item_draw   <- as.data.frame(t(apply(P_item[,   1:I, drop=FALSE], 2, qsum))); item_draw$item <- 1:I
        overall_draw <- qsum(P_over)
        
        # ----- (2) Plug-in + Delta method -----
        a0_bar <- mean(a0);   a0_var <- var(a0)
        A_bar  <- colMeans(A); A_var <- apply(A, 2, var)
        B_bar  <- colMeans(B); B_var <- apply(B, 2, var)
        
        # Covariances with alpha0
        cov_a0_A <- apply(A, 2, function(z) cov(a0, z))
        cov_a0_B <- apply(B, 2, function(z) cov(a0, z))
        
        # Optional Cov(a_n, b_i) across draws (N x I)
        if (include_cross) {
                A_c <- sweep(A, 2, A_bar, "-")
                B_c <- sweep(B, 2, B_bar, "-")
                cov_AB <- crossprod(A_c, B_c) / (S - 1)  # N x I
        } else {
                cov_AB <- matrix(0, N, I)
        }
        
        # Means and variances of η_ni
        MU <- outer(A_bar, B_bar, "+") + a0_bar                   # N x I
        V  <- (a0_var + outer(A_var, rep(1,I)) + outer(rep(1,N), B_var) +
                       2*outer(cov_a0_A, rep(1,I)) + 2*outer(rep(1,N), cov_a0_B) + 2*cov_AB)  # N x I
        
        # Delta-corrected E[σ(η_ni)]
        plug_ni  <- .sigmoid(MU)
        delta_ni <- plug_ni + 0.5 * .sigmoid2(MU) * V
        
        # Assemble outputs
        person_delta <- data.frame(
                id     = 1:N,
                mean   = rowMeans(plug_ni),      # plug-in (for reference)
                delta  = rowMeans(delta_ni),     # delta-corrected (use this)
                mean_draw = person_draw$mean_draw,
                q025   = person_draw$q025,
                q500   = person_draw$q500,
                q975   = person_draw$q975
        )
        item_delta <- data.frame(
                item   = 1:I,
                mean   = colMeans(plug_ni),
                delta  = colMeans(delta_ni),
                mean_draw = item_draw$mean_draw,
                q025   = item_draw$q025,
                q500   = item_draw$q500,
                q975   = item_draw$q975
        )
        overall_delta <- c(
                mean      = mean(plug_ni),
                delta     = mean(delta_ni),
                mean_draw = overall_draw["mean_draw"],
                q025      = overall_draw["q025"],
                q500      = overall_draw["q500"],
                q975      = overall_draw["q975"]
        )
        
        list(person = person_delta, item = item_delta, overall = overall_delta)
}


# Recover posterior Pr(Tree A) per (n,i) and optional hard classes
# M : posterior draws as a matrix (e.g., as_mat(fit$samples))
# Y : observed N x I matrix (only used to infer N,I if needed)
# threshold : hard classify to 1 if posterior mean >= threshold
# simulate_alloc : if TRUE (marginal model), simulate SPI draws from p to get a posterior distribution over classes
recover_gate_classes_G2 <- function(M, Y = NULL, threshold = 0.5,
                                    simulate_alloc = FALSE, seed = NULL) {
        if (!is.null(seed)) set.seed(seed)
        stopifnot(is.matrix(M))
        
        has_SPI <- any(grepl("^SPI\\[", colnames(M)))
        # Infer N,I from available blocks if Y not given
        if (!is.null(Y)) {
                N <- nrow(Y); I <- ncol(Y)
        } else {
                # try alpha_gate and beta_gate
                ag <- try(get_vec(M, "alpha_gate", 1), silent = TRUE)
                bg <- try(get_vec(M, "beta_gate", 1),  silent = TRUE)
                if (!inherits(ag, "try-error")) {
                        N <- ncol(get_vec(M, "alpha_gate", 1))  # will be 1, but next line fixes
                        N <- ncol(get_vec(M, "alpha_gate", sum(grepl("^alpha_gate\\[", colnames(M)))))
                } else if (has_SPI) {
                        # count from SPI columns
                        spi_cols <- sum(grepl("^SPI\\[", colnames(M)))
                        # try to guess N and I by matching theta_dir/bA1 if present
                        N <- ncol(get_vec(M, "theta_dir", sum(grepl("^theta_dir\\[", colnames(M)))))
                        I <- spi_cols / N
                } else stop("Please provide Y (to get N and I).")
                if (is.null(I)) {
                        I <- ncol(get_vec(M, "beta_gate", sum(grepl("^beta_gate\\[", colnames(M)))))
                }
        }
        
        res <- list(N = N, I = I)
        
        if (has_SPI) {
                # ----- SPI model: read SPI draws directly -----
                SPI_draws <- get_vec_disc(M, "SPI", N * I)   # S x (N*I)
                S <- nrow(SPI_draws)
                SPI_arr <- array(SPI_draws, dim = c(S, N, I))
                p_mean  <- apply(SPI_arr, c(2,3), mean)
                p_q025  <- apply(SPI_arr, c(2,3), function(x) quantile(x, 0.025))
                p_q975  <- apply(SPI_arr, c(2,3), function(x) quantile(x, 0.975))
                hard    <- (p_mean >= threshold) * 1L
                res$posterior_prob_A <- p_mean
                res$ci95 <- list(q025 = p_q025, q975 = p_q975)
                res$hard_class <- hard
                # also report posterior SD as uncertainty
                res$sd <- apply(SPI_arr, c(2,3), sd)
                res$source <- "SPI"
                return(res)
        }
        
        # ----- Marginal model: compute p_gate = logistic(alpha0 + a_n + b_i) per draw -----
        a0 <- as.numeric(M[, "alpha0"])                    # S
        aN <- get_vec(M, "alpha_gate", N)                  # S x N
        bI <- get_vec(M, "beta_gate", I)                   # S x I
        S  <- nrow(aN)
        
        # p_{s,n,i} = logistic(a0_s + aN_{s,n} + bI_{s,i})
        # We'll compute posterior mean and optionally simulate allocations
        p_mean <- matrix(NA_real_, N, I)
        p_q025 <- matrix(NA_real_, N, I)
        p_q975 <- matrix(NA_real_, N, I)
        p_sd   <- matrix(NA_real_, N, I)
        
        if (simulate_alloc) {
                SPI_count <- matrix(0L, N, I)
        }
        
        for (n in 1:N) {
                for (i in 1:I) {
                        eta_s <- a0 + aN[, n] + bI[, i]
                        p_s   <- plogis(eta_s)
                        p_mean[n, i] <- mean(p_s)
                        p_q025[n, i] <- quantile(p_s, 0.025)
                        p_q975[n, i] <- quantile(p_s, 0.975)
                        p_sd[n, i]   <- sd(p_s)
                        if (simulate_alloc) {
                                SPI_count[n, i] <- sum(rbinom(S, size = 1, prob = p_s))
                        }
                }
        }
        
        hard <- (p_mean >= threshold) * 1L
        res$posterior_prob_A <- p_mean
        res$ci95 <- list(q025 = p_q025, q975 = p_q975)
        res$sd   <- p_sd
        res$hard_class <- hard
        if (simulate_alloc) {
                res$SPI_sim_mean <- SPI_count / S          # freq classified as A across simulated allocations
        }
        res$source <- "marginal"
        res
}
