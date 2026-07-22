library(nimble)
library(loo)

## Proper DIC (needs loglik at posterior mean parameters)
dic_from_ll <- function(loglik_mat, ll_at_mean) {
        mean_ll <- mean(rowSums(loglik_mat))
        llhat   <- sum(ll_at_mean)
        pD <- 2 * (llhat - mean_ll)
        DIC <- 2 * (-mean_ll*2) / 2  # placeholder; compute correctly:
        # DIC = 2*E[D] - D_hat = 2*(-2*mean_ll) - (-2*llhat) = -4*mean_ll + 2*llhat
        DIC <- -4*mean_ll + 2*llhat
        list(pD=pD, DIC=DIC, ll_mean = mean_ll, ll_hat = llhat)
}

# 3-category mixture density: x in {1,2,3}; returns scalar double
dMix3catTree <- nimbleFunction(
        run = function(x          = integer(0),   # must be named `x`
                       theta_dir  = double(0),
                       kappa_ext  = double(0),
                       bA1        = double(0),
                       bA2        = double(0),
                       bB1        = double(0),
                       bB2        = double(0),
                       p_gate     = double(0),
                       log        = integer(0, default = 0)) {
                returnType(double(0))
                
                # Tree A: Not-No -> Yes | Not-No
                pA1      <- 1/(1+exp(-(theta_dir - bA1)))   # not-No
                pA2c     <- 1/(1+exp(-(theta_dir - bA2)))   # Yes | not-No
                pA_no    <- 1 - pA1
                pA_yes   <- pA1 * pA2c
                pA_maybe <- 1 - pA_no - pA_yes
                
                # Tree B: Extreme -> Yes | Extreme
                pB1e     <- 1/(1+exp(-(kappa_ext - bB1)))   # Extreme
                pB2c     <- 1/(1+exp(-(theta_dir - bB2)))  # Yes | Extreme
                pB_maybe <- 1 - pB1e
                pB_yes   <- pB1e * pB2c
                pB_no    <- 1 - pB_maybe - pB_yes
                
                # Mixture probabilities for each category
                pmix1 <- p_gate * pA_no    + (1 - p_gate) * pB_no
                pmix2 <- p_gate * pA_maybe + (1 - p_gate) * pB_maybe
                pmix3 <- p_gate * pA_yes   + (1 - p_gate) * pB_yes
                
                # Select probability with control-flow (not an if-expression)
                prob <- 0.0
                if (x == 1L) {
                        prob <- pmix1
                } else if (x == 2L) {
                        prob <- pmix2
                } else {                      # assume x == 3L
                        prob <- pmix3
                }
                
                # numerical guard
                if (prob < 1e-12) prob <- 1e-12
                if (prob > 1 - 1e-12) prob <- 1 - 1e-12
                
                if (log) return(log(prob)) else return(prob)
        }
)

# (sampler is optional)
rMix3catTree <- nimbleFunction(
        run = function(n          = integer(0),
                       theta_dir  = double(0),
                       kappa_ext  = double(0),
                       bA1        = double(0),
                       bA2        = double(0),
                       bB1        = double(0),
                       bB2        = double(0),
                       p_gate     = double(0)) {
                returnType(integer(0))
                return(1L)
        }
)

# Re-register (restart R if you previously registered a broken version)
registerDistributions(list(
        dMix3catTree = list(
                BUGSdist = "dMix3catTree(theta_dir, kappa_ext, bA1, bA2, bB1, bB2, p_gate)",
                types    = c("value = integer(0)"),
                discrete = TRUE
        )
))


#get_vec
## ---------- Helpers: pulls ----------
as_mat <- function(samples) {
        if (inherits(samples, "mcmc.list")) {
                do.call(rbind, lapply(samples, as.matrix))
        } else {
                as.matrix(samples)
        }
}

get_vec <- function(M, base, expected_len = NULL) {
        if (is.null(M)) return(NULL)
        cn  <- colnames(M)
        pat <- paste0("^", base, "\\[(\\d+)\\]$")   # exact: base[<idx>]
        idx <- grep(pat, cn)
        if (!length(idx)) return(NULL)
        
        ord <- as.integer(sub(pat, "\\1", cn[idx]))
        # keep only indices ≤ expected_len when provided
        if (!is.null(expected_len)) {
                keep <- ord <= expected_len
                if (!all(keep)) {
                        warning(sprintf("get_vec('%s'): dropping %d columns with index > %d.",
                                        base, sum(!keep), expected_len))
                }
                idx <- idx[keep]
                ord <- ord[keep]
        }
        o <- order(ord)
        M[, idx[o], drop = FALSE]
}


# -------- Combine chains to a single draw-by-parameter matrix (safe for discrete nodes) -----
as_mat_disc <- function(samples) {
        if (inherits(samples, "mcmc.list")) {
                do.call(rbind, lapply(samples, as.matrix))
        } else {
                as.matrix(samples)
        }
}

# -------- Extract S x (N*I) matrix for a discrete array param like SPI[n,i] ------------------
# name: parameter base name, e.g., "SPI"
# N, I: dimensions of the array
# order: "i_n" (default) matches k stepping in your llmat loops: for (i in 1:I) for (n in 1:N)
#        "n_i" gives n outer, i inner (use if you prefer that convention)
get_vec_disc <- function(samples, name = "SPI", N, I, order = c("i_n", "n_i")) {
        order <- match.arg(order)
        M <- as_mat_disc(samples)
        cn <- colnames(M)
        if (is.null(cn)) stop("Samples object has no column names; cannot locate discrete nodes.")
        
        # Regex to capture "name[n,i]" columns
        pat <- paste0("^", name, "\\[([0-9]+),([0-9]+)\\]$")
        hits <- grep(pat, cn)
        if (length(hits) == 0L) stop("No columns found for ", name, " using pattern ", pat)
        
        # Parse n,i indices from column names
        # strcapture is robust and fast
        idx_df <- utils::strcapture(pat, cn[hits],
                                    proto = data.frame(n = integer(), i = integer()))
        if (nrow(idx_df) != length(hits)) stop("Parsing error for ", name, " column indices.")
        
        # Sanity check & grid completion
        full_grid <- expand.grid(n = 1:N, i = 1:I)
        # Which of the full (n,i) cells exist in samples?
        key <- function(df) paste(df$n, df$i, sep = ",")
        have <- key(idx_df)
        need <- key(full_grid)
        
        # Warn if any cells are missing; we will insert NA columns to keep shape/alignment
        if (!all(need %in% have)) {
                missing_keys <- setdiff(need, have)
                warning("get_vec_disc: missing draws for ", length(missing_keys), " cells of ", name,
                        ". Inserting NA columns for: ",
                        paste0(missing_keys, collapse = " "))
        }
        
        # Build the column order that matches the requested loop order
        if (order == "i_n") {
                # for (i in 1:I) for (n in 1:N)  -> i major, n minor
                ord_df <- full_grid[order(full_grid$i, full_grid$n), , drop = FALSE]
        } else {
                # for (n in 1:N) for (i in 1:I)  -> n major, i minor
                ord_df <- full_grid[order(full_grid$n, full_grid$i), , drop = FALSE]
        }
        
        # For each desired (n,i) in ord_df, pick the matching column or make NA
        lookup <- match(key(ord_df), have)
        # Initialize result with NA and then fill the present columns
        S <- nrow(M)
        out <- matrix(NA_real_, nrow = S, ncol = N * I)
        
        present <- which(!is.na(lookup))
        out[, present] <- M[, hits[lookup[present]], drop = FALSE]
        
        # Optional: name columns for debugging/consistency with llmat colnames
        colnames(out) <- paste0(name, "[", ord_df$n, ",", ord_df$i, "]")
        
        out
}

#summarize_params

## ---------- ICs: pD/DIC/WAIC/LOO ----------
# loglik_mat: S x K (rows = draws, cols = observations or your chosen unit)
# chain_id   : length S (which chain each draw came from). If NULL → assume 1 chain.
# r_eff      : optional precomputed loo::relative_eff(...) if you already have it
# auto_mm    : automatically run moment-matching if any Pareto-k > k_threshold
info_criteria <- function(loglik_mat, chain_id = NULL, r_eff = NULL,
                          auto_mm = TRUE, k_threshold = 0.7,
                          mm_post_draws = NULL) {
        
        loglik_mat <- as.matrix(loglik_mat)
        S <- nrow(loglik_mat); K <- ncol(loglik_mat)
        
        # helper
        if (!exists("log_mean_exp", mode = "function")) {
                log_mean_exp <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
        }
        
        # ---- WAIC ----
        lppd_i   <- apply(loglik_mat, 2, log_mean_exp)
        var_ll_i <- apply(loglik_mat, 2, stats::var)
        lppd       <- sum(lppd_i)
        p_waic     <- sum(var_ll_i)
        elpd_waic  <- sum(lppd_i - var_ll_i)
        WAIC       <- -2 * elpd_waic
        
        out <- list(lppd = lppd, p_waic = p_waic, elpd_waic = elpd_waic, WAIC = WAIC)
        
        # ---- LOO ----
        if (!requireNamespace("loo", quietly = TRUE)) {
                # crude IS-LOO fallback
                elpd_loo_i <- apply(loglik_mat, 2, function(v) { mx <- max(-v); -log(mean(exp(-(v) - mx))) - mx })
                out$LOO_elpd <- sum(elpd_loo_i)
                out$LOO_se   <- sqrt(K * stats::var(elpd_loo_i))
                out$p_loo    <- sum(lppd_i - elpd_loo_i)
                out$LOOIC    <- -2 * out$LOO_elpd
                out$Pareto_k <- NA
                out$n_k_gt_0_7 <- NA_integer_
                out$n_k_gt_1   <- NA_integer_
                return(out)
        }
        
        if (is.null(r_eff)) {
                if (is.null(chain_id)) chain_id <- rep(1L, S)
                r_eff <- loo::relative_eff(exp(loglik_mat), chain_id = chain_id)
        }
        loo_res <- loo::loo(loglik_mat, r_eff = r_eff)
        
        # ---- Optional: moment matching (version-proof) ----
        if (isTRUE(auto_mm)) {
                pk <- tryCatch(loo_res$diagnostics$pareto_k, error = function(...) NULL)
                if (!is.null(pk) && any(pk > k_threshold, na.rm = TRUE)) {
                        # Build args positionally: first arg = loo object (covers both 'loo' and 'x')
                        mm_args <- list(loo_res)
                        # Prefer post_draws if this version requires it; otherwise use log_lik (+ r_eff if accepted)
                        mm_fun <- get("loo_moment_match", envir = asNamespace("loo"))
                        mm_formals <- names(formals(mm_fun))
                        if ("post_draws" %in% mm_formals && !is.null(mm_post_draws)) {
                                mm_args$post_draws <- mm_post_draws
                        } else if ("log_lik" %in% mm_formals) {
                                mm_args$log_lik <- loglik_mat
                                if ("r_eff" %in% mm_formals) mm_args$r_eff <- r_eff
                        }
                        if ("k_threshold" %in% mm_formals) mm_args$k_threshold <- k_threshold
                        
                        mm_try <- tryCatch(do.call(loo::loo_moment_match, mm_args),
                                           error = function(e) {
                                                   warning(sprintf("Moment matching failed: %s ; falling back to original LOO.", e$message))
                                                   NULL
                                           })
                        if (!is.null(mm_try)) loo_res <- mm_try
                }
        }
        
        # ---- Collect ----
        out$LOO_elpd <- loo_res$estimates["elpd_loo", "Estimate"]
        out$LOO_se   <- loo_res$estimates["elpd_loo", "SE"]
        out$p_loo    <- loo_res$estimates["p_loo",  "Estimate"]
        out$LOOIC    <- -2 * out$LOO_elpd
        
        pk <- tryCatch(loo_res$diagnostics$pareto_k, error = function(...) NULL)
        out$Pareto_k    <- pk
        out$n_k_gt_0_7  <- if (is.null(pk)) NA_integer_ else sum(pk > 0.7, na.rm = TRUE)
        out$n_k_gt_1    <- if (is.null(pk)) NA_integer_ else sum(pk > 1,   na.rm = TRUE)
        out$max_k       <- if (is.null(pk)) NA_real_     else max(pk, na.rm = TRUE)
        out$r_eff       <- r_eff
        out$loo_obj     <- loo_res
        
        out
}

## ---------- (optional) DIC using the helper above ----------
# mean_D  = -2 * E_s[ sum_k loglik_{s,k} ]
# D_hat   = -2 * sum_k loglik at posterior-mean params
# pD      = mean_D - D_hat
# DIC     = D_hat + 2*pD
dic_from_loglik_and_mean <- function(loglik_mat, ll_at_mean_vec) {
        mean_D <- -2 * mean(rowSums(loglik_mat))
        D_hat  <- -2 * sum(ll_at_mean_vec)
        pD     <- mean_D - D_hat
        DIC    <- D_hat + 2 * pD
        list(mean_D = mean_D, D_hat = D_hat, pD = pD, DIC = DIC)
}

summarize_params <- function(samples, N, I,
                             model = c("A","B","G1","G2","G3","G4","G5","G6","G7"),
                             D = NULL,
                             return_G4_wz_draws = TRUE,
                             reshape_wz = TRUE) {
        model <- match.arg(model)
        M <- as_mat(samples)
        qs <- function(v) stats::quantile(v, c(.025, .5, .975), na.rm = TRUE)
        
        # --- helpers for Rhat ---
        .has_coda <- requireNamespace("coda", quietly = TRUE)
        
        .rhat_vec <- function(samples, base, dim) {
                if (!.has_coda || !inherits(samples, "mcmc.list") || length(samples) < 2) return(NULL)
                sel <- paste0(base, "[", seq_len(dim), "]")
                ml <- lapply(samples, function(ch) {
                        mat <- as.matrix(ch)
                        cols <- sel[sel %in% colnames(mat)]
                        if (!length(cols)) return(NULL)
                        coda::as.mcmc(mat[, cols, drop = FALSE])
                })
                ml <- Filter(Negate(is.null), ml)
                if (!length(ml)) return(NULL)
                gd <- coda::gelman.diag(coda::as.mcmc.list(ml), autoburnin = FALSE, multivariate = FALSE)
                R <- gd$psrf[, "Point est."]
                out <- rep(NA_real_, dim); names(out) <- sel
                out[names(R)] <- R
                as.numeric(out)
        }
        
        .rhat_scalar <- function(samples, name) {
                if (!.has_coda || !inherits(samples, "mcmc.list") || length(samples) < 2) return(NA_real_)
                ml <- lapply(samples, function(ch) {
                        mat <- as.matrix(ch)
                        if (!(name %in% colnames(mat))) return(NULL)
                        coda::as.mcmc(mat[, name, drop = FALSE])
                })
                ml <- Filter(Negate(is.null), ml)
                if (!length(ml)) return(NA_real_)
                gd <- coda::gelman.diag(coda::as.mcmc.list(ml), autoburnin = FALSE, multivariate = FALSE)
                as.numeric(gd$psrf[, "Point est."])
        }
        
        # Rhat for a 2D discrete node like SPI[n,i]; returns vector in "i_n" order (for i=1..I, n=1..N)
        .rhat_SPI <- function(samples, N, I, name = "SPI", order = c("i_n","n_i")) {
                order <- match.arg(order)
                if (!.has_coda || !inherits(samples, "mcmc.list") || length(samples) < 2) return(NULL)
                pat <- paste0("^", name, "\\[([0-9]+),([0-9]+)\\]$")
                # Common SPI columns across chains (ensure alignment)
                spi_cols_list <- lapply(samples, function(ch) grep(pat, colnames(as.matrix(ch)), value = TRUE))
                common <- Reduce(intersect, spi_cols_list)
                if (!length(common)) return(NULL)
                
                # Build desired full order
                full <- expand.grid(n = 1:N, i = 1:I)
                if (order == "i_n") ord_df <- full[order(full$i, full$n), ] else ord_df <- full[order(full$n, full$i), ]
                desired_names <- paste0(name, "[", ord_df$n, ",", ord_df$i, "]")
                
                # Restrict to common columns & order them
                keep <- desired_names[desired_names %in% common]
                if (!length(keep)) return(NULL)
                
                ml <- lapply(samples, function(ch) {
                        mat <- as.matrix(ch)
                        cols <- keep[keep %in% colnames(mat)]
                        coda::as.mcmc(mat[, cols, drop = FALSE])
                })
                ml <- Filter(Negate(is.null), ml)
                if (!length(ml)) return(NULL)
                
                gd <- coda::gelman.diag(coda::as.mcmc.list(ml), autoburnin = FALSE, multivariate = FALSE)
                R <- gd$psrf[, "Point est."]  # named by keep
                
                # Place into full N*I vector (i_n order) with NAs for missing
                out <- rep(NA_real_, N * I); names(out) <- desired_names
                out[names(R)] <- R
                as.numeric(out)
        }
        
        out <- list()
        
        # -------- Persons --------
        th <- get_vec(M, "theta_dir", N)
        if (!is.null(th)) {
                df <- data.frame(id = 1:N, mean = colMeans(th), t(apply(th, 2, qs)))
                Rh <- .rhat_vec(samples, "theta_dir", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$theta_dir <- df
        }
        
        # include kappa_ext for all models that are not pure A (B, G1..G7)
        if (model != "A") {
                kp <- get_vec(M, "kappa_ext", N)
                if (!is.null(kp)) {
                        df <- data.frame(id = 1:N, mean = colMeans(kp), t(apply(kp, 2, qs)))
                        Rh <- .rhat_vec(samples, "kappa_ext", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$kappa_ext <- df
                }
        }
        
        # G6 person indicators + alpha_person
        if (model == "G6") {
                sp <- get_vec(M, "SP", N)
                if (!is.null(sp)) {
                        df <- data.frame(id = 1:N, mean = colMeans(sp), t(apply(sp, 2, qs)))
                        Rh <- .rhat_vec(samples, "SP", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$SP <- df
                }
                if (!is.null(colnames(M)) && "alpha_person" %in% colnames(M)) {
                        ap <- M[, "alpha_person"]
                        df <- data.frame(mean = mean(ap), t(qs(ap)))
                        Rh <- .rhat_scalar(samples, "alpha_person"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$alpha_person <- df
                }
        }
        
        # G3: r[n], a[i], delta[i]
        if (model == "G3") {
                r <- get_vec(M, "r", N)
                if (!is.null(r)) {
                        df <- data.frame(id = 1:N, mean = colMeans(r), t(apply(r, 2, qs)))
                        Rh <- .rhat_vec(samples, "r", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$r <- df
                }
                a <- get_vec(M, "a", I)
                if (!is.null(a)) {
                        df <- data.frame(item = 1:I, mean = colMeans(a), t(apply(a, 2, qs)))
                        Rh <- .rhat_vec(samples, "a", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$a <- df
                }
                dl <- get_vec(M, "delta", I)
                if (!is.null(dl)) {
                        df <- data.frame(item = 1:I, mean = colMeans(dl), t(apply(dl, 2, qs)))
                        Rh <- .rhat_vec(samples, "delta", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$delta <- df
                }
        }
        
        # G4: u[n], v[i]; optional z/w draws
        if (model == "G4") {
                u <- get_vec(M, "u", N)
                if (!is.null(u)) {
                        df <- data.frame(id = 1:N, mean = colMeans(u), t(apply(u, 2, qs)))
                        Rh <- .rhat_vec(samples, "u", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$u <- df
                }
                v <- get_vec(M, "v", I)
                if (!is.null(v)) {
                        df <- data.frame(item = 1:I, mean = colMeans(v), t(apply(v, 2, qs)))
                        Rh <- .rhat_vec(samples, "v", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$v <- df
                }
                
                if (isTRUE(return_G4_wz_draws)) {
                        if (is.null(D)) stop("For model G4, please provide D (latent space dimension).")
                        z_mat <- get_vec(M, "z", N * D)  # S x (N*D)
                        w_mat <- get_vec(M, "w", I * D)  # S x (I*D)
                        out$draws <- out$draws %||% list()
                        out$draws$z_mat <- z_mat
                        out$draws$w_mat <- w_mat
                        
                        if (isTRUE(reshape_wz) && !is.null(z_mat) && !is.null(w_mat)) {
                                reshape_SxKD <- function(X, K, D) {
                                        if (is.null(X)) return(NULL)
                                        S <- nrow(X)
                                        arr <- array(NA_real_, c(S, K, D))
                                        for (d in 1:D) arr[ , , d] <- X[ , ((d - 1) * K + 1):(d * K), drop = FALSE]
                                        arr
                                }
                                out$draws$z_arr <- reshape_SxKD(z_mat, N, D)  # S x N x D
                                out$draws$w_arr <- reshape_SxKD(w_mat, I, D)  # S x I x D
                        }
                }
        }
        
        # G7: u[n], v[i] and scalars mu, sigma_u, sigma_v, phi
        if (model == "G7") {
                u <- get_vec(M, "u", N)
                if (!is.null(u)) {
                        df <- data.frame(id = 1:N, mean = colMeans(u), t(apply(u, 2, qs)))
                        Rh <- .rhat_vec(samples, "u", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$u <- df
                }
                v <- get_vec(M, "v", I)
                if (!is.null(v)) {
                        df <- data.frame(item = 1:I, mean = colMeans(v), t(apply(v, 2, qs)))
                        Rh <- .rhat_vec(samples, "v", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$v <- df
                }
                if ("mu" %in% colnames(M)) {
                        x <- M[, "mu"]; df <- data.frame(mean = mean(x), t(qs(x)))
                        Rh <- .rhat_scalar(samples, "mu"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$mu <- df
                }
                if ("sigma_u" %in% colnames(M)) {
                        x <- M[, "sigma_u"]; df <- data.frame(mean = mean(x), t(qs(x)))
                        Rh <- .rhat_scalar(samples, "sigma_u"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$sigma_u <- df
                }
                if ("sigma_v" %in% colnames(M)) {
                        x <- M[, "sigma_v"]; df <- data.frame(mean = mean(x), t(qs(x)))
                        Rh <- .rhat_scalar(samples, "sigma_v"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$sigma_v <- df
                }
                if ("phi" %in% colnames(M)) {
                        x <- M[, "phi"]; df <- data.frame(mean = mean(x), t(qs(x)))
                        Rh <- .rhat_scalar(samples, "phi"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$phi <- df
                }
        }
        
        # -------- Items (trees) --------
        if (model == "A") {
                b1 <- get_vec(M, "bA1", I); b2 <- get_vec(M, "bA2", I)
                if (!is.null(b1)) { df <- data.frame(item = 1:I, mean = colMeans(b1), t(apply(b1, 2, qs)))
                Rh <- .rhat_vec(samples, "bA1", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bA1 <- df }
                if (!is.null(b2)) { df <- data.frame(item = 1:I, mean = colMeans(b2), t(apply(b2, 2, qs)))
                Rh <- .rhat_vec(samples, "bA2", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bA2 <- df }
        } else if (model == "B") {
                b1 <- get_vec(M, "bB1", I); b2 <- get_vec(M, "bB2", I)
                if (!is.null(b1)) { df <- data.frame(item = 1:I, mean = colMeans(b1), t(apply(b1, 2, qs)))
                Rh <- .rhat_vec(samples, "bB1", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bB1 <- df }
                if (!is.null(b2)) { df <- data.frame(item = 1:I, mean = colMeans(b2), t(apply(b2, 2, qs)))
                Rh <- .rhat_vec(samples, "bB2", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bB2 <- df }
        } else {  # mixtures/gates (incl. updated G1, G2, G5, G6, G7)
                bA1 <- get_vec(M, "bA1", I); bA2 <- get_vec(M, "bA2", I)
                bB1 <- get_vec(M, "bB1", I); bB2 <- get_vec(M, "bB2", I)
                if (!is.null(bA1)) { df <- data.frame(item = 1:I, mean = colMeans(bA1), t(apply(bA1, 2, qs)))
                Rh <- .rhat_vec(samples, "bA1", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bA1 <- df }
                if (!is.null(bA2)) { df <- data.frame(item = 1:I, mean = colMeans(bA2), t(apply(bA2, 2, qs)))
                Rh <- .rhat_vec(samples, "bA2", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bA2 <- df }
                if (!is.null(bB1)) { df <- data.frame(item = 1:I, mean = colMeans(bB1), t(apply(bB1, 2, qs)))
                Rh <- .rhat_vec(samples, "bB1", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bB1 <- df }
                if (!is.null(bB2)) { df <- data.frame(item = 1:I, mean = colMeans(bB2), t(apply(bB2, 2, qs)))
                Rh <- .rhat_vec(samples, "bB2", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$bB2 <- df }
                
                # G5 item indicators
                if (model == "G5") {
                        si <- get_vec(M, "SI", I)
                        if (!is.null(si)) {
                                df <- data.frame(item = 1:I, mean = colMeans(si), t(apply(si, 2, qs)))
                                Rh <- .rhat_vec(samples, "SI", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                                out$SI <- df
                        }
                }
        }
        
        # -------- Scalars / hyperparameters for gates --------
        if (model %in% c("G1","G2","G3","G4") && "alpha0" %in% colnames(M)) {
                a0 <- M[, "alpha0"]; df <- data.frame(mean = mean(a0), t(qs(a0)))
                Rh <- .rhat_scalar(samples, "alpha0"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                out$alpha0 <- df
        }
        if (model == "G2") {
                ag <- get_vec(M, "alpha_gate", N)
                if (!is.null(ag)) {
                        df <- data.frame(id = 1:N, mean = colMeans(ag), t(apply(ag, 2, qs)))
                        Rh <- .rhat_vec(samples, "alpha_gate", N); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$alpha_gate <- df
                }
                bg <- get_vec(M, "beta_gate", I)
                if (!is.null(bg)) {
                        df <- data.frame(item = 1:I, mean = colMeans(bg), t(apply(bg, 2, qs)))
                        Rh <- .rhat_vec(samples, "beta_gate", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$beta_gate <- df
                }
                if ("sigma_alpha" %in% colnames(M)) {
                        x <- M[, "sigma_alpha"]; df <- data.frame(mean = mean(x), t(qs(x)))
                        Rh <- .rhat_scalar(samples, "sigma_alpha"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$sigma_alpha <- df
                }
                if ("sigma_beta" %in% colnames(M)) {
                        x <- M[, "sigma_beta"]; df <- data.frame(mean = mean(x), t(qs(x)))
                        Rh <- .rhat_scalar(samples, "sigma_beta"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$sigma_beta <- df
                }
        }
        
        # d2 hierarchy shared by mixture/gate models
        if (model %in% c("G1","G2","G3","G4","G5","G6","G7")) {
                if ("mu_d2" %in% colnames(M)) {
                        md2 <- M[, "mu_d2"]; df <- data.frame(mean = mean(md2), t(qs(md2)))
                        Rh <- .rhat_scalar(samples, "mu_d2"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$mu_d2 <- df
                }
                if ("sigma_d2" %in% colnames(M)) {
                        sd2 <- M[, "sigma_d2"]; df <- data.frame(mean = mean(sd2), t(qs(sd2)))
                        Rh <- .rhat_scalar(samples, "sigma_d2"); if (is.finite(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$sigma_d2 <- df
                }
                d2 <- get_vec(M, "d2", I)
                if (!is.null(d2)) {
                        df <- data.frame(item = 1:I, mean = colMeans(d2), t(apply(d2, 2, qs)))
                        Rh <- .rhat_vec(samples, "d2", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$d2 <- df
                }
                tau_vec <- get_vec(M, "tau", I)
                if (!is.null(d2)) {
                        df <- data.frame(item = 1:I, mean = colMeans(tau_vec), t(apply(tau_vec, 2, qs)))
                        Rh <- .rhat_vec(samples, "tau", I); if (!is.null(Rh)) { df$Rhat <- Rh; df$Rhat2 <- Rh^2 }
                        out$tau <- df
                }
        }
        
        # -------- SPI summaries (if present) --------
        if (any(grepl("^SPI\\[", colnames(M)))) {
                # S x (N*I) matrix, columns ordered i-major, n-minor to match your llmat loops
                SPI_mat <- tryCatch(get_vec_disc(M, "SPI", N, I, order = "i_n"),
                                    error = function(e) NULL)
                if (!is.null(SPI_mat)) {
                        S <- nrow(SPI_mat)
                        # Cell-level posterior summaries
                        cell_mean <- colMeans(SPI_mat, na.rm = TRUE)
                        cell_qs   <- t(apply(SPI_mat, 2, qs))
                        MAP       <- as.integer(cell_mean >= 0.5)
                        
                        # Optional R-hat per cell (NA if not available)
                        Rhat_spi <- .rhat_SPI(samples, N, I, name = "SPI", order = "i_n")
                        Rhat2_spi <- if (is.null(Rhat_spi)) NULL else Rhat_spi^2
                        
                        # (n,i) indices in i-major order
                        idx_df <- expand.grid(n = 1:N, i = 1:I)
                        idx_df <- idx_df[order(idx_df$i, idx_df$n), , drop = FALSE]
                        
                        SPI_cell <- data.frame(
                                n = idx_df$n, i = idx_df$i,
                                mean = cell_mean,
                                `2.5%` = cell_qs[, 1],
                                `50%`  = cell_qs[, 2],
                                `97.5%`= cell_qs[, 3],
                                MAP = MAP
                        )
                        if (!is.null(Rhat_spi)) {
                                SPI_cell$Rhat  <- Rhat_spi
                                SPI_cell$Rhat2 <- Rhat2_spi
                        }
                        out$SPI_cell <- SPI_cell
                        
                        # Per-draw reshape: S x N x I array
                        SPI_arr <- array(SPI_mat, dim = c(S, N, I))  # i-major order preserved
                        
                        # Per-person (avg over items) and per-item (avg over persons) per draw
                        SPI_person_draw <- apply(SPI_arr, c(1, 2), mean, na.rm = TRUE)  # S x N
                        SPI_item_draw   <- apply(SPI_arr, c(1, 3), mean, na.rm = TRUE)  # S x I
                        SPI_overall_draw<- rowMeans(SPI_mat, na.rm = TRUE)              # length S
                        
                        out$SPI_overall <- data.frame(mean = mean(SPI_overall_draw),
                                                      t(qs(SPI_overall_draw)))
                        
                        dfp <- data.frame(id = 1:N,
                                          mean = colMeans(SPI_person_draw),
                                          t(apply(SPI_person_draw, 2, qs)))
                        # (No simple per-person SPI Rhat here; it’s an aggregated derived node.)
                        out$SPI_person <- dfp
                        
                        dfi <- data.frame(item = 1:I,
                                          mean = colMeans(SPI_item_draw),
                                          t(apply(SPI_item_draw, 2, qs)))
                        out$SPI_item <- dfi
                }
        }
        
        out
}


report_gof <- function(res){
        return(
                c(
                        LOO_elpd = res$IC$loo_obj$estimates[1,1],
                        LOO_elpd_se = res$IC$loo_obj$estimates[1,2],
                        LOOIC = res$IC$loo_obj$estimates[3,1],
                        pLOO = res$IC$loo_obj$estimates[2,1],
                        WAIC_elpd = res$IC$elpd_waic,
                        pWAIC = res$IC$p_waic,
                        WAIC = res$IC$WAIC,
                        pD = res$IC$pD,
                        DIC = res$IC$DIC,
                        Pareto_gt_0_7 = res$IC$n_k_gt_0_7,
                        Pareto_gt_1 = res$IC$n_k_gt_1,
                        max_k = res$IC$max_k   
                )
        )
}