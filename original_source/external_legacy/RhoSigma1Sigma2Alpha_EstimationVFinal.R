library(yuima)

env_num <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) return(default)
  out <- suppressWarnings(as.numeric(val))
  if (is.finite(out)) out else default
}

env_bool <- function(name, default = FALSE) {
  val <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(val)) return(default)
  val %in% c("1", "true", "yes", "y")
}

#--------------------- Simulation parameters
Terminal   <- env_num("TERMINAL", 1)
n_steps    <- as.integer(env_num("N_STEPS", 10000000))
alpha_true <- env_num("ALPHA_TRUE", 1.25)
sigma2_true<- env_num("SIGMA2_TRUE", 1)
sigma1_true<- env_num("SIGMA1_TRUE", 1.2)
rho_true   <- env_num("RHO_TRUE", 1.45)
init_v_stationary <- env_bool("INIT_V_STATIONARY", FALSE)
seed_val   <- Sys.getenv("SEED", unset = "")
if (nzchar(seed_val)) set.seed(as.integer(seed_val))
Uprime <- function(x) 16*x^3 - 16*x

V0 <- if (init_v_stationary) rt(1, df = 3) / sqrt(3) else 0

#------------ Yuima model
mod2D <- setModel(
  drift = c("-(16*x^3 - 16*x) - sigma1*(rho^2)*v",
            "-rho^2*v"),
  diffusion = matrix(c("sigma1*rho*sqrt(1+v^2)",
                       "rho*sqrt(1+v^2)"),
                     nrow = 2, byrow = TRUE),
  jump.coeff = matrix(c("sigma2","0"), nrow = 2),
  measure.type = "code",
  measure = list(df = "rstable(z, alpha, 0, 1, 0)"),
  state.variable = c("x","v"),
  solve.variable = c("x","v"),
  xinit = c("x"=0,"v"=V0)
)

samp <- setSampling(Terminal = Terminal, n = n_steps)
yu   <- setYuima(model = mod2D, sampling = samp)
sim  <- simulate(yu, true.par = list(
  sigma1 = sigma1_true, rho = rho_true,
  sigma2 = sigma2_true, alpha = alpha_true
))
z      <- yuima::get.zoo.data(sim)
time   <- as.numeric(index(z[[1]]))
X      <- as.numeric(z[[1]])
V      <- as.numeric(z[[2]])

dt <- diff(time)
dX <- diff(X)
R  <- dX + Uprime(X[-length(X)]) * dt
A  <- sort(abs(R[is.finite(R) & R != 0]), decreasing = TRUE)
m  <- length(A)

# --------------------- Hill estimator
hill_alpha <- function(y_sorted_desc) {
  m <- length(y_sorted_desc)
  if (m < 5) stop("Not enough tail points for Hill.")
  logY <- log(y_sorted_desc)
  csum_logY <- cumsum(logY[1:(m-1)])
  logY_kp1  <- logY[2:m]
  k_vec     <- 1:(m-1)
  gamma_k <- (csum_logY - k_vec * logY_kp1) / k_vec
  1 / gamma_k
}

runmed1 <- function(x, k = 5L) {
  n <- length(x); if (n == 0L) return(x)
  k <- max(1L, as.integer(k)); if (k %% 2L == 0L) k <- k + 1L
  half <- k %/% 2L
  xm <- x
  for (i in seq_len(n)) {
    lo <- max(1L, i - half); hi <- min(n, i + half)
    xm[i] <- median(x[lo:hi], na.rm = TRUE)
  }
  xm
}


choose_k_via_plateau <- function(hill_vals, alpha_cap = 1.98) {
  m1 <- length(hill_vals)

  k_min <- 20L
  k_max <- min(as.integer(m1 / 100L), 20000L, m1 - 10L)

  if (!is.finite(k_max) || k_max < k_min) {
    fallback_k <- min(max(2L, as.integer(floor(m1 / 10L))), m1)
    warning(sprintf(
      "Hill plateau selection is underpowered with only %d tail candidates; using fallback k=%d for a mechanical run only.",
      m1, fallback_k
    ))
    return(list(k = fallback_k, score = NA_real_, range = c(1L, m1),
                small_sample = TRUE))
  }

  # cumulative sums for O(1) segment variance
  cs  <- c(0, cumsum(hill_vals))
  cs2 <- c(0, cumsum(hill_vals^2))
  
  # log-spaced candidate k-grid
  k_grid <- unique(as.integer(round(exp(seq(log(k_min), log(k_max), length.out = 60L)))))
  k_grid <- k_grid[k_grid >= k_min & k_grid <= k_max]
  
  best_k <- NA_integer_
  best_score <- Inf
  
  for (k in k_grid) {
    a <- as.integer(k / 2L)
    if (a < 1L) a <- 1L
    len <- k - a + 1L
    sum_seg   <- cs[k + 1L]  - cs[a]
    sumsq_seg <- cs2[k + 1L] - cs2[a]
    mean_seg  <- sum_seg / len
    var_seg   <- max(sumsq_seg / len - mean_seg * mean_seg, 0)
    sd_seg    <- sqrt(var_seg)
    
    if (!is.finite(sd_seg)) next
    if (hill_vals[k] > alpha_cap) next
    
    if (sd_seg < best_score) {
      best_score <- sd_seg
      best_k     <- k
    }
  }
  
  if (is.na(best_k)) best_k <- k_min
  list(k = best_k, score = best_score, range = range(k_grid))
}

# ----- Alpha estimation
hill_vec_raw <- hill_alpha(A)
hill_vec     <- runmed1(hill_vec_raw, k = 5L)
sel      <- choose_k_via_plateau(hill_vec)
k_hat    <- sel$k
alpha_hat<- hill_vec[k_hat]

q_idx <- k_hat + 1L
q_val <- if (q_idx <= m) A[q_idx] else A[m]
q_prob<- (q_idx)/(m+1)
cat(sprintf("\nEstimated alpha  = %.4f\n", alpha_hat))
#cat(sprintf("Tail threshold |R| at k: %.4e  (empirical exceed prob = %.4g)\n", q_val, q_prob))
cat(sprintf("True alpha = %.4f\n", alpha_true))

# ===== Improved, still-simple sigma2 estimator (same formula, less noise) =====
estimate_sigma2_tail_robust <- function(alpha_hat, A_sorted_desc, k_hat, dt_vec,
                                        w_half = 20L, k_bounds = NULL) {
  m <- length(A_sorted_desc)

  dtrob <- median(dt_vec)
  c_alpha <- gamma(alpha_hat) * sin(pi * alpha_hat / 2) / pi
  
  if (is.null(k_bounds)) {
    k_lo <- max(2L, k_hat - w_half)
    k_hi <- min(m-1L, k_hat + w_half)
  } else {
    k_lo <- max(2L, k_bounds[1])
    k_hi <- min(m-1L, k_bounds[2])
  }
  
  k_grid <- k_lo:k_hi
  s2h <- numeric(length(k_grid))
  idx <- 0L
  for (kj in k_grid) {
    idx <- idx + 1L
    u_j   <- A_sorted_desc[kj]
    p_j   <- (kj - 0.5) / m
    num   <- p_j * u_j^alpha_hat
    den   <- 2 * c_alpha * dtrob
    s2a   <- num / den
    s2h[idx] <- if (is.finite(s2a) && s2a > 0) s2a^(1 / alpha_hat) else NA_real_
  }
  sigma2_hat_med <- median(s2h[is.finite(s2h)], na.rm = TRUE)
  
  list(
    sigma2_hat = sigma2_hat_med,
    details = list(k_hat = k_hat, k_lo = k_lo, k_hi = k_hi,
                   dtrob = dtrob, c_alpha = c_alpha, window_size = length(k_grid))
  )
}

sig2_out <- estimate_sigma2_tail_robust(alpha_hat, A, k_hat, dt)
cat(sprintf("Estimated sigma2  = %.4f\n", sig2_out$sigma2_hat))
#cat(sprintf("  window k in [%d,%d] (size=%d) | c(alpha)=%.5f | median(dt)=%.3e\n",
            # sig2_out$details$k_lo,
            # sig2_out$details$k_hi,
            # sig2_out$details$window_size,
            # sig2_out$details$c_alpha,
            # sig2_out$details$dtrob))
cat(sprintf("True sigma2 = %.4f\n", sigma2_true))

#---------Observed-data reconstruction of V_t, then Sigma1 and rho
make_jump_filter <- function(R, dt, alpha_hat, sigma2_hat, k_hat, mode = "hard",
                             tail_prob_floor = 0,
                             tail_prob_cap = 0.01,
                             threshold_mode = "max") {
  mode <- match.arg(mode, c("hard", "soft"))
  threshold_mode <- match.arg(threshold_mode,
                              c("max", "min", "stable", "empirical"))

  n <- length(R)
  if (length(dt) != n) stop("R and dt must have the same length.")
  if (!is.finite(alpha_hat) || alpha_hat <= 0 || alpha_hat >= 2) {
    stop("alpha_hat must be in the stable range (0, 2).")
  }
  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) stop("sigma2_hat must be positive.")

  c_alpha <- gamma(alpha_hat) * sin(pi * alpha_hat / 2) / pi
  p_tail <- (k_hat + 1) / (n + 1)
  p_tail <- max(p_tail, tail_prob_floor, 1 / (n + 1))
  p_tail <- min(p_tail, tail_prob_cap)

  dt_pos <- ifelse(is.finite(dt) & dt > 0, dt, median(dt[is.finite(dt) & dt > 0]))
  absR <- abs(R)
  jump_threshold_stable <- (2 * c_alpha * sigma2_hat^alpha_hat * dt_pos / p_tail)^(1 / alpha_hat)
  empirical_threshold <- as.numeric(quantile(absR[is.finite(absR)],
                                             probs = 1 - p_tail,
                                             na.rm = TRUE,
                                             names = FALSE))
  jump_threshold <- switch(
    threshold_mode,
    max = pmax(jump_threshold_stable, empirical_threshold),
    min = pmin(jump_threshold_stable, empirical_threshold),
    stable = jump_threshold_stable,
    empirical = rep(empirical_threshold, n)
  )
  jump_threshold[!is.finite(jump_threshold) | jump_threshold <= 0] <- Inf

  is_jump <- !is.finite(R) | absR > jump_threshold

  if (mode == "hard") {
    weights <- as.numeric(!is_jump)
  } else {
    weights <- rep(1, n)
    large <- is.finite(absR) & absR > jump_threshold
    weights[large] <- pmax((jump_threshold[large] / absR[large])^alpha_hat, 0)
    weights[!is.finite(weights)] <- 0
  }

  list(
    weights = weights,
    is_jump = is_jump,
    threshold = jump_threshold,
    stable_threshold = jump_threshold_stable,
    empirical_threshold = empirical_threshold,
    p_tail = p_tail,
    c_alpha = c_alpha,
    mode = mode,
    threshold_mode = threshold_mode
  )
}

estimate_sigma1_rho_from_vhat <- function(V_hat, dt, weights,
                                          rho_estimator = "qv") {
  rho_estimator <- match.arg(rho_estimator, c("qv", "contrast"))
  dV_hat <- diff(V_hat)
  V_left_hat <- V_hat[-length(V_hat)]
  good <- is.finite(dV_hat) & is.finite(V_left_hat) &
    is.finite(dt) & dt > 0 & is.finite(weights) & weights > 0

  if (!any(good)) {
    return(list(rho_hat = NA_real_, rho2_hat = NA_real_,
                rho_qv = NA_real_, rho2_qv = NA_real_,
                rho_contrast = NA_real_, lambda_contrast = NA_real_,
                rho_robust = NA_real_, rho2_robust = NA_real_,
                dV_hat = dV_hat))
  }

  num_rho <- sum(weights[good] * dV_hat[good]^2)
  den_rho <- sum(weights[good] * (1 + V_left_hat[good]^2) * dt[good])
  rho2_qv <- if (is.finite(den_rho) && den_rho > 0) num_rho / den_rho else NA_real_
  rho_qv <- if (is.finite(rho2_qv)) sqrt(max(rho2_qv, 0)) else NA_real_

  contrast_num <- sum(weights[good] * dV_hat[good]^2 /
                        ((1 + V_left_hat[good]^2) * dt[good]))
  contrast_den <- sum(weights[good] * V_left_hat[good]^2 * dt[good] /
                        (1 + V_left_hat[good]^2))
  lambda_contrast <- if (is.finite(contrast_num) && is.finite(contrast_den) &&
                         contrast_num >= 0 && contrast_den > 0) {
    sqrt(contrast_num / contrast_den)
  } else {
    NA_real_
  }
  rho_contrast <- if (is.finite(lambda_contrast)) sqrt(max(lambda_contrast, 0)) else NA_real_

  rho_hat <- if (rho_estimator == "contrast") rho_contrast else rho_qv
  rho2_hat <- rho_hat^2
  rho_unit <- dV_hat[good]^2 / ((1 + V_left_hat[good]^2) * dt[good])
  rho2_robust <- median(rho_unit[is.finite(rho_unit)], na.rm = TRUE) /
    qchisq(0.5, df = 1)
  rho_robust <- if (is.finite(rho2_robust)) sqrt(max(rho2_robust, 0)) else NA_real_

  list(rho_hat = rho_hat, rho2_hat = rho2_hat,
       rho_qv = rho_qv, rho2_qv = rho2_qv,
       rho_contrast = rho_contrast, lambda_contrast = lambda_contrast,
       rho_robust = rho_robust, rho2_robust = rho2_robust,
       dV_hat = dV_hat)
}

estimate_sigma1_from_diffusion_qv <- function(Y_centered, dt, weights,
                                              trim_q = 0.95,
                                              n_bins = 25L) {
  dY <- diff(Y_centered)
  Y_left <- Y_centered[-length(Y_centered)]
  good <- is.finite(dY) & is.finite(Y_left) &
    is.finite(dt) & dt > 0 & is.finite(weights) & weights > 0

  if (sum(good) < 50L) {
    return(list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                slope = NA_real_, intercept = NA_real_))
  }

  x <- Y_left[good]^2
  z <- dY[good]^2 / (qchisq(0.5, df = 1) * dt[good])
  finite <- is.finite(x) & is.finite(z) & z >= 0
  x <- x[finite]
  z <- z[finite]

  if (length(z) < 50L || length(unique(x)) < 2L) {
    return(list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                slope = NA_real_, intercept = NA_real_))
  }

  z_cap <- as.numeric(quantile(z, probs = trim_q, na.rm = TRUE, names = FALSE))
  keep <- is.finite(z_cap) & z <= z_cap
  x <- x[keep]
  z <- z[keep]

  if (length(z) < 50L || length(unique(x)) < 2L) {
    return(list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                slope = NA_real_, intercept = NA_real_))
  }

  breaks <- unique(as.numeric(quantile(
    x,
    probs = seq(0, 1, length.out = as.integer(n_bins) + 1L),
    na.rm = TRUE,
    names = FALSE
  )))

  if (length(breaks) >= 3L) {
    bins <- cut(x, breaks = breaks, include.lowest = TRUE)
    bx <- as.numeric(tapply(x, bins, median, na.rm = TRUE))
    bz <- as.numeric(tapply(z, bins, median, na.rm = TRUE))
    bw <- as.numeric(tapply(z, bins, length))
    ok <- is.finite(bx) & is.finite(bz) & is.finite(bw) & bw > 0
    bx <- bx[ok]
    bz <- bz[ok]
    bw <- bw[ok]
  } else {
    bx <- x
    bz <- z
    bw <- rep(1, length(z))
  }

  slope <- intercept <- NA_real_
  if (length(bz) >= 3L && length(unique(bx)) >= 2L) {
    fit <- tryCatch(lm(bz ~ bx, weights = bw), error = function(e) NULL)
    if (!is.null(fit)) {
      co <- coef(fit)
      intercept <- unname(co[1])
      slope <- unname(co[2])
    }
  }

  if (!is.finite(slope) || slope <= 0 ||
      !is.finite(intercept) || intercept <= 0) {
    return(list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                slope = slope, intercept = intercept))
  }

  list(
    sigma1_hat = sqrt(intercept / slope),
    rho_hat = sqrt(slope),
    slope = slope,
    intercept = intercept
  )
}

estimate_center_from_ou_drift <- function(Y_raw, dt, weights) {
  dY <- diff(Y_raw)
  Y_left <- Y_raw[-length(Y_raw)]
  good <- is.finite(dY) & is.finite(Y_left) &
    is.finite(dt) & dt > 0 & is.finite(weights) & weights > 0

  if (sum(good) < 100L || sd(Y_left[good]) <= 0) {
    return(list(center = NA_real_, lambda = NA_real_))
  }

  x <- cbind(dt = dt[good], ydt = Y_left[good] * dt[good])
  y <- dY[good]
  w <- weights[good]
  fit <- tryCatch(lm.wfit(x = x, y = y, w = w), error = function(e) NULL)
  if (is.null(fit) || length(fit$coefficients) < 2L) {
    return(list(center = NA_real_, lambda = NA_real_))
  }

  co <- fit$coefficients
  lambda_hat <- -unname(co[2])
  center_hat <- unname(co[1]) / lambda_hat
  if (!is.finite(lambda_hat) || lambda_hat <= 0 || !is.finite(center_hat)) {
    return(list(center = NA_real_, lambda = lambda_hat))
  }

  y_med <- median(Y_raw[is.finite(Y_raw)], na.rm = TRUE)
  y_mad <- median(abs(Y_raw[is.finite(Y_raw)] - y_med), na.rm = TRUE)
  if (is.finite(y_mad) && y_mad > 0 &&
      abs(center_hat - y_med) > 10 * y_mad) {
    return(list(center = NA_real_, lambda = lambda_hat))
  }

  list(center = center_hat, lambda = lambda_hat)
}

estimate_sigma1_from_ou_likelihood <- function(Y_centered, dt, weights,
                                               scale_center,
                                               lower_mult = 0.25,
                                               upper_mult = 4) {
  dY <- diff(Y_centered)
  Y_left <- Y_centered[-length(Y_centered)]
  good <- is.finite(dY) & is.finite(Y_left) &
    is.finite(dt) & dt > 0 & is.finite(weights) & weights > 0

  if (sum(good) < 100L || !is.finite(scale_center) || scale_center <= 0) {
    return(list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                objective = NA_real_))
  }

  dY <- dY[good]
  Y_left <- Y_left[good]
  dt_good <- dt[good]
  w <- weights[good]

  score_log_scale <- function(log_s) {
    s <- exp(log_s)
    den <- sum(w * (s^2 + Y_left^2) * dt_good)
    num <- sum(w * dY^2)
    rho2 <- if (is.finite(den) && den > 0) num / den else NA_real_
    if (!is.finite(rho2) || rho2 <= 0) return(Inf)

    var <- rho2 * (s^2 + Y_left^2) * dt_good
    resid <- dY + rho2 * Y_left * dt_good
    ok <- is.finite(var) & var > 0 & is.finite(resid) & is.finite(w) & w > 0
    if (sum(ok) < 100L) return(Inf)
    sum(w[ok] * (log(var[ok]) + resid[ok]^2 / var[ok])) / sum(w[ok])
  }

  lo <- log(scale_center * lower_mult)
  hi <- log(scale_center * upper_mult)
  opt <- tryCatch(optimize(score_log_scale, interval = c(lo, hi)),
                  error = function(e) NULL)
  if (is.null(opt) || !is.finite(opt$minimum) || !is.finite(opt$objective)) {
    return(list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                objective = NA_real_))
  }

  sigma1_hat <- exp(opt$minimum)
  rho2_hat <- sum(w * dY^2) / sum(w * (sigma1_hat^2 + Y_left^2) * dt_good)
  rho_hat <- if (is.finite(rho2_hat) && rho2_hat > 0) sqrt(rho2_hat) else NA_real_

  list(sigma1_hat = sigma1_hat,
       rho_hat = rho_hat,
       objective = opt$objective)
}

reconstruct_v_from_observed <- function(R, dt, alpha_hat, sigma2_hat, k_hat,
                                        max_iter = 5L, burn_frac = 0.10,
                                        mode = "hard", tol = 1e-3,
                                        tail_prob_floor = 0,
                                        tail_prob_cap = 0.01,
                                        threshold_mode = "max",
                                        rho_estimator = "qv",
                                        scale_method = "stationary",
                                        scale_trim_q = 0.95,
                                        scale_n_bins = 25L,
                                        scale_abs_probs = c(0.25, 0.35, 0.45,
                                                            0.55, 0.65, 0.75),
                                        scale_likelihood_lower_mult = 0.25,
                                        scale_likelihood_upper_mult = 4,
                                        dynamic_refilter = FALSE,
                                        dynamic_tail_prob = NA_real_,
                                        dynamic_tail_cap = 0.01,
                                        dynamic_max_frac = 0.02,
                                        dynamic_action = "fill",
                                        dynamic_update_once = FALSE,
                                        center_method = "median") {
  scale_method <- match.arg(scale_method,
                            c("stationary", "stationary_quantile",
                              "stationary_rank", "diffusion_qv",
                              "ou_likelihood"))
  center_method <- match.arg(center_method, c("median", "ou_drift"))
  dynamic_action <- match.arg(dynamic_action, c("fill", "winsor"))
  jf <- make_jump_filter(R, dt, alpha_hat, sigma2_hat, k_hat, mode = mode,
                         tail_prob_floor = tail_prob_floor,
                         tail_prob_cap = tail_prob_cap,
                         threshold_mode = threshold_mode)
  n <- length(R)
  burn_n <- min(max(as.integer(floor(burn_frac * (n + 1))), 0L), n)
  active_is_jump <- jf$is_jump
  active_weights <- jf$weights
  active_weights[active_is_jump] <- 0
  dynamic_added <- rep(FALSE, n)
  dynamic_fill_override <- rep(NA_real_, n)
  dynamic_updates <- 0L
  dynamic_tail_prob_used <- NA_real_
  dynamic_cutoff <- NA_real_
  center_method_used <- "median"
  center_ou_lambda <- NA_real_

  stationary_v_cdf <- function(v) {
    0.5 + (atan(v) + v / (1 + v^2)) / pi
  }
  stationary_v_abs_median <- uniroot(
    function(v) stationary_v_cdf(v) - 0.75,
    interval = c(0, 100)
  )$root
  stationary_v_abs_quantile <- function(p) {
    vapply(p, function(pp) {
      uniroot(
        function(v) stationary_v_cdf(v) - ((1 + pp) / 2),
        interval = c(0, 100)
      )$root
    }, numeric(1))
  }

  center_scale_path <- function(Y_raw) {
    idx <- (burn_n + 1L):length(Y_raw)
    y_use <- Y_raw[idx]
    center <- median(y_use[is.finite(y_use)], na.rm = TRUE)
    center_method_used <<- "median"
    center_ou_lambda <<- NA_real_
    if (center_method == "ou_drift") {
      center_ou <- estimate_center_from_ou_drift(Y_raw, dt, active_weights)
      if (is.finite(center_ou$center)) {
        center <- center_ou$center
        center_method_used <<- "ou_drift"
        center_ou_lambda <<- center_ou$lambda
      }
    }
    Y_centered <- Y_raw - center
    y_scale <- Y_centered[idx]
    sigma1_stationary <- median(abs(y_scale[is.finite(y_scale)]), na.rm = TRUE) /
      stationary_v_abs_median
    if (!is.finite(sigma1_stationary) || sigma1_stationary <= 0) {
      sigma1_stationary <- sd(y_scale[is.finite(y_scale)], na.rm = TRUE)
    }
    if (!is.finite(sigma1_stationary) || sigma1_stationary <= 0) {
      sigma1_stationary <- 1
    }

    probs <- scale_abs_probs
    probs <- probs[is.finite(probs) & probs > 0 & probs < 1]
    sigma1_stationary_quantile <- NA_real_
    if (length(probs) > 0L) {
      y_abs <- abs(y_scale[is.finite(y_scale)])
      if (length(y_abs) > 10L) {
        empirical_q <- as.numeric(quantile(y_abs, probs = probs,
                                           na.rm = TRUE, names = FALSE))
        stationary_q <- stationary_v_abs_quantile(probs)
        ratios <- empirical_q / stationary_q
        sigma1_stationary_quantile <- median(
          ratios[is.finite(ratios) & ratios > 0],
          na.rm = TRUE
        )
      }
    }
    if (!is.finite(sigma1_stationary_quantile) ||
        sigma1_stationary_quantile <= 0) {
      sigma1_stationary_quantile <- sigma1_stationary
    }

    scale_qv <- list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                     slope = NA_real_, intercept = NA_real_)
    scale_likelihood <- list(sigma1_hat = NA_real_, rho_hat = NA_real_,
                             objective = NA_real_)
    V_override <- NULL
    sigma1_hat <- sigma1_stationary
    scale_method_used <- "stationary"

    if (scale_method == "stationary_quantile") {
      sigma1_hat <- sigma1_stationary_quantile
      scale_method_used <- "stationary_quantile"
    }

    if (scale_method == "stationary_rank") {
      y_train <- y_scale[is.finite(y_scale)]
      if (length(y_train) > 100L) {
        y_sorted <- sort(y_train)
        p_all <- (findInterval(Y_centered, y_sorted) + 0.5) /
          (length(y_sorted) + 1)
        p_eps <- 0.5 / (length(y_sorted) + 1)
        p_all <- pmin(pmax(p_all, p_eps), 1 - p_eps)
        V_rank <- qt(p_all, df = 3) / sqrt(3)
        idx_fit <- idx[is.finite(Y_centered[idx]) & is.finite(V_rank[idx])]
        den_rank <- sum(V_rank[idx_fit]^2)
        sigma1_rank <- if (is.finite(den_rank) && den_rank > 0) {
          sum(Y_centered[idx_fit] * V_rank[idx_fit]) / den_rank
        } else {
          NA_real_
        }
        if (is.finite(sigma1_rank) && sigma1_rank > 0) {
          sigma1_hat <- sigma1_rank
          V_override <- V_rank
          scale_method_used <- "stationary_rank"
        }
      }
    }

    if (scale_method == "diffusion_qv") {
      scale_qv <- estimate_sigma1_from_diffusion_qv(
        Y_centered = Y_centered,
        dt = dt,
        weights = active_weights,
        trim_q = scale_trim_q,
        n_bins = scale_n_bins
      )
      if (is.finite(scale_qv$sigma1_hat) && scale_qv$sigma1_hat > 0) {
        sigma1_hat <- scale_qv$sigma1_hat
        scale_method_used <- "diffusion_qv"
      }
    }

    if (scale_method == "ou_likelihood") {
      scale_likelihood <- estimate_sigma1_from_ou_likelihood(
        Y_centered = Y_centered,
        dt = dt,
        weights = active_weights,
        scale_center = sigma1_stationary,
        lower_mult = scale_likelihood_lower_mult,
        upper_mult = scale_likelihood_upper_mult
      )
      if (is.finite(scale_likelihood$sigma1_hat) &&
          scale_likelihood$sigma1_hat > 0) {
        sigma1_hat <- scale_likelihood$sigma1_hat
        scale_method_used <- "ou_likelihood"
      }
    }

    list(Y_centered = Y_centered,
         sigma1_hat = sigma1_hat,
         sigma1_stationary = sigma1_stationary,
         sigma1_stationary_quantile = sigma1_stationary_quantile,
         center = center,
         scale_method_used = scale_method_used,
         V_override = V_override,
         scale_qv = scale_qv,
         scale_likelihood = scale_likelihood)
  }

  dY <- ifelse(active_is_jump, 0, R)
  dY[!is.finite(dY)] <- 0

  sigma1_prev <- NA_real_
  rho_prev <- NA_real_
  converged <- FALSE
  rho_hat <- NA_real_

  for (iter in seq_len(max_iter)) {
    Y_raw <- c(0, cumsum(dY))
    scaled <- center_scale_path(Y_raw)
    sigma1_hat <- scaled$sigma1_hat
    V_hat <- if (!is.null(scaled$V_override)) {
      scaled$V_override
    } else {
      scaled$Y_centered / sigma1_hat
    }
    est <- estimate_sigma1_rho_from_vhat(V_hat, dt, active_weights,
                                         rho_estimator = rho_estimator)
    rho_hat <- est$rho_hat

    filter_changed <- FALSE
    next_is_jump <- active_is_jump
    can_update_dynamic <- isTRUE(dynamic_refilter) && is.finite(rho_hat) && rho_hat > 0 &&
      (!isTRUE(dynamic_update_once) || dynamic_updates < 1L)
    if (can_update_dynamic) {
      V_left_hat <- V_hat[-length(V_hat)]
      resid <- diff(V_hat) + rho_hat^2 * V_left_hat * dt
      resid_scale <- rho_hat * sqrt((1 + V_left_hat^2) * dt)
      std_resid <- abs(resid) / resid_scale

      p_dyn <- dynamic_tail_prob
      if (!is.finite(p_dyn) || p_dyn <= 0) p_dyn <- jf$p_tail
      p_dyn <- max(p_dyn, 1 / (n + 1))
      p_dyn <- min(p_dyn, dynamic_tail_cap)
      dynamic_tail_prob_used <- p_dyn
      dynamic_cutoff <- qnorm(1 - p_dyn / 2)

      dyn_candidates <- is.finite(std_resid) & std_resid > dynamic_cutoff
      max_dyn <- if (is.finite(dynamic_max_frac) && dynamic_max_frac > 0) {
        max(1L, as.integer(floor(dynamic_max_frac * n)))
      } else {
        n
      }
      if (sum(dyn_candidates) > max_dyn) {
        dyn_vals <- sort(std_resid[dyn_candidates], decreasing = TRUE)
        dyn_floor <- dyn_vals[max_dyn]
        dyn_candidates <- dyn_candidates & std_resid >= dyn_floor
      }

      next_is_jump <- jf$is_jump | dyn_candidates
      filter_changed <- any(next_is_jump != active_is_jump)
      dynamic_added <- next_is_jump & !jf$is_jump

      dynamic_fill_override <- rep(NA_real_, n)
      if (dynamic_action == "winsor" && any(dynamic_added)) {
        capped_resid <- sign(resid) * pmin(abs(resid), dynamic_cutoff * resid_scale)
        dynamic_fill_override[dynamic_added] <-
          sigma1_hat * (-rho_hat^2 * V_left_hat[dynamic_added] * dt[dynamic_added] +
                          capped_resid[dynamic_added])
      }
      dynamic_updates <- dynamic_updates + 1L
    }

    rel_sigma1 <- if (is.finite(sigma1_prev)) {
      abs(sigma1_hat - sigma1_prev) / max(abs(sigma1_prev), 1e-12)
    } else {
      Inf
    }
    rel_rho <- if (is.finite(rho_prev) && is.finite(rho_hat)) {
      abs(rho_hat - rho_prev) / max(abs(rho_prev), 1e-12)
    } else {
      Inf
    }

    if (iter > 1L && rel_sigma1 < tol && rel_rho < tol && !filter_changed) {
      converged <- TRUE
      break
    }

    if (!is.finite(rho_hat)) break

    if (filter_changed) {
      active_is_jump <- next_is_jump
      active_weights <- jf$weights
      active_weights[active_is_jump] <- 0
    }

    V_left_hat <- V_hat[-length(V_hat)]
    fill_jump <- sigma1_hat * (-rho_hat^2 * V_left_hat * dt)
    fill_jump[!is.finite(fill_jump)] <- 0
    dY <- R
    dY[!is.finite(dY)] <- 0
    dY[active_is_jump] <- fill_jump[active_is_jump]
    if (dynamic_action == "winsor" && any(is.finite(dynamic_fill_override))) {
      use_override <- is.finite(dynamic_fill_override)
      dY[use_override] <- dynamic_fill_override[use_override]
    }

    sigma1_prev <- sigma1_hat
    rho_prev <- rho_hat
  }

  jf_out <- jf
  jf_out$initial_is_jump <- jf$is_jump
  jf_out$is_jump <- active_is_jump
  jf_out$weights <- active_weights
  jf_out$dynamic_is_jump <- dynamic_added
  jf_out$dynamic_refilter <- isTRUE(dynamic_refilter)
  jf_out$dynamic_action <- dynamic_action
  jf_out$dynamic_update_once <- isTRUE(dynamic_update_once)
  jf_out$dynamic_updates <- dynamic_updates
  jf_out$dynamic_tail_prob <- dynamic_tail_prob_used
  jf_out$dynamic_cutoff <- dynamic_cutoff

  list(
    V_hat = V_hat,
    dV_hat = diff(V_hat),
    sigma1_hat = sigma1_hat,
    sigma1_stationary = scaled$sigma1_stationary,
    sigma1_stationary_quantile = scaled$sigma1_stationary_quantile,
    rho_hat = rho_hat,
    rho_qv = est$rho_qv,
    rho_contrast = est$rho_contrast,
    rho_robust = est$rho_robust,
    scale_method = scale_method,
    scale_method_used = scaled$scale_method_used,
    center_method = center_method,
    center_method_used = center_method_used,
    center = scaled$center,
    center_ou_lambda = center_ou_lambda,
    scale_rho_qv = scaled$scale_qv$rho_hat,
    scale_qv_slope = scaled$scale_qv$slope,
    scale_qv_intercept = scaled$scale_qv$intercept,
    scale_likelihood_rho = scaled$scale_likelihood$rho_hat,
    scale_likelihood_objective = scaled$scale_likelihood$objective,
    jump_filter = jf_out,
    rho_estimator = rho_estimator,
    converged = converged,
    n_iter = iter
  )
}

reconstruct_v_joint_tail_odds <- function(R, dt, alpha_hat, sigma2_hat, k_hat,
                                          max_iter = 8L, burn_frac = 0.10,
                                          tol = 1e-3,
                                          tail_prob_floor = 0,
                                          tail_prob_cap = 0.01,
                                          threshold_mode = "max",
                                          min_good_frac = 0.50,
                                          damping_eta = 1) {
  n <- length(R)
  if (!is.finite(damping_eta)) damping_eta <- 1
  damping_eta <- min(max(damping_eta, 0), 1)
  jf <- make_jump_filter(R, dt, alpha_hat, sigma2_hat, k_hat,
                         mode = "hard",
                         tail_prob_floor = tail_prob_floor,
                         tail_prob_cap = tail_prob_cap,
                         threshold_mode = threshold_mode)

  base <- reconstruct_v_from_observed(
    R = R,
    dt = dt,
    alpha_hat = alpha_hat,
    sigma2_hat = sigma2_hat,
    k_hat = k_hat,
    max_iter = 5L,
    burn_frac = burn_frac,
    mode = "hard",
    tail_prob_floor = tail_prob_floor,
    tail_prob_cap = tail_prob_cap,
    threshold_mode = threshold_mode,
    rho_estimator = "qv",
    scale_method = "stationary",
    dynamic_refilter = TRUE,
    dynamic_action = "winsor",
    dynamic_update_once = TRUE
  )

  stationary_v_cdf <- function(v) {
    0.5 + (atan(v) + v / (1 + v^2)) / pi
  }
  stationary_v_abs_median <- uniroot(
    function(v) stationary_v_cdf(v) - 0.75,
    interval = c(0, 100)
  )$root
  center_scale_y_path <- function(Y_raw) {
    burn_n <- min(max(as.integer(floor(burn_frac * length(Y_raw))), 0L),
                  length(Y_raw) - 2L)
    idx <- (burn_n + 1L):length(Y_raw)
    y_use <- Y_raw[idx]
    center <- median(y_use[is.finite(y_use)], na.rm = TRUE)
    Y_centered <- Y_raw - center
    y_scale <- Y_centered[idx]
    sigma1_hat <- median(abs(y_scale[is.finite(y_scale)]), na.rm = TRUE) /
      stationary_v_abs_median
    if (!is.finite(sigma1_hat) || sigma1_hat <= 0) {
      sigma1_hat <- sd(y_scale[is.finite(y_scale)], na.rm = TRUE)
    }
    if (!is.finite(sigma1_hat) || sigma1_hat <= 0) sigma1_hat <- 1
    list(Y_centered = Y_centered,
         V_hat = Y_centered / sigma1_hat,
         sigma1_hat = sigma1_hat,
         center = center)
  }

  V_hat <- base$V_hat
  sigma1_hat <- base$sigma1_hat
  rho_hat <- base$rho_hat
  dY_base <- sigma1_hat * diff(V_hat)
  dY_base[!is.finite(dY_base)] <- 0
  c_alpha <- gamma(alpha_hat) * sin(pi * alpha_hat / 2) / pi
  dt_pos <- ifelse(is.finite(dt) & dt > 0, dt,
                   median(dt[is.finite(dt) & dt > 0], na.rm = TRUE))
  tiny <- sqrt(.Machine$double.eps)
  z_cut <- qnorm(1 - jf$p_tail / 2)
  clip <- jf$is_jump
  jump_like <- rep(FALSE, n)
  jump_resid <- rep(0, n)
  cap <- jf$threshold
  scaled <- list(center = base$center)
  converged <- FALSE
  prev_sigma1 <- sigma1_hat
  prev_rho <- rho_hat

  for (iter in seq_len(max_iter)) {
    if (!is.finite(sigma1_hat) || sigma1_hat <= 0 ||
        !is.finite(rho_hat) || rho_hat <= 0) {
      break
    }

    V_left <- V_hat[-length(V_hat)]
    mu_y <- sigma1_hat * (-rho_hat^2 * V_left * dt_pos)
    sd_y <- sigma1_hat * rho_hat * sqrt(pmax(1 + V_left^2, 0) * dt_pos)
    sd_y[!is.finite(sd_y) | sd_y <= 0] <- tiny
    brownian_cap <- z_cut * sd_y

    resid <- R - mu_y
    abs_resid <- abs(resid)
    stable_tail <- 2 * c_alpha * sigma2_hat^alpha_hat * dt_pos /
      pmax(abs_resid, tiny)^alpha_hat
    stable_tail <- pmin(pmax(stable_tail, 0), 1)
    gaussian_tail <- 2 * pnorm(-abs_resid / sd_y)
    gaussian_tail <- pmax(gaussian_tail, .Machine$double.xmin)

    jump_like <- is.finite(abs_resid) &
      (abs_resid > brownian_cap) &
      (stable_tail >= gaussian_tail)
    clip <- jf$is_jump | jump_like

    cap <- brownian_cap
    bad_cap <- !is.finite(cap) | cap <= 0
    cap[bad_cap] <- jf$threshold[bad_cap]
    cap[!is.finite(cap) | cap <= 0] <- Inf
    latent_resid <- resid
    latent_resid[clip] <- sign(resid[clip]) *
      pmin(abs_resid[clip], cap[clip])

    dY_joint <- mu_y + latent_resid
    dY_joint[!is.finite(dY_joint)] <- 0
    dY <- (1 - damping_eta) * dY_base + damping_eta * dY_joint
    dY[!is.finite(dY)] <- 0
    jump_resid <- R - dY
    jump_resid[!is.finite(jump_resid)] <- 0

    scaled <- center_scale_y_path(c(0, cumsum(dY)))
    V_new <- scaled$V_hat
    weights <- as.numeric(!clip & is.finite(R))
    if (mean(weights > 0) < min_good_frac) {
      weights <- base$jump_filter$weights
    }
    est <- estimate_sigma1_rho_from_vhat(V_new, dt, weights,
                                         rho_estimator = "qv")

    sigma1_new <- scaled$sigma1_hat
    rho_new <- est$rho_hat
    rel_sigma1 <- abs(sigma1_new - prev_sigma1) / max(abs(prev_sigma1), 1e-12)
    rel_rho <- if (is.finite(rho_new) && is.finite(prev_rho)) {
      abs(rho_new - prev_rho) / max(abs(prev_rho), 1e-12)
    } else {
      Inf
    }

    V_hat <- V_new
    sigma1_hat <- sigma1_new
    if (is.finite(rho_new) && rho_new > 0) rho_hat <- rho_new
    if (iter > 1L && rel_sigma1 < tol && rel_rho < tol) {
      converged <- TRUE
      break
    }
    prev_sigma1 <- sigma1_hat
    prev_rho <- rho_hat
  }

  weights <- as.numeric(!clip & is.finite(R))
  if (mean(weights > 0) < min_good_frac) {
    weights <- base$jump_filter$weights
  }
  est <- estimate_sigma1_rho_from_vhat(V_hat, dt, weights,
                                       rho_estimator = "qv")

  jf_out <- jf
  jf_out$initial_is_jump <- jf$is_jump
  jf_out$is_jump <- clip
  jf_out$weights <- weights
  jf_out$dynamic_is_jump <- rep(FALSE, n)
  jf_out$dynamic_refilter <- FALSE
  jf_out$dynamic_action <- "none"
  jf_out$dynamic_update_once <- TRUE
  jf_out$dynamic_updates <- 0L
  jf_out$dynamic_tail_prob <- NA_real_
  jf_out$dynamic_cutoff <- NA_real_
  jf_out$joint_tail_odds_is_jump <- jump_like
  jf_out$joint_tail_odds_clip_frac <- mean(jump_like, na.rm = TRUE)
  jf_out$joint_clip_frac <- mean(clip, na.rm = TRUE)
  jf_out$joint_brownian_cutoff <- z_cut
  jf_out$joint_jump_resid <- jump_resid
  jf_out$joint_damping_eta <- damping_eta
  jf_out$base_dynamic_added_frac <- mean(base$jump_filter$dynamic_is_jump,
                                         na.rm = TRUE)

  list(
    V_hat = V_hat,
    dV_hat = diff(V_hat),
    sigma1_hat = sigma1_hat,
    sigma1_stationary = sigma1_hat,
    sigma1_stationary_quantile = NA_real_,
    rho_hat = est$rho_hat,
    rho_qv = est$rho_qv,
    rho_contrast = est$rho_contrast,
    rho_robust = est$rho_robust,
    scale_method = "joint_tail_odds_cap",
    scale_method_used = "stationary",
    center_method = "median",
    center_method_used = "median",
    center = scaled$center,
    center_ou_lambda = NA_real_,
    scale_rho_qv = NA_real_,
    scale_qv_slope = NA_real_,
    scale_qv_intercept = NA_real_,
    scale_likelihood_rho = NA_real_,
    scale_likelihood_objective = NA_real_,
    jump_filter = jf_out,
    rho_estimator = "qv",
    converged = converged,
    n_iter = iter
  )
}

alpha_v_input <- alpha_hat
sigma2_v_input <- sig2_out$sigma2_hat
using_truth_jump_params <- FALSE

if (env_bool("USE_TRUE_JUMP_PARAMS_FOR_V", FALSE)) {
  alpha_v_input <- alpha_true
  sigma2_v_input <- sigma2_true
  using_truth_jump_params <- TRUE
} else if (!is.finite(alpha_v_input) || alpha_v_input <= 0 || alpha_v_input >= 2 ||
           !is.finite(sigma2_v_input) || sigma2_v_input <= 0) {
  stop("Estimated alpha/sigma2 are unusable for observed-data V reconstruction. Rerun with better first-stage estimates or explicitly set USE_TRUE_JUMP_PARAMS_FOR_V=TRUE for simulation diagnostics only.")
}

v_recon_method <- Sys.getenv("V_RECON_METHOD", unset = "dynamic_winsor")
if (v_recon_method == "joint_tail_odds_cap") {
  v_recon_out <- reconstruct_v_joint_tail_odds(
    R = R,
    dt = dt,
    alpha_hat = alpha_v_input,
    sigma2_hat = sigma2_v_input,
    k_hat = k_hat,
    max_iter = as.integer(env_num("V_RECON_JOINT_MAX_ITER", 8)),
    burn_frac = env_num("V_RECON_BURN_FRAC", 0.10),
    tail_prob_floor = env_num("TAIL_PROB_FLOOR", 0),
    tail_prob_cap = env_num("TAIL_PROB_CAP", 0.01),
    threshold_mode = Sys.getenv("TAIL_THRESHOLD_MODE", unset = "max"),
    min_good_frac = env_num("V_RECON_JOINT_MIN_GOOD_FRAC", 0.50),
    damping_eta = env_num("V_RECON_JOINT_DAMPING_ETA", 1)
  )
} else {
  if (v_recon_method != "dynamic_winsor") {
    stop("Unknown V_RECON_METHOD. Use 'dynamic_winsor' or 'joint_tail_odds_cap'.")
  }
  v_recon_out <- reconstruct_v_from_observed(
    R = R,
    dt = dt,
    alpha_hat = alpha_v_input,
    sigma2_hat = sigma2_v_input,
    k_hat = k_hat,
    max_iter = as.integer(env_num("V_RECON_MAX_ITER", 5)),
    burn_frac = env_num("V_RECON_BURN_FRAC", 0.10),
    mode = Sys.getenv("V_RECON_MODE", unset = "hard"),
    tail_prob_floor = env_num("TAIL_PROB_FLOOR", 0),
    tail_prob_cap = env_num("TAIL_PROB_CAP", 0.01),
    threshold_mode = Sys.getenv("TAIL_THRESHOLD_MODE", unset = "max"),
    rho_estimator = Sys.getenv("RHO_ESTIMATOR", unset = "qv"),
    scale_method = Sys.getenv("V_RECON_SCALE_METHOD", unset = "stationary"),
    scale_trim_q = env_num("V_RECON_SCALE_TRIM_Q", 0.95),
    scale_n_bins = as.integer(env_num("V_RECON_SCALE_N_BINS", 25)),
    scale_likelihood_lower_mult = env_num("V_RECON_LIK_SCALE_LOWER_MULT", 0.25),
    scale_likelihood_upper_mult = env_num("V_RECON_LIK_SCALE_UPPER_MULT", 4),
    dynamic_refilter = env_bool("V_RECON_DYNAMIC_REFILTER", TRUE),
    dynamic_tail_prob = env_num("V_RECON_DYNAMIC_TAIL_PROB", NA_real_),
    dynamic_tail_cap = env_num("V_RECON_DYNAMIC_TAIL_CAP", 0.01),
    dynamic_max_frac = env_num("V_RECON_DYNAMIC_MAX_FRAC", 0.01),
    dynamic_action = Sys.getenv("V_RECON_DYNAMIC_ACTION", unset = "winsor"),
    dynamic_update_once = env_bool("V_RECON_DYNAMIC_UPDATE_ONCE", TRUE),
    center_method = Sys.getenv("V_RECON_CENTER_METHOD", unset = "median")
  )
}

sigma1_hat <- v_recon_out$sigma1_hat
rho_hat <- v_recon_out$rho_hat

cat(sprintf("V reconstruction method = %s\n", v_recon_method))
cat(sprintf("Estimated Sigma1 (V reconstructed) = %.6f\n", sigma1_hat))
cat(sprintf("True sigma1 = %.4f\n", sigma1_true))
cat(sprintf("Estimated Rho (V reconstructed)   = %.6f\n", rho_hat))
cat(sprintf("Estimated Rho QV diagnostic       = %.6f\n", v_recon_out$rho_qv))
cat(sprintf("Estimated Rho contrast diagnostic = %.6f\n", v_recon_out$rho_contrast))
cat(sprintf("Estimated Rho robust diagnostic   = %.6f\n", v_recon_out$rho_robust))
cat(sprintf("True rho = %.4f\n", rho_true))
cat(sprintf("V reconstruction scale method = %s (used %s)\n",
            v_recon_out$scale_method,
            v_recon_out$scale_method_used))
cat(sprintf("V reconstruction center method = %s (used %s, center %.6f)\n",
            v_recon_out$center_method,
            v_recon_out$center_method_used,
            v_recon_out$center))
cat(sprintf("Stationary-scale Sigma1 diagnostic = %.6f\n",
            v_recon_out$sigma1_stationary))
cat(sprintf("Stationary-quantile Sigma1 diagnostic = %.6f\n",
            v_recon_out$sigma1_stationary_quantile))
if (is.finite(v_recon_out$scale_qv_intercept) &&
    is.finite(v_recon_out$scale_qv_slope) &&
    v_recon_out$scale_qv_slope > 0) {
  cat(sprintf("Diffusion-QV Sigma1/Rho diagnostics = %.6f / %.6f\n",
              sqrt(v_recon_out$scale_qv_intercept /
                     v_recon_out$scale_qv_slope),
              v_recon_out$scale_rho_qv))
} else {
  cat("Diffusion-QV Sigma1/Rho diagnostics = unavailable\n")
}
if (is.finite(v_recon_out$scale_likelihood_rho)) {
  cat(sprintf("OU-likelihood scale Rho/objective diagnostics = %.6f / %.6f\n",
              v_recon_out$scale_likelihood_rho,
              v_recon_out$scale_likelihood_objective))
} else {
  cat("OU-likelihood scale diagnostics = unavailable\n")
}
cat(sprintf("Filtered jump increments = %d / %d (%.4f%%)\n",
            sum(v_recon_out$jump_filter$is_jump),
            length(R),
            100 * mean(v_recon_out$jump_filter$is_jump)))
cat(sprintf("Dynamic refilter added = %d / %d (%.4f%%) | enabled = %s | action = %s | updates = %d\n",
            sum(v_recon_out$jump_filter$dynamic_is_jump),
            length(R),
            100 * mean(v_recon_out$jump_filter$dynamic_is_jump),
            v_recon_out$jump_filter$dynamic_refilter,
            v_recon_out$jump_filter$dynamic_action,
            v_recon_out$jump_filter$dynamic_updates))
if (!is.null(v_recon_out$jump_filter$joint_tail_odds_is_jump)) {
  cat(sprintf("Joint tail-odds residual clips = %d / %d (%.4f%%) | total clipped = %.4f%% | Brownian cutoff = %.4f\n",
              sum(v_recon_out$jump_filter$joint_tail_odds_is_jump),
              length(R),
              100 * mean(v_recon_out$jump_filter$joint_tail_odds_is_jump),
              100 * mean(v_recon_out$jump_filter$is_jump),
              v_recon_out$jump_filter$joint_brownian_cutoff))
  cat(sprintf("Joint method base dynamic-added diagnostic = %.4f%%\n",
              100 * v_recon_out$jump_filter$base_dynamic_added_frac))
  cat(sprintf("Joint damping eta = %.4f\n",
              v_recon_out$jump_filter$joint_damping_eta))
}
cat(sprintf("Jump threshold mode = %s | Rho estimator = %s\n",
            v_recon_out$jump_filter$threshold_mode,
            v_recon_out$rho_estimator))
cat(sprintf("V reconstruction iterations = %d | converged = %s\n",
            v_recon_out$n_iter,
            v_recon_out$converged))
if (using_truth_jump_params) {
  cat("V reconstruction jump parameters: true values used for simulation diagnostic only.\n")
}

# Oracle diagnostic only: this block uses simulated V and is not part of estimation.
if (exists("V") && length(V) == length(X)) {
  dV_diag <- diff(V)
  V_left_diag <- V[-length(V)]
  u_thr <- if ((k_hat + 1L) <= m) A[k_hat + 1L] else A[m]
  S <- which(abs(R) <= u_thr)
  num_rho_oracle <- sum(dV_diag[S]^2)
  den_rho_oracle <- sum((1 + V_left_diag[S]^2) * dt[S])
  rho2_oracle <- num_rho_oracle / den_rho_oracle
  rho_oracle <- sqrt(max(rho2_oracle, 0))
  rcov_RV_oracle <- sum(R[S] * dV_diag[S])
  rvar_V_oracle  <- sum(dV_diag[S]^2)
  sigma1_oracle <- rcov_RV_oracle / rvar_V_oracle

  cat(sprintf("Oracle diagnostic Sigma1 = %.6f\n", sigma1_oracle))
  cat(sprintf("Oracle diagnostic Rho   = %.6f\n", rho_oracle))
  burn_diag <- min(max(as.integer(floor(0.10 * length(V))), 0L), length(V) - 2L)
  diag_idx <- (burn_diag + 1L):min(length(V), length(v_recon_out$V_hat))
  v_corr <- suppressWarnings(cor(v_recon_out$V_hat[diag_idx], V[diag_idx],
                                 use = "complete.obs"))
  cat(sprintf("Oracle diagnostic cor(V_hat, V) = %.6f\n", v_corr))
} else {
  cat("Oracle diagnostic skipped: simulated V path is unavailable.\n")
}
