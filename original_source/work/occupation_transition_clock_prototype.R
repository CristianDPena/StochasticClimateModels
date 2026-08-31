# Prototype for the generator-calibrated occupation clock method.
#
# This file is intentionally standalone. It does not reconstruct V_t, does not
# call the particle likelihood, and does not use beta-rho profile logic. It
# estimates beta from the occupation distribution of jump-deconvolved block
# activity and estimates rho from generator/martingale estimating equations.

symmetric_stable_r <- function(n, alpha) {
  if (abs(alpha - 1) < 1e-8) {
    u <- stats::runif(n, -pi / 2, pi / 2)
    w <- stats::rexp(n)
    return((2 / pi) * ((pi / 2 + u) * tan(u) -
                         log((pi / 2 * w * cos(u)) / (pi / 2 + u))))
  }
  u <- stats::runif(n, -pi / 2, pi / 2)
  w <- stats::rexp(n)
  sin(alpha * u) / (cos(u) ^ (1 / alpha)) *
    (cos((1 - alpha) * u) / w) ^ ((1 - alpha) / alpha)
}

simulate_residual_euler <- function(n_steps, terminal, alpha, sigma2,
                                    sigma1, rho, seed = 1L,
                                    keep_v = FALSE) {
  set.seed(seed)
  n <- as.integer(n_steps)
  dt_val <- terminal / n
  dt <- rep(dt_val, n)
  v <- numeric(n + 1L)
  r <- numeric(n)
  d_w <- stats::rnorm(n, sd = sqrt(dt_val))
  d_l <- dt_val ^ (1 / alpha) * symmetric_stable_r(n, alpha)
  rho2 <- rho ^ 2
  for (i in seq_len(n)) {
    d_v <- -rho2 * v[i] * dt_val +
      rho * sqrt(pmax(1 + v[i] ^ 2, 1e-12)) * d_w[i]
    v[i + 1L] <- v[i] + d_v
    r[i] <- sigma1 * d_v + sigma2 * d_l[i]
  }
  if (keep_v) {
    list(R = r, dt = dt, V = v)
  } else {
    list(R = r, dt = dt)
  }
}

hill_alpha_path <- function(y_sorted_desc) {
  m <- length(y_sorted_desc)
  if (m < 5L) stop("Not enough tail points for Hill estimation.")
  log_y <- log(y_sorted_desc)
  k <- seq_len(m - 1L)
  (cumsum(log_y[seq_len(m - 1L)]) - k * log_y[2:m]) / k
}

choose_hill_plateau <- function(alpha_path, alpha_cap = 1.98,
                                k_min = 20L, k_max_cap = 20000L,
                                n_candidates = 60L) {
  m <- length(alpha_path)
  k_max <- min(as.integer(m / 100L), as.integer(k_max_cap), m - 10L)
  if (!is.finite(k_max) || k_max < k_min) {
    fallback <- min(max(2L, as.integer(floor(m / 10L))), m)
    return(list(k = fallback, score = NA_real_, range = c(1L, m),
                small_sample = TRUE))
  }
  cs <- c(0, cumsum(alpha_path))
  cs2 <- c(0, cumsum(alpha_path ^ 2))
  k_grid <- unique(as.integer(round(exp(seq(log(k_min), log(k_max),
                                            length.out = n_candidates)))))
  k_grid <- k_grid[k_grid >= k_min & k_grid <= k_max]
  best_k <- NA_integer_
  best_score <- Inf
  for (k in k_grid) {
    lo <- max(1L, as.integer(k / 2L))
    len <- k - lo + 1L
    seg_sum <- cs[k + 1L] - cs[lo]
    seg_sum2 <- cs2[k + 1L] - cs2[lo]
    seg_mean <- seg_sum / len
    seg_sd <- sqrt(max(seg_sum2 / len - seg_mean ^ 2, 0))
    if (!is.finite(seg_sd) || alpha_path[k] > alpha_cap) next
    if (seg_sd < best_score) {
      best_score <- seg_sd
      best_k <- k
    }
  }
  if (is.na(best_k)) best_k <- k_min
  list(k = best_k, score = best_score, range = range(k_grid),
       small_sample = FALSE)
}

estimate_sigma2_tail_robust <- function(alpha_hat, A_sorted_desc, k_hat,
                                        dt_vec, w_half = 20L,
                                        k_bounds = NULL) {
  m <- length(A_sorted_desc)
  dt_robust <- stats::median(dt_vec)
  c_alpha <- gamma(alpha_hat) * sin(pi * alpha_hat / 2) / pi
  if (is.null(k_bounds)) {
    k_lo <- max(2L, k_hat - as.integer(w_half))
    k_hi <- min(m - 1L, k_hat + as.integer(w_half))
  } else {
    k_lo <- max(2L, as.integer(k_bounds[1L]))
    k_hi <- min(m - 1L, as.integer(k_bounds[2L]))
  }
  k_grid <- k_lo:k_hi
  u <- A_sorted_desc[k_grid]
  p <- (k_grid - 0.5) / m
  sigma2_alpha <- p * u ^ alpha_hat / (2 * c_alpha * dt_robust)
  estimates <- ifelse(is.finite(sigma2_alpha) & sigma2_alpha > 0,
                      sigma2_alpha ^ (1 / alpha_hat), NA_real_)
  list(
    sigma2_hat = stats::median(estimates[is.finite(estimates)], na.rm = TRUE),
    details = list(k_hat = k_hat, k_lo = k_lo, k_hi = k_hi,
                   dtrob = dt_robust, c_alpha = c_alpha,
                   window_size = length(k_grid))
  )
}

estimate_alpha_sigma2_hill <- function(R, dt, alpha_cap = 1.98,
                                       median_width = 5L) {
  n <- min(length(R), length(dt))
  finite <- is.finite(R[seq_len(n)]) & R[seq_len(n)] != 0 &
    is.finite(dt[seq_len(n)]) & dt[seq_len(n)] > 0
  A <- sort(abs(R[seq_len(n)][finite]), decreasing = TRUE)
  if (length(A) < 100L) stop("Not enough finite residuals for tail estimation.")
  gamma_raw <- hill_alpha_path(A)
  alpha_raw <- 1 / gamma_raw
  median_width <- max(1L, as.integer(median_width))
  if (median_width %% 2L == 0L) median_width <- median_width + 1L
  alpha_path <- if (median_width > 1L && length(alpha_raw) >= median_width) {
    stats::runmed(alpha_raw, k = median_width, endrule = "median")
  } else {
    alpha_raw
  }
  selection <- choose_hill_plateau(alpha_path, alpha_cap = alpha_cap)
  k_hat <- selection$k
  alpha_hat <- alpha_path[k_hat]
  if (!is.finite(alpha_hat) || alpha_hat <= 0 || alpha_hat >= 2) {
    stop(sprintf("Hill alpha estimate outside stable range: %.6g", alpha_hat))
  }
  sigma2 <- estimate_sigma2_tail_robust(alpha_hat, A, k_hat,
                                        dt[seq_len(n)][finite])
  if (!is.finite(sigma2$sigma2_hat) || sigma2$sigma2_hat <= 0) {
    stop(sprintf("Invalid sigma2 estimate: %.6g", sigma2$sigma2_hat))
  }
  list(
    alpha_hat = alpha_hat,
    sigma2_hat = sigma2$sigma2_hat,
    k_hat = k_hat,
    hill_score = selection$score,
    hill_range = selection$range,
    hill_small_sample = selection$small_sample,
    sigma2_details = sigma2$details
  )
}

stationary_v_sample <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- numeric(0)
  while (length(out) < n) {
    m <- ceiling(1.5 * (n - length(out)) + 100)
    theta <- stats::runif(m, -pi / 2, pi / 2)
    accept <- stats::runif(m) <= cos(theta) ^ 2
    out <- c(out, tan(theta[accept]))
  }
  out[seq_len(n)]
}

stationary_q_sample <- function(n, seed = NULL) {
  1 + stationary_v_sample(n, seed = seed) ^ 2
}

stationary_q_quantile <- function(prob) {
  prob <- as.numeric(prob)
  if (length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1) {
    stop("prob must be a scalar strictly between zero and one.")
  }
  target_v_cdf <- (1 + prob) / 2
  v_cdf <- function(v) {
    0.5 + (atan(v) + v / (1 + v ^ 2)) / pi
  }
  upper <- 1
  while (v_cdf(upper) < target_v_cdf) upper <- 2 * upper
  v <- stats::uniroot(function(x) v_cdf(x) - target_v_cdf,
                       interval = c(0, upper))$root
  1 + v ^ 2
}

ecf_frequency_multipliers <- function(block_size,
                                      constants = c(0.75, 1.50, 2.25),
                                      exponent = 1 / 6) {
  block_size <- as.numeric(block_size)
  if (!is.finite(block_size) || block_size < 10) {
    stop("block_size must be at least 10.")
  }
  constants * block_size ^ (-exponent)
}

block_cf_activity <- function(R, dt, alpha, sigma2, block_size,
                              freq_mults = c(0.50, 0.75, 1.00),
                              min_mod = 1e-8) {
  n <- min(length(R), length(dt))
  R <- R[seq_len(n)]
  dt <- dt[seq_len(n)]
  block_size <- as.integer(block_size)
  if (!is.finite(block_size) || block_size < 10L) {
    stop("block_size must be at least 10.")
  }
  n_blocks <- floor(n / block_size)
  if (n_blocks < 5L) stop("Not enough blocks for occupation statistics.")
  R <- R[seq_len(n_blocks * block_size)]
  dt <- dt[seq_len(n_blocks * block_size)]

  Z <- R / sqrt(dt)
  scale_z <- stats::median(abs(Z[is.finite(Z)]), na.rm = TRUE)
  if (!is.finite(scale_z) || scale_z <= 0) {
    scale_z <- stats::mad(Z, constant = 1.4826, na.rm = TRUE)
  }
  if (!is.finite(scale_z) || scale_z <= 0) {
    scale_z <- stats::sd(Z, na.rm = TRUE)
  }
  if (!is.finite(scale_z) || scale_z <= 0) {
    stop("Could not determine a finite residual scale.")
  }
  freqs <- freq_mults / scale_z
  freqs <- freqs[is.finite(freqs) & freqs > 0]
  if (!length(freqs)) stop("No valid characteristic-function frequencies.")

  activity <- rep(NA_real_, n_blocks)
  activity_by_freq <- matrix(NA_real_, nrow = n_blocks, ncol = length(freqs))
  for (j in seq_len(n_blocks)) {
    idx <- ((j - 1L) * block_size + 1L):(j * block_size)
    z_j <- Z[idx]
    jump_power_j <- mean(dt[idx] ^ (1 - alpha / 2), na.rm = TRUE)
    vals <- rep(NA_real_, length(freqs))
    for (m in seq_along(freqs)) {
      u <- freqs[m]
      phi <- mean(exp(1i * u * z_j))
      log_mod <- log(max(Mod(phi), min_mod))
      jump_attenuation <- sigma2 ^ alpha * abs(u) ^ alpha * jump_power_j
      vals[m] <- -2 * (log_mod + jump_attenuation) / (u ^ 2)
    }
    activity_by_freq[j, ] <- vals
    vals_ok <- vals[is.finite(vals) & vals > 0]
    if (length(vals_ok)) activity[j] <- stats::median(vals_ok)
  }
  rel_iqr <- apply(activity_by_freq, 1L, function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) < 2L) return(NA_real_)
    med <- stats::median(x)
    if (!is.finite(med) || med <= 0) return(NA_real_)
    stats::IQR(x) / med
  })

  list(
    activity = activity,
    activity_by_freq = activity_by_freq,
    freq_rel_iqr = rel_iqr,
    freq_rel_iqr_median = stats::median(rel_iqr, na.rm = TRUE),
    block_h = mean(dt) * block_size,
    block_size = block_size,
    freqs = freqs,
    n_blocks = n_blocks
  )
}

block_cf_log_observations <- function(R, dt, alpha, sigma2, block_size,
                                      freq_mults = c(0.50, 0.75, 1.00),
                                      min_mod = 1e-8,
                                      min_cov_var = 1e-10) {
  n <- min(length(R), length(dt))
  R <- R[seq_len(n)]
  dt <- dt[seq_len(n)]
  block_size <- as.integer(block_size)
  if (!is.finite(block_size) || block_size < 10L) {
    stop("block_size must be at least 10.")
  }
  n_blocks <- floor(n / block_size)
  if (n_blocks < 5L) stop("Not enough blocks for ECF observations.")
  R <- R[seq_len(n_blocks * block_size)]
  dt <- dt[seq_len(n_blocks * block_size)]

  Z <- R / sqrt(dt)
  scale_z <- stats::median(abs(Z[is.finite(Z)]), na.rm = TRUE)
  if (!is.finite(scale_z) || scale_z <= 0) {
    scale_z <- stats::mad(Z, constant = 1.4826, na.rm = TRUE)
  }
  if (!is.finite(scale_z) || scale_z <= 0) {
    scale_z <- stats::sd(Z, na.rm = TRUE)
  }
  if (!is.finite(scale_z) || scale_z <= 0) {
    stop("Could not determine a finite residual scale.")
  }
  freqs <- freq_mults / scale_z
  freqs <- freqs[is.finite(freqs) & freqs > 0]
  if (!length(freqs)) stop("No valid characteristic-function frequencies.")

  n_freq <- length(freqs)
  y_raw <- matrix(NA_real_, nrow = n_blocks, ncol = n_freq)
  y_dejumped <- matrix(NA_real_, nrow = n_blocks, ncol = n_freq)
  jump_attenuation <- matrix(NA_real_, nrow = n_blocks, ncol = n_freq)
  log_mod_bias <- matrix(NA_real_, nrow = n_blocks, ncol = n_freq)
  phi_hat <- matrix(NA_complex_, nrow = n_blocks, ncol = n_freq)
  obs_cov <- array(NA_real_, dim = c(n_blocks, n_freq, n_freq))

  for (j in seq_len(n_blocks)) {
    idx <- ((j - 1L) * block_size + 1L):(j * block_size)
    z_j <- Z[idx]
    dt_j <- dt[idx]
    jump_power_j <- mean(dt_j ^ (1 - alpha / 2), na.rm = TRUE)
    w_mat <- matrix(NA_complex_, nrow = block_size, ncol = n_freq)
    for (m in seq_along(freqs)) {
      u <- freqs[m]
      w <- exp(1i * u * z_j)
      w_mat[, m] <- w
      phi <- mean(w)
      phi_hat[j, m] <- phi
      log_mod <- log(max(Mod(phi), min_mod))
      jump_att <- sigma2 ^ alpha * abs(u) ^ alpha * jump_power_j
      y_raw[j, m] <- log_mod
      jump_attenuation[j, m] <- jump_att
      y_dejumped[j, m] <- log_mod + jump_att
    }

    reim <- cbind(Re(w_mat), Im(w_mat))
    cov_reim <- stats::cov(reim, use = "pairwise.complete.obs") / block_size
    grad <- matrix(0, nrow = n_freq, ncol = 2L * n_freq)
    for (m in seq_len(n_freq)) {
      pr <- Re(phi_hat[j, m])
      pi <- Im(phi_hat[j, m])
      den <- max(pr ^ 2 + pi ^ 2, min_mod ^ 2)
      grad[m, m] <- pr / den
      grad[m, n_freq + m] <- pi / den
    }
    cov_y <- grad %*% cov_reim %*% t(grad)
    cov_y <- (cov_y + t(cov_y)) / 2
    for (m in seq_len(n_freq)) {
      pr <- Re(phi_hat[j, m])
      pi <- Im(phi_hat[j, m])
      r2 <- max(pr ^ 2 + pi ^ 2, min_mod ^ 2)
      hess <- matrix(c(
        (r2 - 2 * pr ^ 2) / r2 ^ 2,
        -2 * pr * pi / r2 ^ 2,
        -2 * pr * pi / r2 ^ 2,
        (r2 - 2 * pi ^ 2) / r2 ^ 2
      ), nrow = 2L, byrow = TRUE)
      cov_m <- cov_reim[c(m, n_freq + m), c(m, n_freq + m), drop = FALSE]
      log_mod_bias[j, m] <- 0.5 * sum(hess * cov_m)
    }
    d <- diag(cov_y)
    d[!is.finite(d) | d < min_cov_var] <- min_cov_var
    diag(cov_y) <- d
    obs_cov[j, , ] <- cov_y
  }

  colnames(y_raw) <- colnames(y_dejumped) <- paste0("u", seq_along(freqs))
  colnames(jump_attenuation) <- colnames(log_mod_bias) <-
    paste0("u", seq_along(freqs))
  list(
    y_raw = y_raw,
    y_dejumped = y_dejumped,
    y_dejumped_bias_corrected = y_dejumped - log_mod_bias,
    jump_attenuation = jump_attenuation,
    log_mod_bias = log_mod_bias,
    phi_hat = phi_hat,
    obs_cov = obs_cov,
    block_h = mean(dt) * block_size,
    block_size = block_size,
    n_blocks = n_blocks,
    freqs = freqs,
    scale_z = scale_z
  )
}

select_clock_block_size <- function(R, dt, alpha, sigma2,
                                    candidate_block_sizes = c(500L, 750L,
                                                              1000L),
                                    freq_mults = c(0.50, 0.75, 1.00)) {
  rows <- list()
  k <- 0L
  for (block_size in as.integer(candidate_block_sizes)) {
    activity <- try(block_cf_activity(
      R = R,
      dt = dt,
      alpha = alpha,
      sigma2 = sigma2,
      block_size = block_size,
      freq_mults = freq_mults
    ), silent = TRUE)
    if (inherits(activity, "try-error")) next
    k <- k + 1L
    rows[[k]] <- data.frame(
      block_size = activity$block_size,
      block_h = activity$block_h,
      n_blocks = activity$n_blocks,
      activity_freq_rel_iqr_median = activity$freq_rel_iqr_median,
      generator_h_proxy = activity$block_h,
      transition_se_proxy = 1 / sqrt(activity$n_blocks)
    )
  }
  if (!length(rows)) stop("No candidate block size produced valid activity.")
  tab <- do.call(rbind, rows)
  tab$activity_rank <- rank(tab$activity_freq_rel_iqr_median,
                            ties.method = "average", na.last = "keep")
  tab$generator_rank <- rank(tab$generator_h_proxy,
                             ties.method = "average", na.last = "keep")
  tab$transition_rank <- rank(tab$transition_se_proxy,
                              ties.method = "average", na.last = "keep")
  tab$tradeoff_score <- tab$activity_rank + tab$generator_rank +
    tab$transition_rank
  ord <- order(tab$tradeoff_score, tab$generator_h_proxy,
               tab$activity_freq_rel_iqr_median, na.last = NA)
  selected <- tab[ord[1L], , drop = FALSE]
  list(
    block_size = selected$block_size,
    selected = selected,
    table = tab
  )
}

estimate_beta_from_activity <- function(activity, q_sample,
                                        probs = c(0.25, 0.50, 0.75),
                                        method = c("wls",
                                                   "median_ratio"),
                                        weights = NULL) {
  method <- match.arg(method)
  a <- activity[is.finite(activity) & activity > 0]
  q <- q_sample[is.finite(q_sample) & q_sample > 0]
  if (length(a) < 10L || length(q) < 100L) {
    stop("Not enough activity or stationary samples.")
  }
  aq <- stats::quantile(a, probs = probs, names = FALSE, type = 8,
                        na.rm = TRUE)
  qq <- stats::quantile(q, probs = probs, names = FALSE, type = 8,
                        na.rm = TRUE)
  beta2_vals <- aq / qq
  beta2_vals <- beta2_vals[is.finite(beta2_vals) & beta2_vals > 0]
  if (!length(beta2_vals)) stop("Could not estimate beta from quantiles.")
  if (is.null(weights)) weights <- rep(1, length(probs))
  weights <- as.numeric(weights)
  weights <- weights[seq_along(probs)]
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) weights <- rep(1, length(probs))
  ok_wls <- is.finite(aq) & is.finite(qq) & qq > 0 &
    is.finite(weights) & weights > 0
  beta2_wls <- if (any(ok_wls)) {
    sum(weights[ok_wls] * qq[ok_wls] * aq[ok_wls]) /
      sum(weights[ok_wls] * qq[ok_wls] ^ 2)
  } else {
    NA_real_
  }
  beta2_median <- stats::median(beta2_vals)
  beta2_hat <- if (identical(method, "wls") &&
                   is.finite(beta2_wls) && beta2_wls > 0) {
    beta2_wls
  } else {
    beta2_median
  }
  list(
    beta_hat = sqrt(beta2_hat),
    beta2_hat = beta2_hat,
    beta_method = method,
    beta2_wls = beta2_wls,
    beta2_median_ratio = beta2_median,
    beta2_values = beta2_vals,
    probs = probs,
    activity_quantiles = aq,
    q_quantiles = qq
  )
}

estimate_beta_from_activity_matrix <- function(activity_by_freq, q_sample,
                                               probs = c(0.25, 0.50,
                                                         0.75),
                                               min_positive = 10L) {
  a_mat <- as.matrix(activity_by_freq)
  q <- q_sample[is.finite(q_sample) & q_sample > 0]
  if (length(q) < 100L) stop("Not enough stationary Q samples.")
  qq <- stats::quantile(q, probs = probs, names = FALSE, type = 8,
                        na.rm = TRUE)
  rows <- list()
  k <- 0L
  for (m in seq_len(ncol(a_mat))) {
    a <- a_mat[, m]
    a <- a[is.finite(a) & a > 0]
    if (length(a) < min_positive) next
    aq <- stats::quantile(a, probs = probs, names = FALSE, type = 8,
                          na.rm = TRUE)
    ok <- is.finite(aq) & is.finite(qq) & qq > 0
    if (!any(ok)) next
    k <- k + 1L
    rows[[k]] <- data.frame(
      frequency_index = m,
      prob = probs[ok],
      activity_quantile = aq[ok],
      q_quantile = qq[ok]
    )
  }
  tab <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(tab)) stop("No valid frequency-level activity quantiles.")
  beta2_hat <- sum(tab$q_quantile * tab$activity_quantile) /
    sum(tab$q_quantile ^ 2)
  if (!is.finite(beta2_hat) || beta2_hat <= 0) {
    stop("Could not estimate beta from activity matrix.")
  }
  tab$beta2_component <- tab$activity_quantile / tab$q_quantile
  list(
    beta_hat = sqrt(beta2_hat),
    beta2_hat = beta2_hat,
    beta_method = "frequency_vector_wls",
    frequency_quantiles = tab,
    probs = probs,
    q_quantiles = qq,
    n_frequencies = length(unique(tab$frequency_index))
  )
}

activity_regime_stats <- function(activity, beta_hat, q_thresholds,
                                  lags = c(1L, 2L, 4L)) {
  q_hat <- activity / beta_hat ^ 2
  keep <- is.finite(q_hat)
  q_hat <- q_hat[keep]
  regimes <- cut(
    q_hat,
    breaks = c(-Inf, q_thresholds, Inf),
    labels = FALSE,
    right = TRUE
  )
  n_reg <- length(q_thresholds) + 1L
  occupation <- tabulate(regimes, nbins = n_reg) / length(regimes)
  transitions <- list()
  flat <- occupation
  names(flat) <- paste0("occ_", seq_len(n_reg))
  for (lag in as.integer(lags)) {
    if (!is.finite(lag) || lag < 1L || length(regimes) <= lag) next
    from <- regimes[seq_len(length(regimes) - lag)]
    to <- regimes[(lag + 1L):length(regimes)]
    mat <- matrix(0, n_reg, n_reg)
    for (i in seq_along(from)) {
      if (is.finite(from[i]) && is.finite(to[i])) {
        mat[from[i], to[i]] <- mat[from[i], to[i]] + 1
      }
    }
    denom <- rowSums(mat)
    prob <- mat
    for (a in seq_len(n_reg)) {
      if (denom[a] > 0) prob[a, ] <- mat[a, ] / denom[a]
    }
    transitions[[paste0("lag", lag)]] <- prob
    vals <- as.numeric(prob)
    names(vals) <- paste0("lag", lag, "_p", rep(seq_len(n_reg), each = n_reg),
                          "_", rep(seq_len(n_reg), times = n_reg))
    flat <- c(flat, vals)
  }
  list(
    q_hat = q_hat,
    regimes = regimes,
    occupation = occupation,
    transitions = transitions,
    summary = flat
  )
}

generator_test_eval <- function(q, test) {
  q <- pmax(as.numeric(q), 1 + 1e-10)
  if (identical(test, "log")) {
    f <- log(q)
    a1 <- 4 / q - 3
  } else if (identical(test, "inv1")) {
    c0 <- 1
    f <- 1 / (q + c0)
    fp <- -1 / (q + c0) ^ 2
    fpp <- 2 / (q + c0) ^ 3
    a1 <- (2 - q) * fp + 2 * q * (q - 1) * fpp
  } else if (identical(test, "inv4")) {
    c0 <- 4
    f <- 1 / (q + c0)
    fp <- -1 / (q + c0) ^ 2
    fpp <- 2 / (q + c0) ^ 3
    a1 <- (2 - q) * fp + 2 * q * (q - 1) * fpp
  } else if (identical(test, "ratio1")) {
    c0 <- 1
    f <- q / (q + c0)
    fp <- c0 / (q + c0) ^ 2
    fpp <- -2 * c0 / (q + c0) ^ 3
    a1 <- (2 - q) * fp + 2 * q * (q - 1) * fpp
  } else if (identical(test, "exp025")) {
    s <- 0.25
    f <- exp(-s * q)
    fp <- -s * f
    fpp <- s ^ 2 * f
    a1 <- (2 - q) * fp + 2 * q * (q - 1) * fpp
  } else {
    stop("Unknown generator test function: ", test)
  }
  list(f = f, a1 = a1)
}

generator_instrument_eval <- function(q, instrument, q_thresholds = NULL) {
  q <- pmax(as.numeric(q), 1 + 1e-10)
  center <- function(x) x - mean(x[is.finite(x)], na.rm = TRUE)
  if (identical(instrument, "center_log")) {
    w <- center(log(q))
  } else if (identical(instrument, "bounded_poly1")) {
    x <- q / (q + 1)
    w <- center(x)
  } else if (identical(instrument, "bounded_poly2")) {
    x <- q / (q + 1)
    xc <- center(x)
    w <- center(xc ^ 2)
  } else if (identical(instrument, "tanh_log")) {
    x <- tanh(log(q) - stats::median(log(q), na.rm = TRUE))
    w <- center(x)
  } else if (identical(instrument, "logistic_low")) {
    thr <- if (length(q_thresholds)) q_thresholds[1L] else
      stats::quantile(q, 0.40, names = FALSE, na.rm = TRUE)
    s <- stats::mad(log(q), constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- 1
    x <- 1 / (1 + exp((log(q) - log(thr)) / s))
    w <- center(x)
  } else if (identical(instrument, "logistic_high")) {
    thr <- if (length(q_thresholds) >= 2L) q_thresholds[2L] else
      stats::quantile(q, 0.75, names = FALSE, na.rm = TRUE)
    s <- stats::mad(log(q), constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- 1
    x <- 1 / (1 + exp(-(log(q) - log(thr)) / s))
    w <- center(x)
  } else if (identical(instrument, "low_bin")) {
    thr <- if (length(q_thresholds)) q_thresholds[1L] else
      stats::quantile(q, 0.40, names = FALSE, na.rm = TRUE)
    w <- center(as.numeric(q <= thr))
  } else if (identical(instrument, "mid_bin")) {
    if (length(q_thresholds) < 2L) {
      q_thresholds <- stats::quantile(q, c(0.40, 0.75), names = FALSE,
                                      na.rm = TRUE)
    }
    w <- center(as.numeric(q > q_thresholds[1L] & q <= q_thresholds[2L]))
  } else if (identical(instrument, "high_bin")) {
    thr <- if (length(q_thresholds) >= 2L) q_thresholds[2L] else
      stats::quantile(q, 0.75, names = FALSE, na.rm = TRUE)
    w <- center(as.numeric(q > thr))
  } else {
    stop("Unknown generator instrument: ", instrument)
  }
  w
}

estimate_rho_generator_gmm <- function(q_hat, block_h,
                                       tests = c("log", "inv1",
                                                 "ratio1", "exp025"),
                                       instruments = "a1",
                                       q_thresholds = NULL,
                                       lags = c(1L, 2L, 4L),
                                       standardize = FALSE,
                                       rho_grid = exp(seq(log(0.40),
                                                          log(3.50),
                                                          length.out = 31L))) {
  q_hat <- q_hat[is.finite(q_hat) & q_hat >= 1]
  if (length(q_hat) < 20L) stop("Not enough normalized activity values.")
  rows <- list()
  k <- 0L
  for (test in tests) {
    ev <- generator_test_eval(q_hat, test)
    for (instrument in instruments) {
      w_all <- if (identical(instrument, "a1")) {
        ev$a1
      } else {
        generator_instrument_eval(q_hat, instrument,
                                  q_thresholds = q_thresholds)
      }
      for (lag in as.integer(lags)) {
        if (!is.finite(lag) || lag < 1L || length(q_hat) <= lag) next
        idx <- seq_len(length(q_hat) - lag)
        df <- ev$f[idx + lag] - ev$f[idx]
        a1 <- ev$a1[idx]
        w <- w_all[idx]
        good <- is.finite(df) & is.finite(a1) & is.finite(w)
        if (!any(good)) next
        h <- lag * block_h
        c_moment <- mean(w[good] * df[good])
        b_moment <- h * mean(w[good] * a1[good])
        moment_scale <- sqrt(mean((w[good] * df[good]) ^ 2))
        if (!is.finite(moment_scale) || moment_scale <= 0) {
          moment_scale <- sqrt(mean((h * w[good] * a1[good]) ^ 2))
        }
        if (!is.finite(c_moment) || !is.finite(b_moment) ||
            !is.finite(moment_scale) || moment_scale <= 0) next
        if (abs(b_moment) <= .Machine$double.eps) next
        k <- k + 1L
        rows[[k]] <- data.frame(
          test = test,
          instrument = instrument,
          lag = lag,
          h = h,
          c_moment = c_moment,
          b_moment = b_moment,
          moment_scale = moment_scale,
          c_scaled = c_moment / moment_scale,
          b_scaled = b_moment / moment_scale,
          theta_component = c_moment / b_moment,
          sensitivity_abs = abs(b_moment)
        )
      }
    }
  }
  moments <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(moments)) stop("No valid generator-GMM moments.")
  c_use <- if (isTRUE(standardize)) moments$c_scaled else moments$c_moment
  b_use <- if (isTRUE(standardize)) moments$b_scaled else moments$b_moment
  theta_unconstrained <- sum(b_use * c_use) / sum(b_use ^ 2)
  curve <- data.frame(rho = rho_grid)
  curve$theta <- curve$rho ^ 2
  curve$objective <- vapply(curve$theta, function(theta) {
    mean((c_use - theta * b_use) ^ 2)
  }, numeric(1))
  best <- curve[which.min(curve$objective), , drop = FALSE]
  rho_closed <- if (is.finite(theta_unconstrained) &&
                    theta_unconstrained > 0) {
    sqrt(theta_unconstrained)
  } else {
    NA_real_
  }
  rho_hat <- if (is.finite(rho_closed) &&
                 rho_closed >= min(rho_grid) &&
                 rho_closed <= max(rho_grid)) rho_closed else best$rho
  theta_hat <- rho_hat ^ 2
  list(
    rho_hat = rho_hat,
    theta_hat = theta_hat,
    theta_unconstrained = theta_unconstrained,
    objective = if (is.finite(rho_closed)) {
      mean((c_use - theta_hat * b_use) ^ 2)
    } else {
      best$objective
    },
    curve = curve,
    moments = moments,
    tests = paste(tests, collapse = ","),
    instruments = paste(instruments, collapse = ","),
    lags = paste(lags, collapse = ","),
    standardize = isTRUE(standardize),
    theta_curvature = 2 * mean(b_use ^ 2)
  )
}

estimate_rho_vector_ecf_generator_gmm <- function(activity_by_freq, beta_hat,
                                                  block_h,
                                                  tests = c("log", "inv1",
                                                            "inv4",
                                                            "ratio1",
                                                            "exp025"),
                                                  instruments = c(
                                                    "a1", "center_log",
                                                    "bounded_poly1",
                                                    "bounded_poly2",
                                                    "tanh_log",
                                                    "logistic_low",
                                                    "logistic_high",
                                                    "low_bin", "mid_bin",
                                                    "high_bin"
                                                  ),
                                                  q_thresholds = NULL,
                                                  lags = c(1L, 2L),
                                                  standardize = TRUE,
                                                  rho_grid = exp(seq(
                                                    log(0.40), log(3.50),
                                                    length.out = 31L))) {
  a_mat <- as.matrix(activity_by_freq)
  if (!is.finite(beta_hat) || beta_hat <= 0) stop("Invalid beta_hat.")
  rows <- list()
  k <- 0L
  for (m in seq_len(ncol(a_mat))) {
    q_m <- a_mat[, m] / beta_hat ^ 2
    q_m <- q_m[is.finite(q_m) & q_m > 0]
    if (length(q_m) < 20L) next
    q_m <- pmax(q_m, 1 + 1e-10)
    for (test in tests) {
      ev <- generator_test_eval(q_m, test)
      for (instrument in instruments) {
        w_all <- if (identical(instrument, "a1")) {
          ev$a1
        } else {
          generator_instrument_eval(q_m, instrument,
                                    q_thresholds = q_thresholds)
        }
        for (lag in as.integer(lags)) {
          if (!is.finite(lag) || lag < 1L || length(q_m) <= lag) next
          idx <- seq_len(length(q_m) - lag)
          df <- ev$f[idx + lag] - ev$f[idx]
          a1 <- ev$a1[idx]
          w <- w_all[idx]
          good <- is.finite(df) & is.finite(a1) & is.finite(w)
          if (!any(good)) next
          h <- lag * block_h
          c_moment <- mean(w[good] * df[good])
          b_moment <- h * mean(w[good] * a1[good])
          moment_scale <- sqrt(mean((w[good] * df[good]) ^ 2))
          if (!is.finite(moment_scale) || moment_scale <= 0) {
            moment_scale <- sqrt(mean((h * w[good] * a1[good]) ^ 2))
          }
          if (!is.finite(c_moment) || !is.finite(b_moment) ||
              !is.finite(moment_scale) || moment_scale <= 0) next
          if (abs(b_moment) <= .Machine$double.eps) next
          k <- k + 1L
          rows[[k]] <- data.frame(
            frequency_index = m,
            test = test,
            instrument = instrument,
            lag = lag,
            h = h,
            c_moment = c_moment,
            b_moment = b_moment,
            moment_scale = moment_scale,
            c_scaled = c_moment / moment_scale,
            b_scaled = b_moment / moment_scale,
            theta_component = c_moment / b_moment,
            sensitivity_abs = abs(b_moment)
          )
        }
      }
    }
  }
  moments <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(moments)) stop("No valid vector empirical-CF GMM moments.")
  c_use <- if (isTRUE(standardize)) moments$c_scaled else moments$c_moment
  b_use <- if (isTRUE(standardize)) moments$b_scaled else moments$b_moment
  theta_unconstrained <- sum(b_use * c_use) / sum(b_use ^ 2)
  curve <- data.frame(rho = rho_grid)
  curve$theta <- curve$rho ^ 2
  curve$objective <- vapply(curve$theta, function(theta) {
    mean((c_use - theta * b_use) ^ 2)
  }, numeric(1))
  best <- curve[which.min(curve$objective), , drop = FALSE]
  rho_closed <- if (is.finite(theta_unconstrained) &&
                    theta_unconstrained > 0) {
    sqrt(theta_unconstrained)
  } else {
    NA_real_
  }
  rho_hat <- if (is.finite(rho_closed) &&
                 rho_closed >= min(rho_grid) &&
                 rho_closed <= max(rho_grid)) rho_closed else best$rho
  theta_hat <- rho_hat ^ 2
  list(
    rho_hat = rho_hat,
    theta_hat = theta_hat,
    theta_unconstrained = theta_unconstrained,
    objective = mean((c_use - theta_hat * b_use) ^ 2),
    curve = curve,
    moments = moments,
    tests = paste(tests, collapse = ","),
    instruments = paste(instruments, collapse = ","),
    lags = paste(lags, collapse = ","),
    standardize = isTRUE(standardize),
    theta_curvature = 2 * mean(b_use ^ 2),
    method = "vector_ecf_generator_gmm"
  )
}

estimate_sigma1_rho_vector_ecf_clock <- function(R, dt, alpha, sigma2,
                                                 block_size,
                                                 rho_grid = exp(seq(
                                                   log(0.40), log(3.50),
                                                   length.out = 31L)),
                                                 freq_mults = c(0.50, 0.75,
                                                                1.00),
                                                 q_probs = c(0.25, 0.50,
                                                             0.75),
                                                 generator_tests = c(
                                                   "log", "inv1", "inv4",
                                                   "ratio1", "exp025"
                                                 ),
                                                 generator_instruments = c(
                                                   "a1", "center_log",
                                                   "bounded_poly1",
                                                   "bounded_poly2",
                                                   "tanh_log",
                                                   "logistic_low",
                                                   "logistic_high",
                                                   "low_bin", "mid_bin",
                                                   "high_bin"
                                                 ),
                                                 regime_probs = c(0.40,
                                                                  0.75),
                                                 lags = c(1L, 2L),
                                                 q_sample_n = 200000L,
                                                 seed = 1L) {
  q_sample <- stationary_q_sample(q_sample_n, seed = seed)
  activity <- block_cf_activity(
    R = R,
    dt = dt,
    alpha = alpha,
    sigma2 = sigma2,
    block_size = block_size,
    freq_mults = freq_mults
  )
  beta <- estimate_beta_from_activity_matrix(
    activity_by_freq = activity$activity_by_freq,
    q_sample = q_sample,
    probs = q_probs
  )
  q_thresholds <- stats::quantile(q_sample, probs = regime_probs,
                                  names = FALSE, type = 8)
  rho <- estimate_rho_vector_ecf_generator_gmm(
    activity_by_freq = activity$activity_by_freq,
    beta_hat = beta$beta_hat,
    block_h = activity$block_h,
    tests = generator_tests,
    instruments = generator_instruments,
    q_thresholds = q_thresholds,
    lags = lags,
    standardize = TRUE,
    rho_grid = rho_grid
  )
  list(
    sigma1_hat = beta$beta_hat / rho$rho_hat,
    rho_hat = rho$rho_hat,
    beta_hat = beta$beta_hat,
    beta_details = beta,
    rho_curve = rho$curve,
    rho_moments = rho$moments,
    activity = activity,
    q_thresholds = q_thresholds,
    block_h = activity$block_h,
    block_size = activity$block_size,
    n_blocks = activity$n_blocks,
    method = "vector_empirical_cf_generator_clock"
  )
}

transition_moment_vector <- function(q_hat,
                                     tests = c("log", "inv1", "inv4",
                                               "ratio1", "exp025"),
                                     instruments = c(
                                       "a1", "center_log",
                                       "bounded_poly1", "bounded_poly2",
                                       "tanh_log", "logistic_low",
                                       "logistic_high", "low_bin",
                                       "mid_bin", "high_bin"
                                     ),
                                     q_thresholds = NULL,
                                     lags = c(1L, 2L)) {
  q_hat <- q_hat[is.finite(q_hat) & q_hat > 0]
  if (length(q_hat) < 20L) stop("Not enough activity values.")
  rows <- list()
  k <- 0L
  for (test in tests) {
    ev <- generator_test_eval(q_hat, test)
    for (instrument in instruments) {
      w_all <- if (identical(instrument, "a1")) {
        ev$a1
      } else {
        generator_instrument_eval(q_hat, instrument,
                                  q_thresholds = q_thresholds)
      }
      for (lag in as.integer(lags)) {
        if (!is.finite(lag) || lag < 1L || length(q_hat) <= lag) next
        idx <- seq_len(length(q_hat) - lag)
        y <- w_all[idx] * (ev$f[idx + lag] - ev$f[idx])
        y <- y[is.finite(y)]
        if (!length(y)) next
        scale <- sqrt(mean(y ^ 2))
        if (!is.finite(scale) || scale <= 0) scale <- 1
        k <- k + 1L
        rows[[k]] <- data.frame(
          test = test,
          instrument = instrument,
          lag = lag,
          moment = mean(y),
          scale = scale,
          n_pairs = length(y)
        )
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(out)) stop("No valid transition moments.")
  out$key <- paste(out$test, out$instrument, out$lag, sep = "|")
  out
}

estimate_activity_log_noise_model <- function(activity_obj,
                                              truth_activity = NULL,
                                              q_ref = NULL,
                                              mode = c("replicate_median",
                                                       "oracle_log_error",
                                                       "none"),
                                              center = TRUE) {
  mode <- match.arg(mode)
  if (identical(mode, "none")) {
    return(list(mode = "none", errors = 0, n_freq = 1L, bias = 0,
                sd = 0, mad = 0, n_error = 1L))
  }
  activity <- activity_obj$activity
  by_freq <- activity_obj$activity_by_freq
  if (identical(mode, "oracle_log_error")) {
    if (is.null(truth_activity)) {
      stop("truth_activity is required for oracle_log_error mode.")
    }
    ok <- is.finite(activity) & activity > 0 &
      is.finite(truth_activity) & truth_activity > 0
    err <- log(activity[ok] / truth_activity[ok])
    if (isTRUE(center)) err <- err - stats::median(err, na.rm = TRUE)
    err <- err[is.finite(err)]
    if (!length(err)) stop("No valid oracle log errors.")
    return(list(
      mode = "oracle_log_error",
      errors = err,
      q_ref = if (is.null(q_ref)) NULL else q_ref[ok][is.finite(err)],
      n_freq = 1L,
      bias = mean(err, na.rm = TRUE),
      median = stats::median(err, na.rm = TRUE),
      sd = stats::sd(err, na.rm = TRUE),
      mad = stats::mad(err, constant = 1.4826, na.rm = TRUE),
      n_error = length(err)
    ))
  }
  if (is.null(by_freq)) {
    stop("activity_by_freq is required for replicate_median mode.")
  }
  if (is.null(activity)) {
    activity <- apply(by_freq, 1L, function(x) {
      x <- x[is.finite(x) & x > 0]
      if (length(x)) stats::median(x) else NA_real_
    })
  }
  ratio <- sweep(by_freq, 1L, activity, "/")
  err <- as.numeric(log(ratio))
  err <- err[is.finite(err)]
  if (isTRUE(center)) err <- err - stats::median(err, na.rm = TRUE)
  err <- err[is.finite(err)]
  if (!length(err)) stop("No valid replicate log errors.")
  list(
    mode = "replicate_median",
    errors = err,
    n_freq = max(1L, ncol(by_freq)),
    bias = mean(err, na.rm = TRUE),
    median = stats::median(err, na.rm = TRUE),
    sd = stats::sd(err, na.rm = TRUE),
    mad = stats::mad(err, constant = 1.4826, na.rm = TRUE),
    n_error = length(err)
  )
}

estimate_activity_bootstrap_log_noise_model <- function(R, dt, alpha, sigma2,
                                                        activity_obj,
                                                        q_ref = NULL,
                                                        n_boot = 20L,
                                                        min_mod = 1e-8,
                                                        seed = 1L) {
  set.seed(seed)
  n_boot <- as.integer(n_boot)
  if (n_boot < 1L) stop("n_boot must be positive.")
  block_size <- as.integer(activity_obj$block_size)
  n_blocks <- as.integer(activity_obj$n_blocks)
  freqs <- activity_obj$freqs
  if (!length(freqs)) stop("activity_obj has no frequencies.")
  n <- n_blocks * block_size
  R <- R[seq_len(n)]
  dt <- dt[seq_len(n)]
  Z <- R / sqrt(dt)
  errors <- numeric(0)
  error_q_ref <- numeric(0)
  for (j in seq_len(n_blocks)) {
    idx <- ((j - 1L) * block_size + 1L):(j * block_size)
    z_j <- Z[idx]
    dt_j <- dt[idx]
    base <- activity_obj$activity[j]
    if (!is.finite(base) || base <= 0) next
    for (b in seq_len(n_boot)) {
      boot_idx <- sample.int(block_size, size = block_size, replace = TRUE)
      z_b <- z_j[boot_idx]
      dt_b <- dt_j[boot_idx]
      jump_power_b <- mean(dt_b ^ (1 - alpha / 2), na.rm = TRUE)
      vals <- rep(NA_real_, length(freqs))
      for (m in seq_along(freqs)) {
        u <- freqs[m]
        phi <- mean(exp(1i * u * z_b))
        log_mod <- log(max(Mod(phi), min_mod))
        jump_attenuation <- sigma2 ^ alpha * abs(u) ^ alpha * jump_power_b
        vals[m] <- -2 * (log_mod + jump_attenuation) / (u ^ 2)
      }
      vals <- vals[is.finite(vals) & vals > 0]
      if (length(vals)) {
        boot_activity <- stats::median(vals)
        if (is.finite(boot_activity) && boot_activity > 0) {
          errors <- c(errors, log(boot_activity / base))
          if (!is.null(q_ref) && length(q_ref) >= j) {
            error_q_ref <- c(error_q_ref, q_ref[j])
          }
        }
      }
    }
  }
  errors <- errors[is.finite(errors)]
  if (!length(errors)) stop("No valid bootstrap activity errors.")
  errors <- errors - stats::median(errors, na.rm = TRUE)
  list(
    mode = "bootstrap_log_error",
    errors = errors,
    q_ref = if (length(error_q_ref) == length(errors)) error_q_ref else NULL,
    n_freq = 1L,
    n_boot = n_boot,
    bias = mean(errors, na.rm = TRUE),
    median = stats::median(errors, na.rm = TRUE),
    sd = stats::sd(errors, na.rm = TRUE),
    mad = stats::mad(errors, constant = 1.4826, na.rm = TRUE),
    n_error = length(errors)
  )
}

condition_activity_noise_model <- function(noise_model,
                                           probs = c(1 / 3, 2 / 3)) {
  q_ref <- noise_model$q_ref
  err <- noise_model$errors
  ok <- is.finite(q_ref) & is.finite(err)
  q_ref <- q_ref[ok]
  err <- err[ok]
  if (length(err) < 10L) {
    stop("Not enough errors with q_ref for conditional noise model.")
  }
  breaks <- stats::quantile(q_ref, probs = probs, names = FALSE,
                            type = 8, na.rm = TRUE)
  bin <- cut(q_ref, breaks = c(-Inf, breaks, Inf), labels = FALSE,
             right = TRUE)
  errors_by_bin <- split(err, bin)
  errors_by_bin <- lapply(errors_by_bin, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) x else err
  })
  noise_model$conditional <- TRUE
  noise_model$breaks <- breaks
  noise_model$errors_by_bin <- errors_by_bin
  noise_model$mode <- paste0(noise_model$mode, "_conditional")
  noise_model
}

draw_activity_log_noise <- function(noise_model, n, q_ref = NULL) {
  n <- as.integer(n)
  if (is.null(noise_model) || identical(noise_model$mode, "none")) {
    return(rep(0, n))
  }
  err <- noise_model$errors
  err <- err[is.finite(err)]
  if (!length(err)) return(rep(0, n))
  if (isTRUE(noise_model$conditional) && !is.null(q_ref)) {
    bin <- cut(q_ref, breaks = c(-Inf, noise_model$breaks, Inf),
               labels = FALSE, right = TRUE)
    out <- numeric(n)
    for (i in seq_len(n)) {
      pool <- noise_model$errors_by_bin[[as.character(bin[i])]]
      if (is.null(pool) || !length(pool)) pool <- err
      out[i] <- sample(pool, size = 1L, replace = TRUE)
    }
    return(out)
  }
  if (identical(noise_model$mode, "replicate_median")) {
    n_freq <- max(1L, as.integer(noise_model$n_freq))
    draws <- sample(err, size = n * n_freq, replace = TRUE)
    draw_mat <- matrix(draws, nrow = n, ncol = n_freq)
    apply(draw_mat, 1L, stats::median)
  } else {
    sample(err, size = n, replace = TRUE)
  }
}

simulate_block_average_q <- function(rho, block_h, n_blocks,
                                     substeps_per_block = 20L,
                                     burn_blocks = 1000L,
                                     seed = 1L) {
  set.seed(seed)
  n_blocks <- as.integer(n_blocks)
  substeps_per_block <- as.integer(substeps_per_block)
  if (n_blocks < 10L) stop("n_blocks must be at least 10.")
  if (substeps_per_block < 1L) stop("substeps_per_block must be positive.")
  dt_sub <- block_h / substeps_per_block
  v <- stationary_v_sample(1L)
  rho2 <- rho ^ 2
  total_blocks <- n_blocks + as.integer(burn_blocks)
  q_bar <- numeric(n_blocks)
  out_i <- 0L
  for (b in seq_len(total_blocks)) {
    q_sum <- 0
    for (s in seq_len(substeps_per_block)) {
      q_now <- 1 + v ^ 2
      q_sum <- q_sum + q_now
      d_v <- -rho2 * v * dt_sub +
        rho * sqrt(pmax(q_now, 1e-12)) * sqrt(dt_sub) *
        stats::rnorm(1L)
      v <- v + d_v
    }
    if (b > burn_blocks) {
      out_i <- out_i + 1L
      q_bar[out_i] <- q_sum / substeps_per_block
    }
  }
  q_bar
}

log_sum_exp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

regularized_chol <- function(S, min_var = 1e-10) {
  S <- as.matrix(S)
  S <- (S + t(S)) / 2
  d <- diag(S)
  d[!is.finite(d) | d < min_var] <- min_var
  diag(S) <- d
  base <- stats::median(d[is.finite(d) & d > 0], na.rm = TRUE)
  if (!is.finite(base) || base <= 0) base <- min_var
  for (ridge_mult in c(0, 1e-8, 1e-6, 1e-4, 1e-2)) {
    out <- try(chol(S + diag(base * ridge_mult, nrow(S))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  chol(diag(pmax(d, min_var), nrow(S)))
}

build_qbar_hmm_transition <- function(rho, block_h, n_states = 31L,
                                      n_sim_blocks = 10000L,
                                      substeps_per_block = 10L,
                                      burn_blocks = 1000L,
                                      transition_smoothing = 1e-4,
                                      seed = 1L) {
  n_states <- as.integer(n_states)
  n_sim_blocks <- as.integer(n_sim_blocks)
  if (n_states < 5L) stop("n_states must be at least 5.")
  if (n_sim_blocks < 10L * n_states) {
    stop("n_sim_blocks is too small for the requested state grid.")
  }
  q_sim <- simulate_block_average_q(
    rho = rho,
    block_h = block_h,
    n_blocks = n_sim_blocks + 1L,
    substeps_per_block = substeps_per_block,
    burn_blocks = burn_blocks,
    seed = seed
  )
  probs <- seq(0, 1, length.out = n_states + 1L)
  breaks <- stats::quantile(q_sim, probs = probs, names = FALSE, type = 8,
                            na.rm = TRUE)
  breaks <- cummax(breaks)
  for (i in 2:length(breaks)) {
    if (breaks[i] <= breaks[i - 1L]) {
      breaks[i] <- breaks[i - 1L] + .Machine$double.eps
    }
  }
  cut_breaks <- breaks
  cut_breaks[1L] <- -Inf
  cut_breaks[length(cut_breaks)] <- Inf
  state <- cut(q_sim, breaks = cut_breaks, labels = FALSE, right = TRUE)
  centers <- numeric(n_states)
  for (s in seq_len(n_states)) {
    vals <- q_sim[state == s]
    centers[s] <- if (length(vals)) {
      stats::median(vals, na.rm = TRUE)
    } else {
      stats::median(breaks[c(s, s + 1L)])
    }
  }
  centers <- pmax(centers, 1 + 1e-10)
  init_counts <- tabulate(state[-length(state)], nbins = n_states) +
    transition_smoothing
  init_prob <- init_counts / sum(init_counts)
  trans_counts <- matrix(transition_smoothing, n_states, n_states)
  from <- state[-length(state)]
  to <- state[-1L]
  good <- is.finite(from) & is.finite(to)
  for (i in which(good)) {
    trans_counts[from[i], to[i]] <- trans_counts[from[i], to[i]] + 1
  }
  transition <- trans_counts / rowSums(trans_counts)
  list(
    state_centers = centers,
    transition = transition,
    init_prob = init_prob,
    breaks = breaks,
    q_sim = q_sim,
    rho = rho,
    block_h = block_h,
    n_states = n_states,
    n_sim_blocks = n_sim_blocks,
    substeps_per_block = substeps_per_block
  )
}

ecf_hmm_emission_log_matrix <- function(y_dejumped, obs_cov, state_centers,
                                        beta, freqs,
                                        cov_inflation = 1,
                                        min_cov_var = 1e-10) {
  y <- as.matrix(y_dejumped)
  n_blocks <- nrow(y)
  n_freq <- ncol(y)
  n_states <- length(state_centers)
  if (dim(obs_cov)[1L] != n_blocks || dim(obs_cov)[2L] != n_freq ||
      dim(obs_cov)[3L] != n_freq) {
    stop("obs_cov dimensions do not match y_dejumped.")
  }
  mean_by_state <- outer(-0.5 * freqs ^ 2 * beta ^ 2, state_centers)
  out <- matrix(NA_real_, nrow = n_blocks, ncol = n_states)
  log_2pi <- log(2 * pi)
  for (j in seq_len(n_blocks)) {
    y_j <- y[j, ]
    if (any(!is.finite(y_j))) next
    S <- obs_cov[j, , ] * cov_inflation
    R <- regularized_chol(S, min_var = min_cov_var)
    diff_mat <- sweep(mean_by_state, 1L, y_j, "-")
    z <- backsolve(R, diff_mat, transpose = TRUE)
    qf <- colSums(z ^ 2)
    logdet <- 2 * sum(log(diag(R)))
    out[j, ] <- -0.5 * (n_freq * log_2pi + logdet + qf)
  }
  out
}

hmm_forward_loglik <- function(log_emission, transition, init_prob) {
  log_emission <- as.matrix(log_emission)
  n_blocks <- nrow(log_emission)
  n_states <- ncol(log_emission)
  if (length(init_prob) != n_states ||
      any(dim(transition) != c(n_states, n_states))) {
    stop("HMM transition dimensions do not match emissions.")
  }
  log_transition <- log(pmax(transition, .Machine$double.xmin))
  log_alpha <- log(pmax(init_prob, .Machine$double.xmin)) +
    log_emission[1L, ]
  inc <- log_sum_exp(log_alpha)
  if (!is.finite(inc)) return(NA_real_)
  loglik <- inc
  log_alpha <- log_alpha - inc
  if (n_blocks >= 2L) {
    for (j in 2:n_blocks) {
      pred <- vapply(seq_len(n_states), function(s) {
        log_sum_exp(log_alpha + log_transition[, s])
      }, numeric(1))
      log_alpha <- pred + log_emission[j, ]
      inc <- log_sum_exp(log_alpha)
      if (!is.finite(inc)) return(NA_real_)
      loglik <- loglik + inc
      log_alpha <- log_alpha - inc
    }
  }
  loglik
}

estimate_beta_rho_ecf_hmm <- function(ecf_obs,
                                      beta_grid = exp(seq(log(0.50),
                                                          log(3.50),
                                                          length.out = 15L)),
                                      rho_grid = exp(seq(log(0.40),
                                                         log(3.50),
                                                         length.out = 17L)),
                                      n_states = 31L,
                                      qbar_sim_blocks = 10000L,
                                      qbar_substeps_per_block = 10L,
                                      burn_blocks = 1000L,
                                      cov_inflation = 1,
                                      min_cov_var = 1e-10,
                                      transition_models = NULL,
                                      seed = 1L) {
  y <- as.matrix(ecf_obs$y_dejumped)
  keep <- is.finite(rowSums(y))
  if (!any(keep)) stop("No finite ECF observation rows.")
  y <- y[keep, , drop = FALSE]
  obs_cov <- ecf_obs$obs_cov[keep, , , drop = FALSE]
  freqs <- ecf_obs$freqs
  block_h <- ecf_obs$block_h
  rows <- vector("list", length(rho_grid) * length(beta_grid))
  trans_store <- vector("list", length(rho_grid))
  k <- 0L
  for (i in seq_along(rho_grid)) {
    rho <- rho_grid[i]
    if (!is.null(transition_models)) {
      q_hmm <- transition_models[[i]]
      if (is.null(q_hmm)) stop("transition_models has a missing entry.")
    } else {
      q_hmm <- build_qbar_hmm_transition(
        rho = rho,
        block_h = block_h,
        n_states = n_states,
        n_sim_blocks = qbar_sim_blocks,
        substeps_per_block = qbar_substeps_per_block,
        burn_blocks = burn_blocks,
        seed = seed + 7919L * i
      )
    }
    trans_store[[i]] <- q_hmm
    for (beta in beta_grid) {
      log_emission <- ecf_hmm_emission_log_matrix(
        y_dejumped = y,
        obs_cov = obs_cov,
        state_centers = q_hmm$state_centers,
        beta = beta,
        freqs = freqs,
        cov_inflation = cov_inflation,
        min_cov_var = min_cov_var
      )
      loglik <- hmm_forward_loglik(
        log_emission = log_emission,
        transition = q_hmm$transition,
        init_prob = q_hmm$init_prob
      )
      k <- k + 1L
      rows[[k]] <- data.frame(
        beta = beta,
        rho = rho,
        sigma1 = beta / rho,
        loglik = loglik,
        n_blocks = nrow(y),
        n_states = n_states,
        qbar_sim_blocks = qbar_sim_blocks
      )
    }
  }
  grid <- do.call(rbind, rows)
  grid <- grid[is.finite(grid$loglik), , drop = FALSE]
  if (!nrow(grid)) stop("No finite ECF-HMM likelihood evaluations.")
  best <- grid[which.max(grid$loglik), , drop = FALSE]
  profile <- aggregate(loglik ~ rho, data = grid, FUN = max)
  profile$deviance <- -2 * (profile$loglik - max(profile$loglik))
  support <- profile[profile$deviance <= stats::qchisq(0.95, df = 1), ,
                     drop = FALSE]
  list(
    beta_hat = best$beta,
    rho_hat = best$rho,
    sigma1_hat = best$sigma1,
    loglik = best$loglik,
    grid = grid,
    rho_profile = profile,
    rho_support = support,
    rho_support_range = range(support$rho, finite = TRUE),
    transition_models = trans_store,
    block_h = block_h,
    n_blocks = nrow(y),
    n_states = n_states,
    qbar_sim_blocks = qbar_sim_blocks,
    cov_inflation = cov_inflation,
    method = "direct_ecf_hmm"
  )
}

estimate_sigma1_rho_ecf_hmm_clock <- function(R, dt, alpha, sigma2,
                                              block_size,
                                              beta_grid = exp(seq(log(0.50),
                                                                  log(3.50),
                                                                  length.out = 15L)),
                                              rho_grid = exp(seq(log(0.40),
                                                                 log(3.50),
                                                                 length.out = 17L)),
                                              freq_mults = c(0.50, 0.75,
                                                             1.00),
                                              n_states = 31L,
                                              qbar_sim_blocks = 10000L,
                                              qbar_substeps_per_block = 10L,
                                              seed = 1L) {
  ecf_obs <- block_cf_log_observations(
    R = R,
    dt = dt,
    alpha = alpha,
    sigma2 = sigma2,
    block_size = block_size,
    freq_mults = freq_mults
  )
  est <- estimate_beta_rho_ecf_hmm(
    ecf_obs = ecf_obs,
    beta_grid = beta_grid,
    rho_grid = rho_grid,
    n_states = n_states,
    qbar_sim_blocks = qbar_sim_blocks,
    qbar_substeps_per_block = qbar_substeps_per_block,
    seed = seed
  )
  est$ecf_obs <- ecf_obs
  est$block_size <- ecf_obs$block_size
  est$method <- "direct_empirical_cf_observation_hmm_clock"
  est
}

estimate_rho_block_semigroup_gmm <- function(q_hat, block_h,
                                             tests = c("log", "inv1",
                                                       "inv4", "ratio1",
                                                       "exp025"),
                                             instruments = c(
                                               "a1", "center_log",
                                               "bounded_poly1",
                                               "bounded_poly2",
                                               "tanh_log",
                                               "logistic_low",
                                               "logistic_high", "low_bin",
                                               "mid_bin", "high_bin"
                                             ),
                                             q_thresholds = NULL,
                                             lags = c(1L, 2L),
                                             standardize = TRUE,
                                             rho_grid = exp(seq(log(0.40),
                                                                log(3.50),
                                                                length.out = 31L)),
                                             n_sim_blocks = NULL,
                                             n_sim_reps = 1L,
                                             substeps_per_block = 20L,
                                             burn_blocks = 1000L,
                                             noise_model = NULL,
                                             seed = 1L) {
  q_hat <- q_hat[is.finite(q_hat) & q_hat > 0]
  if (length(q_hat) < 20L) stop("Not enough activity values.")
  obs <- transition_moment_vector(
    q_hat = q_hat,
    tests = tests,
    instruments = instruments,
    q_thresholds = q_thresholds,
    lags = lags
  )
  if (is.null(n_sim_blocks)) {
    n_sim_blocks <- max(5000L, min(50000L, 10L * length(q_hat)))
  }
  n_sim_blocks <- as.integer(n_sim_blocks)
  n_sim_reps <- as.integer(n_sim_reps)
  if (n_sim_reps < 1L) n_sim_reps <- 1L
  obs_scale <- if (isTRUE(standardize)) obs$scale else rep(1, nrow(obs))
  obs_scale[!is.finite(obs_scale) | obs_scale <= 0] <- 1
  rows <- vector("list", length(rho_grid))
  model_store <- vector("list", length(rho_grid))
  for (i in seq_along(rho_grid)) {
    rho <- rho_grid[i]
    sim_reps <- vector("list", n_sim_reps)
    for (rep_id in seq_len(n_sim_reps)) {
      sim_q <- simulate_block_average_q(
        rho = rho,
        block_h = block_h,
        n_blocks = n_sim_blocks + max(as.integer(lags)),
        substeps_per_block = substeps_per_block,
        burn_blocks = burn_blocks,
        seed = seed + 1009L * i + 104729L * rep_id
      )
      if (!is.null(noise_model) &&
          !identical(noise_model$mode, "none")) {
        sim_q <- pmax(
          sim_q * exp(draw_activity_log_noise(noise_model, length(sim_q),
                                              q_ref = sim_q)),
          .Machine$double.eps
        )
      }
      sim_m <- transition_moment_vector(
        q_hat = sim_q,
        tests = tests,
        instruments = instruments,
        q_thresholds = q_thresholds,
        lags = lags
      )
      sim_reps[[rep_id]] <- sim_m[, c("key", "moment")]
    }
    sim_all <- Reduce(function(x, y) merge(x, y, by = "key", all = TRUE),
                      sim_reps)
    sim_vals <- as.matrix(sim_all[, -1L, drop = FALSE])
    sim_mean <- rowMeans(sim_vals, na.rm = TRUE)
    model <- data.frame(key = sim_all$key, model_moment = sim_mean)
    merged <- merge(obs, model, by = "key", all = FALSE)
    diff <- (merged$moment - merged$model_moment) / obs_scale[
      match(merged$key, obs$key)
    ]
    rows[[i]] <- data.frame(
      rho = rho,
      theta = rho ^ 2,
      objective = mean(diff ^ 2, na.rm = TRUE),
      n_moments = length(diff)
    )
    model_store[[i]] <- model
  }
  curve <- do.call(rbind, rows)
  best_i <- which.min(curve$objective)
  best <- curve[best_i, , drop = FALSE]
  best_model <- model_store[[best_i]]
  moments <- merge(obs, best_model, by = "key", all = FALSE)
  moments$diff_at_hat <- moments$moment - moments$model_moment
  moments$scale_used <- obs_scale[match(moments$key, obs$key)]
  list(
    rho_hat = best$rho,
    theta_hat = best$theta,
    objective = best$objective,
    curve = curve,
    moments = moments,
    tests = paste(tests, collapse = ","),
    instruments = paste(instruments, collapse = ","),
    lags = paste(lags, collapse = ","),
    standardize = isTRUE(standardize),
    n_sim_blocks = n_sim_blocks,
    n_sim_reps = n_sim_reps,
    substeps_per_block = substeps_per_block,
    noise_model = if (is.null(noise_model)) "none" else noise_model$mode,
    method = "block_semigroup_gmm"
  )
}

simulate_regime_summary <- function(rho, n_blocks, block_h, q_thresholds,
                                    lags, burn_blocks = 1000L,
                                    seed = 1L) {
  set.seed(seed)
  total <- n_blocks + burn_blocks
  v <- stationary_v_sample(1L)
  rho2 <- rho ^ 2
  q <- numeric(total)
  for (i in seq_len(total)) {
    d_v <- -rho2 * v * block_h +
      rho * sqrt(pmax(1 + v ^ 2, 1e-12)) * sqrt(block_h) *
      stats::rnorm(1L)
    v <- v + d_v
    q[i] <- 1 + v ^ 2
  }
  q <- q[(burn_blocks + 1L):total]
  activity_regime_stats(q, beta_hat = 1, q_thresholds, lags)$summary
}

estimate_rho_from_regime_clock <- function(obs_stats, rho_grid, n_blocks,
                                           block_h, q_thresholds, lags,
                                           n_sim = 3L, seed = 100L) {
  obs <- obs_stats$summary
  rows <- vector("list", length(rho_grid))
  for (i in seq_along(rho_grid)) {
    rho <- rho_grid[i]
    sim_summaries <- lapply(seq_len(n_sim), function(rep_id) {
      simulate_regime_summary(
        rho = rho,
        n_blocks = n_blocks,
        block_h = block_h,
        q_thresholds = q_thresholds,
        lags = lags,
        seed = seed + 1009L * i + 104729L * rep_id
      )
    })
    common <- Reduce(intersect, c(list(names(obs)), lapply(sim_summaries,
                                                           names)))
    common <- common[is.finite(obs[common])]
    sim_mat <- do.call(rbind, lapply(sim_summaries, function(x) x[common]))
    sim_mean <- colMeans(sim_mat, na.rm = TRUE)
    diff <- obs[common] - sim_mean
    rows[[i]] <- data.frame(
      rho = rho,
      objective = mean(diff ^ 2, na.rm = TRUE),
      n_features = length(common)
    )
  }
  curve <- do.call(rbind, rows)
  best <- curve[which.min(curve$objective), , drop = FALSE]
  list(rho_hat = best$rho, objective = best$objective, curve = curve)
}

estimate_sigma1_rho_clock <- function(R, dt, alpha, sigma2, block_size,
                                      rho_grid = exp(seq(log(0.40),
                                                         log(3.50),
                                                         length.out = 31L)),
                                      freq_mults = c(0.50, 0.75, 1.00),
                                      q_probs = c(0.25, 0.50, 0.75),
                                      beta_method = "wls",
                                      rho_method = c("generator_gmm",
                                                     "block_semigroup_gmm",
                                                     "eiv_block_semigroup_gmm",
                                                     "transition_clock"),
                                      generator_tests = c("log", "inv1",
                                                          "inv4", "ratio1",
                                                          "exp025"),
                                      generator_instruments = c(
                                        "a1", "center_log",
                                        "bounded_poly1", "bounded_poly2",
                                        "tanh_log", "logistic_low",
                                        "logistic_high", "low_bin",
                                        "mid_bin", "high_bin"
                                      ),
                                      generator_standardize = TRUE,
                                      regime_probs = c(0.40, 0.75),
                                      lags = c(1L, 2L),
                                      q_sample_n = 200000L,
                                      n_sim = 3L,
                                      semigroup_sim_blocks = NULL,
                                      semigroup_sim_reps = 1L,
                                      semigroup_substeps_per_block = 20L,
                                      seed = 1L) {
  rho_method <- match.arg(rho_method)
  q_sample <- stationary_q_sample(q_sample_n, seed = seed)
  block_tradeoff <- NULL
  if (length(block_size) > 1L) {
    block_choice <- select_clock_block_size(
      R = R,
      dt = dt,
      alpha = alpha,
      sigma2 = sigma2,
      candidate_block_sizes = block_size,
      freq_mults = freq_mults
    )
    block_tradeoff <- block_choice$table
    block_size <- block_choice$block_size
  }
  activity <- block_cf_activity(
    R = R,
    dt = dt,
    alpha = alpha,
    sigma2 = sigma2,
    block_size = block_size,
    freq_mults = freq_mults
  )
  beta <- estimate_beta_from_activity(activity$activity, q_sample,
                                      probs = q_probs,
                                      method = beta_method)
  q_thresholds <- stats::quantile(q_sample, probs = regime_probs,
                                  names = FALSE, type = 8)
  obs_stats <- activity_regime_stats(activity$activity, beta$beta_hat,
                                     q_thresholds, lags = lags)
  if (identical(rho_method, "generator_gmm")) {
    rho <- estimate_rho_generator_gmm(
      q_hat = obs_stats$q_hat,
      block_h = activity$block_h,
      tests = generator_tests,
      instruments = generator_instruments,
      q_thresholds = q_thresholds,
      lags = lags,
      standardize = generator_standardize,
      rho_grid = rho_grid
    )
  } else if (identical(rho_method, "block_semigroup_gmm")) {
    rho <- estimate_rho_block_semigroup_gmm(
      q_hat = obs_stats$q_hat,
      block_h = activity$block_h,
      tests = generator_tests,
      instruments = generator_instruments,
      q_thresholds = q_thresholds,
      lags = lags,
      standardize = generator_standardize,
      rho_grid = rho_grid,
      n_sim_blocks = semigroup_sim_blocks,
      n_sim_reps = semigroup_sim_reps,
      substeps_per_block = semigroup_substeps_per_block,
      seed = seed + 777L
    )
  } else if (identical(rho_method, "eiv_block_semigroup_gmm")) {
    noise_model <- estimate_activity_log_noise_model(
      activity_obj = activity,
      mode = "replicate_median"
    )
    rho <- estimate_rho_block_semigroup_gmm(
      q_hat = obs_stats$q_hat,
      block_h = activity$block_h,
      tests = generator_tests,
      instruments = generator_instruments,
      q_thresholds = q_thresholds,
      lags = lags,
      standardize = generator_standardize,
      rho_grid = rho_grid,
      n_sim_blocks = semigroup_sim_blocks,
      n_sim_reps = semigroup_sim_reps,
      substeps_per_block = semigroup_substeps_per_block,
      noise_model = noise_model,
      seed = seed + 777L
    )
  } else {
    rho <- estimate_rho_from_regime_clock(
      obs_stats = obs_stats,
      rho_grid = rho_grid,
      n_blocks = activity$n_blocks,
      block_h = activity$block_h,
      q_thresholds = q_thresholds,
      lags = lags,
      n_sim = n_sim,
      seed = seed + 999L
    )
  }
  list(
    sigma1_hat = beta$beta_hat / rho$rho_hat,
    rho_hat = rho$rho_hat,
    beta_hat = beta$beta_hat,
    beta_details = beta,
    rho_curve = rho$curve,
    rho_moments = if (is.null(rho$moments)) NULL else rho$moments,
    rho_method = rho_method,
    activity = activity$activity,
    q_thresholds = q_thresholds,
    obs_summary = obs_stats$summary,
    block_tradeoff = block_tradeoff,
    block_h = activity$block_h,
    block_size = activity$block_size,
    n_blocks = activity$n_blocks,
    method = "generator_calibrated_occupation_clock"
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (any(args == "--demo")) {
    sim <- simulate_residual_euler(
      n_steps = 100000L,
      terminal = 10,
      alpha = 1.25,
      sigma2 = 1.0,
      sigma1 = 1.2,
      rho = 1.45,
      seed = 1301
    )
    est <- estimate_sigma1_rho_clock(
      R = sim$R,
      dt = sim$dt,
      alpha = 1.25,
      sigma2 = 1.0,
      block_size = 500L,
      rho_grid = exp(seq(log(0.40), log(3.50), length.out = 17L)),
      q_sample_n = 50000L,
      n_sim = 2L,
      seed = 99L
    )
    print(data.frame(
      method = est$method,
      beta_hat = est$beta_hat,
      sigma1_hat = est$sigma1_hat,
      rho_hat = est$rho_hat,
      block_size = est$block_size,
      n_blocks = est$n_blocks
    ))
    print(est$rho_curve)
  }
}
