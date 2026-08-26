# Frozen Level-4 oracle-residual log-ECF estimator.
#
# The estimator functions in this file use only observed residual increments
# and the supplied alpha, sigma2, and beta. Oracle state arrays are accepted
# only by the separately named diagnostic functions.

rho_level4_design <- function() {
  regimes <- data.frame(
    regime = c(
      "case2",
      "low_clock_weak_jump",
      "high_clock_strong_jump"
    ),
    alpha = c(1.25, 1.45, 1.55),
    sigma2 = c(1.00, 0.60, 1.40),
    sigma1 = c(1.20, 1.60, 0.70),
    rho = c(1.45, 0.90, 2.20),
    stringsAsFactors = FALSE
  )
  regimes$beta <- regimes$sigma1 * regimes$rho
  list(
    regimes = regimes,
    terminal = 50,
    n_grid = c(250000L, 1000000L, 4000000L, 10000000L),
    stage_replicates = c(8L, 32L, 96L, 96L),
    maximum_stage_c_replicates = 192L,
    maximum_stage_d_replicates = 192L,
    default_chunk_size = 100000L
  )
}

rho_level4_seed <- function(regime_index, stage_index, replication) {
  as.integer(
    740000L + 10000L * regime_index +
      1000L * stage_index + replication
  )
}

rho_level4_stable_exponent <- function(
    u, dt, alpha, sigma2) {
  sigma2^alpha * abs(u)^alpha * dt^(1 - alpha / 2)
}

rho_level4_h_delta <- function(
    u, c_value, lambda_one, alpha) {
  lambda_two <- 2^alpha * lambda_one
  a <- u^2 * c_value
  2 / u^4 * (
    expm1(a + 2 * lambda_one) +
      expm1(-a + 2 * lambda_one - lambda_two)
  )
}

rho_level4_h_log <- function(
    u, c_value, lambda_one, alpha) {
  rho_level4_h_delta(u, c_value, lambda_one, alpha) / c_value^2
}

rho_level4_rolling_means <- function(values, k, central_index) {
  n <- length(values)
  k <- as.integer(k)
  central_index <- as.integer(central_index)
  if (k < 1L || 2L * k >= n) {
    stop("k must be positive and satisfy 2*k < length(values).")
  }
  if (
    any(central_index < k + 1L) ||
      any(central_index > n - k)
  ) {
    stop("central_index is outside the valid rolling range.")
  }
  cumulative <- c(0, cumsum(values))
  backward <- (
    cumulative[central_index] -
      cumulative[central_index - k]
  ) / k
  forward <- (
    cumulative[central_index + k + 1L] -
      cumulative[central_index + 1L]
  ) / k
  list(backward = backward, forward = forward)
}

rho_level4_local_from_cumulative <- function(
    cumulative_cos, central_index, k, deconvolution, u,
    ecf_floor = 1 / sqrt(k)) {
  backward_raw <- (
    cumulative_cos[central_index] -
      cumulative_cos[central_index - k]
  ) / k
  forward_raw <- (
    cumulative_cos[central_index + k + 1L] -
    cumulative_cos[central_index + 1L]
  ) / k
  backward_effective <- pmax(backward_raw, ecf_floor)
  forward_effective <- pmax(forward_raw, ecf_floor)
  backward_tilde <- deconvolution * backward_effective
  forward_tilde <- deconvolution * forward_effective
  backward_c <- rep(NA_real_, length(backward_tilde))
  forward_c <- rep(NA_real_, length(forward_tilde))
  backward_ok <- is.finite(backward_tilde) & backward_tilde > 0
  forward_ok <- is.finite(forward_tilde) & forward_tilde > 0
  backward_c[backward_ok] <-
    -2 * log(backward_tilde[backward_ok]) / u^2
  forward_c[forward_ok] <-
    -2 * log(forward_tilde[forward_ok]) / u^2
  list(
    backward_raw = backward_raw,
    forward_raw = forward_raw,
    backward_effective = backward_effective,
    forward_effective = forward_effective,
    backward_tilde = backward_tilde,
    forward_tilde = forward_tilde,
    backward_c = backward_c,
    forward_c = forward_c
  )
}

rho_level4_estimate <- function(
    dR, dt, alpha, sigma2, beta,
    k = floor(sqrt(length(dR))),
    u = 1 / (beta * sqrt(log(exp(1) + length(dR)))),
    chunk_size = 100000L,
    return_observed_cache = FALSE) {
  started <- proc.time()[["elapsed"]]
  n <- length(dR)
  k <- as.integer(k)
  chunk_size <- as.integer(chunk_size)
  if (n < 3L || any(!is.finite(dR))) {
    stop("dR must contain at least three finite increments.")
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    stop("dt must be one positive scalar.")
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 2) {
    stop("alpha must lie in (0, 2).")
  }
  if (!is.finite(sigma2) || sigma2 < 0) {
    stop("sigma2 must be nonnegative.")
  }
  if (!is.finite(beta) || beta <= 0) {
    stop("beta must be positive.")
  }
  if (!is.finite(u) || u <= 0) stop("u must be positive.")
  if (k < 1L || 2L * k >= n) {
    stop("k must be positive and satisfy 2*k < length(dR).")
  }
  if (chunk_size < 1L) stop("chunk_size must be positive.")

  lambda_one <- rho_level4_stable_exponent(
    u, dt, alpha, sigma2
  )
  deconvolution <- exp(lambda_one)
  standardized <- dR / sqrt(dt)
  cosine <- cos(u * standardized)
  cumulative_cos <- c(0, cumsum(cosine))
  rm(standardized, cosine)

  central_start <- k + 1L
  central_end <- n - k
  n_local <- central_end - central_start + 1L
  beta2 <- beta^2

  ecf_nonpositive_count <- 0
  ecf_floor_activation_count <- 0
  ecf_tilde_gt_one_count <- 0
  c_nonpositive_count <- 0
  c_below_beta2_count <- 0
  c_below_beta2_plus_count <- 0
  nonfinite_h_count <- 0
  numerator <- 0
  denominator_sum <- 0
  squared_log_difference_sum <- 0
  correction_sum <- 0
  raw_ecf_sum <- 0
  raw_ecf_square_sum <- 0
  ecf_floor <- 1 / sqrt(k)

  for (
    chunk_start in seq.int(central_start, central_end, by = chunk_size)
  ) {
    chunk_end <- min(central_end, chunk_start + chunk_size - 1L)
    central <- seq.int(chunk_start, chunk_end)
    local <- rho_level4_local_from_cumulative(
      cumulative_cos, central, k, deconvolution, u
    )
    all_raw <- c(local$backward_raw, local$forward_raw)
    all_tilde <- c(local$backward_tilde, local$forward_tilde)
    all_c <- c(local$backward_c, local$forward_c)

    ecf_nonpositive_count <- ecf_nonpositive_count +
      sum(!is.finite(all_raw) | all_raw <= 0)
    ecf_floor_activation_count <- ecf_floor_activation_count +
      sum(!is.finite(all_raw) | all_raw < ecf_floor)
    ecf_tilde_gt_one_count <- ecf_tilde_gt_one_count +
      sum(is.finite(all_tilde) & all_tilde > 1)
    c_nonpositive_count <- c_nonpositive_count +
      sum(!is.finite(all_c) | all_c <= 0)
    c_below_beta2_count <- c_below_beta2_count +
      sum(is.finite(all_c) & all_c > 0 & all_c < beta2)
    c_below_beta2_plus_count <- c_below_beta2_plus_count +
      sum(
        is.finite(local$forward_c) &
          local$forward_c > 0 &
          local$forward_c < beta2
      )
    raw_ecf_sum <- raw_ecf_sum + sum(all_raw[is.finite(all_raw)])
    raw_ecf_square_sum <- raw_ecf_square_sum +
      sum(all_raw[is.finite(all_raw)]^2)

    valid <- is.finite(local$backward_c) &
      local$backward_c > 0 &
      is.finite(local$forward_c) &
      local$forward_c > 0
    if (any(valid)) {
      c_backward <- local$backward_c[valid]
      c_forward <- local$forward_c[valid]
      log_difference <- log(c_forward) - log(c_backward)
      h_log <- rho_level4_h_log(
        u, c_forward, lambda_one, alpha
      )
      nonfinite_h_count <- nonfinite_h_count +
        sum(!is.finite(h_log) | h_log <= 0)
      finite_h <- is.finite(h_log) & h_log > 0
      if (any(finite_h)) {
        squared_component <-
          3 / (2 * k) * log_difference[finite_h]^2
        correction_component <- 3 / k^2 * h_log[finite_h]
        squared_log_difference_sum <- squared_log_difference_sum +
          sum(squared_component)
        correction_sum <- correction_sum + sum(correction_component)
        numerator <- numerator +
          sum(squared_component - correction_component)
      }
      denominator_sum <- denominator_sum +
        sum(1 - beta2 / c_forward)
    }
  }

  total_local_estimates <- 2 * n_local
  invalid_spot_count <- c_nonpositive_count
  all_local_valid <- invalid_spot_count == 0L &&
    nonfinite_h_count == 0L
  denominator <- if (all_local_valid) {
    4 * dt * denominator_sum
  } else {
    NA_real_
  }
  numerator_final <- if (all_local_valid) numerator else NA_real_

  primary_failure <- "none"
  status <- "success"
  if (c_nonpositive_count > 0L) {
    status <- primary_failure <- "spot_variance_failure"
  } else if (
    nonfinite_h_count > 0L ||
      !is.finite(numerator_final) ||
      !is.finite(denominator)
  ) {
    status <- primary_failure <- "numerical_failure"
  } else if (numerator_final <= 0) {
    status <- primary_failure <- "numerator_failure"
  } else if (denominator <= 0) {
    status <- primary_failure <- "denominator_failure"
  }

  rho_hat <- if (status == "success") {
    sqrt(numerator_final / denominator)
  } else {
    NA_real_
  }
  sigma1_hat <- if (
    is.finite(rho_hat) && rho_hat > 0
  ) {
    beta / rho_hat
  } else {
    NA_real_
  }
  if (
    status == "success" &&
      (!is.finite(rho_hat) || rho_hat <= 0 ||
        !is.finite(sigma1_hat) || sigma1_hat <= 0)
  ) {
    status <- primary_failure <- "numerical_failure"
    rho_hat <- sigma1_hat <- NA_real_
  }

  result <- list(
    estimator = "level4_log_ecf_vov",
    n = n,
    dt = dt,
    k = k,
    window_span = k * dt,
    u = u,
    alpha = alpha,
    sigma2 = sigma2,
    beta = beta,
    lambda_one = lambda_one,
    lambda_two = 2^alpha * lambda_one,
    n_local = n_local,
    i_log_hat = numerator_final,
    i_log_squared_component = squared_log_difference_sum,
    i_log_correction_component = correction_sum,
    denominator_hat = denominator,
    rho_hat = rho_hat,
    sigma1_hat = sigma1_hat,
    status = status,
    primary_failure = primary_failure,
    statistical_resolution_status = if (
      ecf_floor_activation_count > 0L
    ) {
      "ecf_floor_activated"
    } else {
      "resolved"
    },
    ecf_floor_threshold = ecf_floor,
    ecf_floor_activation_count = ecf_floor_activation_count,
    ecf_floor_activation_fraction =
      ecf_floor_activation_count / total_local_estimates,
    ecf_nonpositive_count = ecf_nonpositive_count,
    ecf_nonpositive_fraction =
      ecf_nonpositive_count / total_local_estimates,
    ecf_tilde_gt_one_count = ecf_tilde_gt_one_count,
    ecf_tilde_gt_one_fraction =
      ecf_tilde_gt_one_count / total_local_estimates,
    c_nonpositive_count = c_nonpositive_count,
    c_nonpositive_fraction =
      c_nonpositive_count / total_local_estimates,
    c_below_beta2_count = c_below_beta2_count,
    c_below_beta2_fraction =
      c_below_beta2_count / total_local_estimates,
    c_below_beta2_plus_count = c_below_beta2_plus_count,
    c_below_beta2_plus_fraction =
      c_below_beta2_plus_count / n_local,
    nonfinite_h_count = nonfinite_h_count,
    raw_ecf_mean = raw_ecf_sum / total_local_estimates,
    raw_ecf_variance = (
      raw_ecf_square_sum / total_local_estimates -
        (raw_ecf_sum / total_local_estimates)^2
    ),
    runtime_seconds = proc.time()[["elapsed"]] - started
  )
  if (isTRUE(return_observed_cache)) {
    result$observed_cache <- list(
      cumulative_cos = cumulative_cos,
      central_start = central_start,
      central_end = central_end,
      deconvolution = deconvolution
    )
  }
  result
}

rho_level4_oracle_diagnostics <- function(
    estimate, c_left, oracle_i_log_fine, oracle_d,
    rho_true, sigma1_true, chunk_size = 100000L) {
  if (is.null(estimate$observed_cache)) {
    stop("estimate must contain an observed cache.")
  }
  n <- estimate$n
  k <- estimate$k
  if (length(c_left) != n || any(!is.finite(c_left)) || any(c_left <= 0)) {
    stop("c_left must be a positive finite oracle diagnostic vector.")
  }
  cumulative_c <- c(0, cumsum(c_left))
  cache <- estimate$observed_cache
  log_error_sum <- 0
  log_error_square_sum <- 0
  log_hat_sum <- 0
  log_oracle_sum <- 0
  log_hat_square_sum <- 0
  log_oracle_square_sum <- 0
  log_cross_sum <- 0
  ecf_mean_residual_sum <- 0
  ecf_variance_residual_sum <- 0
  ecf_variance_ratio_sum <- 0
  valid_count <- 0L

  for (
    chunk_start in seq.int(
      cache$central_start, cache$central_end, by = chunk_size
    )
  ) {
    chunk_end <- min(
      cache$central_end, chunk_start + chunk_size - 1L
    )
    central <- seq.int(chunk_start, chunk_end)
    local <- rho_level4_local_from_cumulative(
      cache$cumulative_cos, central, k,
      cache$deconvolution, estimate$u
    )
    oracle_forward <- (
      cumulative_c[central + k + 1L] -
        cumulative_c[central + 1L]
    ) / k
    valid <- is.finite(local$forward_c) &
      local$forward_c > 0 &
      is.finite(oracle_forward) &
      oracle_forward > 0
    if (any(valid)) {
      log_hat <- log(local$forward_c[valid])
      log_oracle <- log(oracle_forward[valid])
      error <- log_hat - log_oracle
      count <- length(error)
      valid_count <- valid_count + count
      log_error_sum <- log_error_sum + sum(error)
      log_error_square_sum <- log_error_square_sum + sum(error^2)
      log_hat_sum <- log_hat_sum + sum(log_hat)
      log_oracle_sum <- log_oracle_sum + sum(log_oracle)
      log_hat_square_sum <- log_hat_square_sum + sum(log_hat^2)
      log_oracle_square_sum <-
        log_oracle_square_sum + sum(log_oracle^2)
      log_cross_sum <- log_cross_sum + sum(log_hat * log_oracle)

      predicted_raw <- exp(
        -0.5 * estimate$u^2 * oracle_forward[valid] -
          estimate$lambda_one
      )
      ecf_mean_residual_sum <- ecf_mean_residual_sum +
        sum(local$forward_raw[valid] - predicted_raw)
      h_oracle <- rho_level4_h_delta(
        estimate$u, oracle_forward[valid],
        estimate$lambda_one, estimate$alpha
      )
      observed_scaled_square <-
        k * (local$forward_c[valid] - oracle_forward[valid])^2
      finite_variance <- is.finite(h_oracle) & h_oracle > 0
      ecf_variance_residual_sum <- ecf_variance_residual_sum +
        sum(
          observed_scaled_square[finite_variance] -
            h_oracle[finite_variance]
        )
      ecf_variance_ratio_sum <- ecf_variance_ratio_sum +
        sum(
          observed_scaled_square[finite_variance] /
            h_oracle[finite_variance]
        )
    }
  }

  covariance_numerator <- log_cross_sum -
    log_hat_sum * log_oracle_sum / valid_count
  variance_hat <- log_hat_square_sum - log_hat_sum^2 / valid_count
  variance_oracle <-
    log_oracle_square_sum - log_oracle_sum^2 / valid_count
  correlation <- covariance_numerator /
    sqrt(variance_hat * variance_oracle)
  rho_error <- estimate$rho_hat - rho_true
  sigma1_error <- estimate$sigma1_hat - sigma1_true
  primary_failure <- estimate$primary_failure
  if (
    primary_failure == "none" &&
      (
        !is.finite(rho_error) ||
          !is.finite(sigma1_error) ||
          abs(rho_error) >= 0.10 ||
          abs(sigma1_error) >= 0.10
      )
  ) {
    primary_failure <- if (
      estimate$ecf_floor_activation_count > 0L
    ) {
      "ecf_resolution_failure"
    } else {
      "rho_statistical_error"
    }
  }

  list(
    oracle_i_log_fine = oracle_i_log_fine,
    oracle_i_log_model = rho_true^2 * oracle_d,
    oracle_d = oracle_d,
    i_log_error_fine = estimate$i_log_hat - oracle_i_log_fine,
    i_log_abs_error_fine =
      abs(estimate$i_log_hat - oracle_i_log_fine),
    i_log_error_model =
      estimate$i_log_hat - rho_true^2 * oracle_d,
    i_log_abs_error_model =
      abs(estimate$i_log_hat - rho_true^2 * oracle_d),
    denominator_error = estimate$denominator_hat - oracle_d,
    denominator_abs_error =
      abs(estimate$denominator_hat - oracle_d),
    rho_true = rho_true,
    rho_error = rho_error,
    rho_abs_error = abs(rho_error),
    sigma1_true = sigma1_true,
    sigma1_error = sigma1_error,
    sigma1_abs_error = abs(sigma1_error),
    sigma1_error_cancellation = FALSE,
    log_c_integrated_mean_error = log_error_sum / valid_count,
    log_c_integrated_rmse =
      sqrt(log_error_square_sum / valid_count),
    log_c_correlation = correlation,
    ecf_mean_residual = ecf_mean_residual_sum / valid_count,
    ecf_variance_residual =
      ecf_variance_residual_sum / valid_count,
    ecf_variance_ratio = ecf_variance_ratio_sum / valid_count,
    diagnostic_valid_fraction = valid_count / estimate$n_local,
    primary_failure = primary_failure
  )
}

rho_level4_drop_cache <- function(estimate) {
  estimate$observed_cache <- NULL
  estimate
}
