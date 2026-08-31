# Model-based finite-window correction for the frozen Level-4 estimator.
#
# This file does not alter the local ECF statistic. It evaluates the
# stationary semigroup transfer derived in finite_window_semigroup_derivation.md
# and inverts the corrected one-dimensional rho map.

rho_fw_default_grid <- function() {
  seq(0.40, 3.50, by = 0.10)
}

rho_fw_validate_grid <- function(rho_grid) {
  rho_grid <- as.numeric(rho_grid)
  if (
    length(rho_grid) < 3L ||
      any(!is.finite(rho_grid)) ||
      any(rho_grid <= 0) ||
      is.unsorted(rho_grid, strictly = TRUE)
  ) {
    stop("rho_grid must contain at least three strictly increasing positives.")
  }
  rho_grid
}

rho_fw_integrand <- function(x, q) {
  # cosh(x)^2 = 1 + sinh(x)^2, but this form stays accurate near zero.
  exp(-0.5 * q * cosh(x)^2)
}

rho_fw_aggregate_brownian <- function(z_fine, divisor, fine_dt) {
  divisor <- as.integer(divisor)
  fine_steps <- ncol(z_fine)
  if (
    divisor < 1L ||
      fine_steps %% divisor != 0L
  ) {
    stop("Each resolution divisor must divide fine_steps.")
  }
  coarse_steps <- fine_steps %/% divisor
  if (divisor == 1L) {
    return(z_fine * sqrt(fine_dt))
  }
  output <- matrix(0, nrow(z_fine), coarse_steps)
  for (index in seq_len(coarse_steps)) {
    columns <- ((index - 1L) * divisor + 1L):(index * divisor)
    output[, index] <- rowSums(
      z_fine[, columns, drop = FALSE]
    ) * sqrt(fine_dt)
  }
  output
}

rho_fw_additive_functionals <- function(
    x_initial, brownian_increments, dt, tau_grid, q) {
  path_count <- length(x_initial)
  step_count <- ncol(brownian_increments)
  if (nrow(brownian_increments) != path_count) {
    stop("Brownian increments and initial states have incompatible sizes.")
  }
  if (
    !is.finite(dt) || dt <= 0 ||
      !is.finite(q) || q <= 0 ||
      any(!is.finite(tau_grid)) ||
      any(tau_grid <= 0) ||
      max(tau_grid) > step_count * dt * (1 + 1e-12)
  ) {
    stop("Invalid additive-functional inputs.")
  }

  tau_grid <- as.numeric(tau_grid)
  output <- matrix(NA_real_, path_count, length(tau_grid))
  target_step <- pmax(1L, ceiling(tau_grid / dt - 1e-12))
  x <- as.numeric(x_initial)
  integral <- numeric(path_count)

  for (step in seq_len(step_count)) {
    left_integrand <- rho_fw_integrand(x, q)
    targets <- which(target_step == step)
    if (length(targets)) {
      interval_start <- (step - 1L) * dt
      for (target in targets) {
        remainder <- tau_grid[target] - interval_start
        remainder <- min(max(remainder, 0), dt)
        average <- (
          integral + remainder * left_integrand
        ) / tau_grid[target]
        output[, target] <- average
      }
    }
    integral <- integral + dt * left_integrand
    x <- x - 1.5 * tanh(x) * dt + brownian_increments[, step]
  }

  if (any(!is.finite(output)) || any(output <= 0) || any(output >= 1)) {
    stop("Nonfinite or out-of-range additive functional encountered.")
  }
  output
}

rho_fw_chunk_statistics <- function(
    pair_count, q, h, rho_grid, fine_steps,
    resolution_divisors, seed, chunk_size) {
  rho_grid <- rho_fw_validate_grid(rho_grid)
  pair_count <- as.integer(pair_count)
  fine_steps <- as.integer(fine_steps)
  resolution_divisors <- as.integer(resolution_divisors)
  chunk_size <- as.integer(chunk_size)
  if (
    pair_count < 2L ||
      fine_steps < 4L ||
      chunk_size < 1L ||
      any(resolution_divisors < 1L) ||
      any(fine_steps %% resolution_divisors != 0L)
  ) {
    stop("Invalid Monte Carlo resolution inputs.")
  }

  tau_grid <- rho_grid^2 * h
  resolution_names <- paste0(
    "steps_", fine_steps %/% resolution_divisors
  )
  grid_count <- length(rho_grid)
  resolution_count <- length(resolution_divisors)

  accumulators <- lapply(seq_len(resolution_count), function(unused) {
    list(
      count = 0L,
      sum_pair = numeric(grid_count),
      sum_den = numeric(grid_count),
      sum_pair2 = numeric(grid_count),
      sum_den2 = numeric(grid_count),
      sum_cross = numeric(grid_count)
    )
  })

  set.seed(seed)
  starts <- seq.int(1L, pair_count, by = chunk_size)
  for (chunk_start in starts) {
    current_pairs <- min(chunk_size, pair_count - chunk_start + 1L)
    z0 <- stats::rt(current_pairs, df = 3) / sqrt(3)
    x0 <- asinh(z0)
    x_initial <- c(x0, x0)
    z_fine <- matrix(
      stats::rnorm(2L * current_pairs * fine_steps),
      nrow = 2L * current_pairs,
      ncol = fine_steps
    )

    for (resolution in seq_len(resolution_count)) {
      divisor <- resolution_divisors[resolution]
      standardized_brownian <- rho_fw_aggregate_brownian(
        z_fine, divisor, 1
      )
      accumulator <- accumulators[[resolution]]
      accumulator$count <- accumulator$count + current_pairs
      first <- seq_len(current_pairs)
      second <- current_pairs + first

      for (grid_index in seq_len(grid_count)) {
        fine_dt <- tau_grid[grid_index] / fine_steps
        coarse_dt <- divisor * fine_dt
        brownian <- standardized_brownian * sqrt(fine_dt)
        additive <- rho_fw_additive_functionals(
          x_initial, brownian, coarse_dt,
          tau_grid[grid_index], q
        )
        q_functional <- -2 / q * log(additive[, 1L])
        l_functional <- log(q_functional)
        pair_stat <- 0.5 * (
          l_functional[first] - l_functional[second]
        )^2
        denominator_stat <- 0.5 * (
          1 - 1 / q_functional[first] +
            1 - 1 / q_functional[second]
        )

        accumulator$sum_pair[grid_index] <-
          accumulator$sum_pair[grid_index] + sum(pair_stat)
        accumulator$sum_den[grid_index] <-
          accumulator$sum_den[grid_index] + sum(denominator_stat)
        accumulator$sum_pair2[grid_index] <-
          accumulator$sum_pair2[grid_index] + sum(pair_stat^2)
        accumulator$sum_den2[grid_index] <-
          accumulator$sum_den2[grid_index] + sum(denominator_stat^2)
        accumulator$sum_cross[grid_index] <-
          accumulator$sum_cross[grid_index] +
            sum(pair_stat * denominator_stat)
      }
      accumulators[[resolution]] <- accumulator
      rm(standardized_brownian)
    }
    rm(z_fine)
  }

  rows <- vector("list", resolution_count)
  for (resolution in seq_len(resolution_count)) {
    accumulator <- accumulators[[resolution]]
    count <- accumulator$count
    mean_pair <- accumulator$sum_pair / count
    mean_den <- accumulator$sum_den / count
    var_pair <- pmax(
      (
        accumulator$sum_pair2 -
          count * mean_pair^2
      ) / (count - 1L),
      0
    )
    var_den <- pmax(
      (
        accumulator$sum_den2 -
          count * mean_den^2
      ) / (count - 1L),
      0
    )
    covariance <- (
      accumulator$sum_cross -
        count * mean_pair * mean_den
    ) / (count - 1L)

    numerator_rate_unit <- 3 / tau_grid * mean_pair
    denominator_rate_unit <- 4 * mean_den
    transfer <- numerator_rate_unit / denominator_rate_unit
    derivative_pair <- (3 / tau_grid) / denominator_rate_unit
    derivative_den <- -numerator_rate_unit * 4 /
      denominator_rate_unit^2
    influence_variance <- pmax(
      derivative_pair^2 * var_pair +
        derivative_den^2 * var_den +
        2 * derivative_pair * derivative_den * covariance,
      0
    )
    transfer_se <- sqrt(influence_variance / count)

    rows[[resolution]] <- data.frame(
      numerical_seed = as.integer(seed),
      pair_count = count,
      resolution = resolution_names[resolution],
      time_steps = fine_steps %/% resolution_divisors[resolution],
      time_step =
        tau_grid / (fine_steps %/% resolution_divisors[resolution]),
      q = q,
      h = h,
      rho = rho_grid,
      tau = tau_grid,
      conditional_variance_mean = mean_pair,
      conditional_variance_se = sqrt(var_pair / count),
      denominator_moment_mean = mean_den,
      denominator_moment_se = sqrt(var_den / count),
      numerator_rate_unit = numerator_rate_unit,
      denominator_rate_unit = denominator_rate_unit,
      transfer = transfer,
      transfer_se = transfer_se,
      corrected_ratio = rho_grid^2 * transfer,
      corrected_ratio_se = rho_grid^2 * transfer_se,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

rho_fw_build_profile <- function(
    q, h,
    rho_grid = rho_fw_default_grid(),
    pair_count = 50000L,
    fine_steps = 1024L,
    resolution_divisors = c(4L, 2L, 1L),
    seed = 941001L,
    chunk_size = 2000L) {
  started <- proc.time()[["elapsed"]]
  if (
    !is.finite(q) || q <= 0 ||
      !is.finite(h) || h <= 0
  ) {
    stop("q and h must be positive finite scalars.")
  }
  result <- rho_fw_chunk_statistics(
    pair_count = pair_count,
    q = q,
    h = h,
    rho_grid = rho_grid,
    fine_steps = fine_steps,
    resolution_divisors = resolution_divisors,
    seed = seed,
    chunk_size = chunk_size
  )
  result$runtime_seconds <- proc.time()[["elapsed"]] - started
  result
}

rho_fw_finest_profile <- function(profile) {
  required <- c(
    "rho", "q", "h", "time_steps", "transfer",
    "transfer_se", "corrected_ratio", "corrected_ratio_se"
  )
  if (!all(required %in% names(profile))) {
    stop("Correction profile is missing required columns.")
  }
  maximum_steps <- max(profile$time_steps)
  finest <- profile[
    profile$time_steps == maximum_steps,
    ,
    drop = FALSE
  ]
  finest <- finest[order(finest$rho), , drop = FALSE]
  if (
    any(!is.finite(finest$corrected_ratio)) ||
      is.unsorted(finest$corrected_ratio, strictly = TRUE)
  ) {
    stop("The corrected rho map is not strictly increasing.")
  }
  finest
}

rho_fw_evaluation_profile <- function(profile) {
  finest <- rho_fw_finest_profile(profile)
  step_grid <- sort(unique(profile$time_steps))
  if (length(step_grid) < 2L) {
    finest$evaluation_method <- "finest_resolution"
    finest$discretization_increment <- NA_real_
    return(finest)
  }

  fine_steps <- tail(step_grid, 1L)
  coarse_steps <- tail(step_grid, 2L)[1L]
  ratio <- fine_steps / coarse_steps
  if (abs(ratio - 2) > 1e-12) {
    stop("Richardson evaluation requires a final 2:1 resolution ratio.")
  }
  coarse <- profile[
    profile$time_steps == coarse_steps,
    ,
    drop = FALSE
  ]
  coarse <- coarse[order(coarse$rho), , drop = FALSE]
  if (
    nrow(coarse) != nrow(finest) ||
      max(abs(coarse$rho - finest$rho)) > 1e-12
  ) {
    stop("Fine and coarse correction grids do not match.")
  }

  evaluation <- finest
  increment <- finest$transfer - coarse$transfer
  evaluation$transfer <- finest$transfer + increment
  evaluation$corrected_ratio <-
    evaluation$rho^2 * evaluation$transfer
  evaluation$discretization_increment <- increment
  evaluation$evaluation_method <- "richardson_order1"
  if (
    any(!is.finite(evaluation$transfer)) ||
      any(evaluation$transfer <= 0) ||
      is.unsorted(evaluation$corrected_ratio, strictly = TRUE)
  ) {
    stop("Richardson-corrected rho map is invalid or nonmonotone.")
  }
  evaluation
}

rho_fw_invert_ratio <- function(observed_ratio, profile) {
  finest <- rho_fw_evaluation_profile(profile)
  if (!is.finite(observed_ratio) || observed_ratio <= 0) {
    return(list(
      rho_hat = NA_real_,
      status = "invalid_observed_ratio",
      boundary = NA_character_
    ))
  }
  map_range <- range(finest$corrected_ratio)
  if (observed_ratio < map_range[1L]) {
    return(list(
      rho_hat = NA_real_,
      status = "correction_domain_failure",
      boundary = "lower"
    ))
  }
  if (observed_ratio > map_range[2L]) {
    return(list(
      rho_hat = NA_real_,
      status = "correction_domain_failure",
      boundary = "upper"
    ))
  }
  rho_hat <- stats::approx(
    x = finest$corrected_ratio,
    y = finest$rho,
    xout = observed_ratio,
    method = "linear",
    ties = "ordered"
  )$y
  transfer_hat <- stats::approx(
    x = finest$rho,
    y = finest$transfer,
    xout = rho_hat,
    method = "linear",
    ties = "ordered"
  )$y
  list(
    rho_hat = rho_hat,
    transfer_hat = transfer_hat,
    corrected_ratio_hat = rho_hat^2 * transfer_hat,
    status = "success",
    boundary = "none"
  )
}

rho_fw_correct_level4 <- function(uncorrected_estimate, profile) {
  required <- c(
    "i_log_hat", "denominator_hat", "beta", "window_span",
    "rho_hat", "sigma1_hat", "status"
  )
  if (!all(required %in% names(uncorrected_estimate))) {
    stop("uncorrected_estimate is missing required Level-4 fields.")
  }
  observed_ratio <- (
    uncorrected_estimate$i_log_hat /
      uncorrected_estimate$denominator_hat
  )
  inversion <- rho_fw_invert_ratio(observed_ratio, profile)
  rho_hat <- inversion$rho_hat
  sigma1_hat <- if (is.finite(rho_hat) && rho_hat > 0) {
    uncorrected_estimate$beta / rho_hat
  } else {
    NA_real_
  }
  list(
    estimator = "level4_log_ecf_vov_finite_window_corrected",
    observed_ratio = observed_ratio,
    uncorrected_rho_hat = uncorrected_estimate$rho_hat,
    uncorrected_sigma1_hat = uncorrected_estimate$sigma1_hat,
    rho_hat = rho_hat,
    sigma1_hat = sigma1_hat,
    transfer_hat = if (!is.null(inversion$transfer_hat)) {
      inversion$transfer_hat
    } else {
      NA_real_
    },
    corrected_expected_numerator = if (
      is.finite(rho_hat) && is.finite(inversion$transfer_hat)
    ) {
      uncorrected_estimate$denominator_hat *
        rho_hat^2 * inversion$transfer_hat
    } else {
      NA_real_
    },
    uncorrected_numerator = uncorrected_estimate$i_log_hat,
    denominator = uncorrected_estimate$denominator_hat,
    beta = uncorrected_estimate$beta,
    h = uncorrected_estimate$window_span,
    q = unique(profile$q),
    status = if (uncorrected_estimate$status != "success") {
      paste0("uncorrected_", uncorrected_estimate$status)
    } else {
      inversion$status
    },
    boundary = inversion$boundary
  )
}
