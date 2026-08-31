# Complete observed-data research pipeline:
# alpha -> sigma2 -> beta -> corrected rho -> sigma1.
#
# The sourced endpoint pipeline supplies the unchanged Hill tail inputs and
# beta coordinate. Its internal rho/sigma1 split is retained only as a
# diagnostic and is replaced in the primary output below.

wrapper_args <- commandArgs(trailingOnly = TRUE)
wrapper_option_value <- function(name, default = NULL) {
  hit <- grep(
    paste0("^--", name, "="),
    wrapper_args,
    value = TRUE
  )
  if (!length(hit)) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[1L])
}
correction_profile_path <- wrapper_option_value(
  "correction-profile",
  "finite_window_correction_numerical_profiles.csv"
)
validation_output_path <- wrapper_option_value(
  "validation-output",
  ""
)

source("rho_level4_oracle_residual_estimator.R")
source("finite_window_corrected_rho_estimator.R")
source("work/run_conditional_cf_endpoint_full_pipeline.R")

profile_rows <- read.csv(
  correction_profile_path,
  stringsAsFactors = FALSE
)
correction_profile <- profile_rows[
  profile_rows$replication == "primary",
  ,
  drop = FALSE
]
expected_q <- 1 / log(exp(1) + n_steps)
expected_h <- floor(sqrt(n_steps)) * terminal / n_steps
if (
  !nrow(correction_profile) ||
    max(abs(correction_profile$q - expected_q)) > 1e-12 ||
    max(abs(correction_profile$h - expected_h)) > 1e-12
) {
  stop(
    "Correction profile does not match the run's frozen q and H."
  )
}

rho_started <- proc.time()[["elapsed"]]
level4 <- rho_level4_estimate(
  dR = sim$R,
  dt = sim$dt[1L],
  alpha = tail$alpha_hat,
  sigma2 = tail$sigma2_hat,
  beta = endpoint$beta_hat,
  k = floor(sqrt(n_steps)),
  u = 1 / (
    endpoint$beta_hat * sqrt(log(exp(1) + n_steps))
  ),
  chunk_size = 100000L,
  return_observed_cache = FALSE
)
corrected <- rho_fw_correct_level4(
  level4,
  correction_profile
)
rho_seconds <- proc.time()[["elapsed"]] - rho_started

corrected_result <- result
corrected_result$beta_stage_internal_sigma1_hat <-
  corrected_result$sigma1_hat
corrected_result$beta_stage_internal_rho_hat <-
  corrected_result$rho_hat
corrected_result$beta_stage_internal_sigma1_abs_err <-
  corrected_result$sigma1_abs_err
corrected_result$beta_stage_internal_rho_abs_err <-
  corrected_result$rho_abs_err

corrected_result$sigma1_hat <- corrected$sigma1_hat
corrected_result$rho_hat <- corrected$rho_hat
corrected_result$sigma1_abs_err <-
  abs(corrected$sigma1_hat - cfg$sigma1)
corrected_result$rho_abs_err <-
  abs(corrected$rho_hat - cfg$rho)
corrected_result$sigma1_under_010 <-
  corrected_result$sigma1_abs_err < 0.10
corrected_result$rho_under_010 <-
  corrected_result$rho_abs_err < 0.10
corrected_result$all_four_under_010 <- all(c(
  corrected_result$alpha_abs_err,
  corrected_result$sigma2_abs_err,
  corrected_result$sigma1_abs_err,
  corrected_result$rho_abs_err
) < 0.10)

corrected_result$estimator_architecture <-
  "hill_tail__endpoint_beta__finite_window_rho__ratio_sigma1"
corrected_result$rho_estimator <-
  "level4_log_ecf_vov_finite_window_corrected"
corrected_result$rho_estimator_status <- level4$status
corrected_result$rho_correction_status <- corrected$status
corrected_result$rho_correction_boundary <- corrected$boundary
corrected_result$uncorrected_rho_hat <- level4$rho_hat
corrected_result$uncorrected_sigma1_hat <- level4$sigma1_hat
corrected_result$rho_transfer_hat <- corrected$transfer_hat
corrected_result$rho_observed_ratio <- corrected$observed_ratio
corrected_result$rho_numerator <- level4$i_log_hat
corrected_result$rho_denominator <- level4$denominator_hat
corrected_result$rho_corrected_expected_numerator <-
  corrected$corrected_expected_numerator
corrected_result$rho_frequency <- level4$u
corrected_result$rho_window_size <- level4$k
corrected_result$rho_window_span <- level4$window_span
corrected_result$rho_ecf_floor_fraction <-
  level4$ecf_floor_activation_fraction
corrected_result$rho_ecf_nonpositive_fraction <-
  level4$ecf_nonpositive_fraction
corrected_result$rho_seconds <- rho_seconds
corrected_result$total_seconds <-
  corrected_result$total_seconds + rho_seconds
corrected_result$peak_r_mb <- max(
  corrected_result$peak_r_mb,
  max(gc(reset = FALSE)[, 6L])
)

tail_diagnostic_started <- proc.time()[["elapsed"]]
tail_finite <- is.finite(sim$R) & sim$R != 0 &
  is.finite(sim$dt) & sim$dt > 0
tail_ordered <- sort(
  abs(sim$R[tail_finite]),
  decreasing = TRUE
)
sigma2_true_alpha <- estimate_sigma2_tail_robust(
  alpha_hat = cfg$alpha,
  A_sorted_desc = tail_ordered,
  k_hat = tail$k_hat,
  dt_vec = sim$dt[tail_finite]
)
corrected_result$tail_observation_count <- length(tail_ordered)
corrected_result$hill_k_min <- tail$hill_range[1L]
corrected_result$hill_k_max <- tail$hill_range[2L]
corrected_result$hill_small_sample <- tail$hill_small_sample
corrected_result$hill_k_at_candidate_boundary <-
  tail$k_hat %in% tail$hill_range
corrected_result$sigma2_k_lo <- tail$sigma2_details$k_lo
corrected_result$sigma2_k_hi <- tail$sigma2_details$k_hi
corrected_result$sigma2_tail_window_size <-
  tail$sigma2_details$window_size
corrected_result$sigma2_hat_true_alpha <-
  sigma2_true_alpha$sigma2_hat
corrected_result$sigma2_true_alpha_abs_err <-
  abs(sigma2_true_alpha$sigma2_hat - cfg$sigma2)
corrected_result$sigma2_alpha_propagation <-
  tail$sigma2_hat - sigma2_true_alpha$sigma2_hat
corrected_result$sigma2_residual_stage_error <-
  sigma2_true_alpha$sigma2_hat - cfg$sigma2
corrected_result$beta_first_order_sigma1_contribution <-
  (endpoint$beta_hat - cfg$sigma1 * cfg$rho) / cfg$rho
corrected_result$rho_first_order_sigma1_contribution <-
  -(cfg$sigma1 * cfg$rho) / cfg$rho^2 *
    (corrected$rho_hat - cfg$rho)
corrected_result$sigma1_first_order_remainder <-
  (corrected$sigma1_hat - cfg$sigma1) -
    corrected_result$beta_first_order_sigma1_contribution -
    corrected_result$rho_first_order_sigma1_contribution
corrected_result$sigma1_error_cancellation <- (
  corrected_result$sigma1_abs_err < 0.10 &&
    (
      corrected_result$beta_abs_err >= 0.10 ||
        corrected_result$rho_abs_err >= 0.10
    )
)
corrected_result$tail_diagnostic_seconds <-
  proc.time()[["elapsed"]] - tail_diagnostic_started
corrected_result$correction_profile_path <-
  normalizePath(
    correction_profile_path,
    winslash = "/",
    mustWork = TRUE
  )

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
stamp_corrected <- format(Sys.time(), "%Y%m%d_%H%M%S")
corrected_path <- file.path(
  "outputs",
  sprintf(
    paste0(
      "finite_window_corrected_full_pipeline_",
      "case%02d_seed%d_%s.csv"
    ),
    cfg$case,
    cfg$seed,
    stamp_corrected
  )
)
write.csv(
  corrected_result,
  corrected_path,
  row.names = FALSE
)
if (nzchar(validation_output_path)) {
  dir.create(
    dirname(validation_output_path),
    showWarnings = FALSE,
    recursive = TRUE
  )
  write.csv(
    corrected_result,
    validation_output_path,
    row.names = FALSE
  )
}

print(corrected_result[, c(
  "case", "alpha_hat", "sigma2_hat", "beta_hat",
  "rho_hat", "sigma1_hat",
  "alpha_abs_err", "sigma2_abs_err", "beta_abs_err",
  "rho_abs_err", "sigma1_abs_err", "all_four_under_010",
  "rho_estimator_status", "rho_correction_status"
)], row.names = FALSE)
cat("Corrected pipeline:", corrected_path, "\n")
