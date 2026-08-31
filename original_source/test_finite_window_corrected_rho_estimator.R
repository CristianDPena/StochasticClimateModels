source("finite_window_corrected_rho_estimator.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

q <- 1 / log(exp(1) + 10000000)
h <- floor(sqrt(10000000)) * 50 / 10000000

profile <- rho_fw_build_profile(
  q = q,
  h = h,
  rho_grid = seq(0.6, 3.4, by = 0.2),
  pair_count = 3000L,
  fine_steps = 128L,
  resolution_divisors = c(2L, 1L),
  seed = 941000L,
  chunk_size = 500L
)
finest <- rho_fw_finest_profile(profile)
evaluation <- rho_fw_evaluation_profile(profile)

assert_true(
  all(is.finite(profile$transfer)) &&
    all(profile$transfer > 0),
  "Transfer profile must be positive and finite."
)
assert_true(
  all(diff(finest$corrected_ratio) > 0),
  "Corrected rho map must be strictly increasing."
)
assert_true(
  all(diff(evaluation$corrected_ratio) > 0) &&
    all(evaluation$evaluation_method == "richardson_order1"),
  "Richardson-evaluated rho map must be strictly increasing."
)

for (rho_value in c(0.9, 1.45, 2.2)) {
  ratio <- stats::approx(
    evaluation$rho,
    evaluation$corrected_ratio,
    xout = rho_value,
    method = "linear"
  )$y
  inverse <- rho_fw_invert_ratio(ratio, profile)
  assert_true(
    inverse$status == "success" &&
      abs(inverse$rho_hat - rho_value) < 1e-12,
    "Correction inversion failed a grid-interpolation fixed point."
  )
}

lower <- rho_fw_invert_ratio(min(finest$corrected_ratio) / 2, profile)
upper <- rho_fw_invert_ratio(max(evaluation$corrected_ratio) * 2, profile)
assert_true(
  lower$status == "correction_domain_failure" &&
    lower$boundary == "lower" &&
    upper$status == "correction_domain_failure" &&
    upper$boundary == "upper",
  "Domain failures must be reported without boundary clipping."
)

dummy <- list(
  i_log_hat = stats::approx(
    evaluation$rho,
    evaluation$corrected_ratio,
    xout = 1.45
  )$y * 10,
  denominator_hat = 10,
  beta = 1.74,
  window_span = h,
  rho_hat = sqrt(
    stats::approx(
      evaluation$rho,
      evaluation$corrected_ratio,
      xout = 1.45
    )$y
  ),
  sigma1_hat = 1.74 / sqrt(
    stats::approx(
      evaluation$rho,
      evaluation$corrected_ratio,
      xout = 1.45
    )$y
  ),
  status = "success"
)
corrected <- rho_fw_correct_level4(dummy, profile)
assert_true(
  corrected$status == "success" &&
    abs(corrected$rho_hat - 1.45) < 1e-12 &&
    abs(corrected$sigma1_hat - 1.2) < 1e-12,
  "Level-4 correction wrapper failed its deterministic fixed point."
)

hidden_state_names <- c("V", "v_path", "C", "c_left", "rho_true")
formal_names <- names(formals(rho_fw_correct_level4))
assert_true(
  !any(hidden_state_names %in% formal_names),
  "Correction wrapper must not accept hidden-state or truth inputs."
)

cat("All finite-window correction implementation tests passed.\n")
