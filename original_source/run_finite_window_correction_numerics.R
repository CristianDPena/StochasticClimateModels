source("finite_window_corrected_rho_estimator.R")

options(stringsAsFactors = FALSE)

n <- 10000000
terminal <- 50
k <- floor(sqrt(n))
h <- k * terminal / n
q <- 1 / log(exp(1) + n)
rho_grid <- rho_fw_default_grid()

reuse <- identical(Sys.getenv("RHO_FW_REUSE"), "1") &&
  file.exists("finite_window_correction_numerical_profiles.csv")
if (reuse) {
  profiles <- read.csv(
    "finite_window_correction_numerical_profiles.csv",
    stringsAsFactors = FALSE
  )
  primary <- profiles[
    profiles$replication == "primary",
    ,
    drop = FALSE
  ]
  independent <- profiles[
    profiles$replication == "independent",
    ,
    drop = FALSE
  ]
} else {
  primary <- rho_fw_build_profile(
    q = q,
    h = h,
    rho_grid = rho_grid,
    pair_count = 50000L,
    fine_steps = 256L,
    resolution_divisors = c(4L, 2L, 1L),
    seed = 941001L,
    chunk_size = 1000L
  )
  primary$replication <- "primary"

  independent <- rho_fw_build_profile(
    q = q,
    h = h,
    rho_grid = rho_grid,
    pair_count = 30000L,
    fine_steps = 256L,
    resolution_divisors = c(4L, 2L, 1L),
    seed = 941002L,
    chunk_size = 1000L
  )
  independent$replication <- "independent"

  profiles <- rbind(primary, independent)
  write.csv(
    profiles,
    "finite_window_correction_numerical_profiles.csv",
    row.names = FALSE
  )
}

primary_fine <- rho_fw_evaluation_profile(primary)
independent_fine <- rho_fw_evaluation_profile(independent)

comparison <- merge(
  primary_fine[, c(
    "rho", "tau", "transfer", "transfer_se", "corrected_ratio"
  )],
  independent_fine[, c(
    "rho", "transfer", "transfer_se", "corrected_ratio"
  )],
  by = "rho",
  suffixes = c("_primary", "_independent")
)
comparison$transfer_difference <-
  comparison$transfer_primary - comparison$transfer_independent
comparison$transfer_combined_se <- sqrt(
  comparison$transfer_se_primary^2 +
    comparison$transfer_se_independent^2
)
comparison$transfer_z_difference <- comparison$transfer_difference /
  comparison$transfer_combined_se

resolution_rows <- lapply(
  split(primary, primary$rho),
  function(group) {
    group <- group[order(group$time_steps), ]
    fine <- group[nrow(group), ]
    previous <- group[nrow(group) - 1L, ]
    data.frame(
      rho = fine$rho,
      tau = fine$tau,
      coarse_steps = previous$time_steps,
      fine_steps = fine$time_steps,
      transfer_coarse = previous$transfer,
      transfer_fine = fine$transfer,
      successive_difference =
        fine$transfer - previous$transfer,
      successive_difference_abs =
        abs(fine$transfer - previous$transfer),
      fine_transfer_se = fine$transfer_se
    )
  }
)
resolution <- do.call(rbind, resolution_rows)
row.names(resolution) <- NULL

fixed_rho <- c(0.90, 1.45, 2.20)
fixed_tests <- lapply(fixed_rho, function(rho_value) {
  expected_ratio <- stats::approx(
    primary_fine$rho,
    primary_fine$corrected_ratio,
    xout = rho_value,
    method = "linear"
  )$y
  inversion <- rho_fw_invert_ratio(expected_ratio, primary)
  data.frame(
    test = "profile_fixed_point",
    rho_input = rho_value,
    expected_ratio = expected_ratio,
    rho_recovered = inversion$rho_hat,
    absolute_error = abs(inversion$rho_hat - rho_value),
    status = inversion$status
  )
})
fixed_tests <- do.call(rbind, fixed_tests)

fixed_tests <- rbind(
  fixed_tests,
  data.frame(
    test = c(
      "corrected_map_monotone",
      "small_tau_transfer_near_one",
      "independent_replication_agreement",
      "successive_resolution_agreement"
    ),
    rho_input = c(NA, min(primary_fine$rho), NA, NA),
    expected_ratio = NA,
    rho_recovered = NA,
    absolute_error = c(
      as.numeric(!all(diff(primary_fine$corrected_ratio) > 0)),
      abs(primary_fine$transfer[1L] - 1),
      max(abs(comparison$transfer_z_difference)),
      max(resolution$successive_difference_abs)
    ),
    status = c(
      if (all(diff(primary_fine$corrected_ratio) > 0)) "pass" else "fail",
      if (abs(primary_fine$transfer[1L] - 1) <= 0.05) {
        "pass"
      } else {
        "fail"
      },
      if (max(abs(comparison$transfer_z_difference)) <= 3) {
        "pass"
      } else {
        "fail"
      },
      if (max(resolution$successive_difference_abs) <= 0.02) {
        "pass"
      } else {
        "fail"
      }
    )
  )
)

write.csv(
  fixed_tests,
  "finite_window_correction_fixed_point_tests.csv",
  row.names = FALSE
)
write.csv(
  comparison,
  "finite_window_correction_independent_replication.csv",
  row.names = FALSE
)
write.csv(
  resolution,
  "finite_window_correction_successive_resolution.csv",
  row.names = FALSE
)

print(fixed_tests)
