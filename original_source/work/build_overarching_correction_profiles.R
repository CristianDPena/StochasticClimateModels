source("finite_window_corrected_rho_estimator.R")

manifest <- read.csv(
  "frozen_seed_manifest.csv",
  stringsAsFactors = FALSE
)
cells <- unique(manifest[, c(
  "n_steps", "terminal", "correction_profile"
)])
cells <- cells[order(cells$n_steps, cells$terminal), , drop = FALSE]
dir.create(
  "overarching_profiles",
  showWarnings = FALSE,
  recursive = TRUE
)

build_profile <- function(cell) {
  n_steps <- as.integer(cell$n_steps)
  terminal <- as.numeric(cell$terminal)
  path <- as.character(cell$correction_profile)
  q <- 1 / log(exp(1) + n_steps)
  h <- floor(sqrt(n_steps)) * terminal / n_steps
  started <- proc.time()[["elapsed"]]

  if (file.exists(path)) {
    profile <- read.csv(path, stringsAsFactors = FALSE)
    if (
      max(abs(profile$q - q)) > 1e-12 ||
        max(abs(profile$h - h)) > 1e-12
    ) {
      stop("Existing correction profile has incompatible q or H: ", path)
    }
    action <- "validated_existing"
  } else {
    primary <- rho_fw_build_profile(
      q = q,
      h = h,
      rho_grid = rho_fw_default_grid(),
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
      rho_grid = rho_fw_default_grid(),
      pair_count = 30000L,
      fine_steps = 256L,
      resolution_divisors = c(4L, 2L, 1L),
      seed = 941002L,
      chunk_size = 1000L
    )
    independent$replication <- "independent"
    profile <- rbind(primary, independent)
    write.csv(profile, path, row.names = FALSE)
    action <- "generated"
  }

  data.frame(
    n_steps = n_steps,
    terminal = terminal,
    q = q,
    h = h,
    correction_profile = path,
    action = action,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE
  )
}

workspace <- normalizePath(".", winslash = "/", mustWork = TRUE)
cluster <- parallel::makeCluster(3L)
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterExport(
  cluster,
  c("workspace"),
  envir = environment()
)
parallel::clusterEvalQ(cluster, {
  setwd(workspace)
  source("finite_window_corrected_rho_estimator.R")
  NULL
})
cell_jobs <- split(cells, seq_len(nrow(cells)))
build_rows <- parallel::parLapplyLB(
  cluster,
  cell_jobs,
  build_profile
)
parallel::stopCluster(cluster)
on.exit(NULL, add = FALSE)
build_summary <- do.call(rbind, build_rows)
write.csv(
  build_summary,
  "overarching_transfer_profile_build_log.csv",
  row.names = FALSE
)

validate_profile <- function(cell) {
  profile <- read.csv(
    cell$correction_profile,
    stringsAsFactors = FALSE
  )
  primary <- profile[
    profile$replication == "primary",
    ,
    drop = FALSE
  ]
  independent <- profile[
    profile$replication == "independent",
    ,
    drop = FALSE
  ]
  primary_evaluation <- rho_fw_evaluation_profile(primary)
  independent_evaluation <- rho_fw_evaluation_profile(independent)
  comparison <- merge(
    primary_evaluation[, c("rho", "transfer", "transfer_se")],
    independent_evaluation[, c("rho", "transfer", "transfer_se")],
    by = "rho",
    suffixes = c("_primary", "_independent")
  )
  z_difference <- (
    comparison$transfer_primary -
      comparison$transfer_independent
  ) / sqrt(
    comparison$transfer_se_primary^2 +
      comparison$transfer_se_independent^2
  )

  time_steps <- sort(unique(primary$time_steps))
  fine <- primary[
    primary$time_steps == max(time_steps),
    ,
    drop = FALSE
  ]
  coarse <- primary[
    primary$time_steps == sort(time_steps, decreasing = TRUE)[2L],
    ,
    drop = FALSE
  ]
  fine <- fine[order(fine$rho), , drop = FALSE]
  coarse <- coarse[order(coarse$rho), , drop = FALSE]
  fixed_rho <- c(0.90, 1.45, 2.20)
  fixed_error <- vapply(fixed_rho, function(rho) {
    ratio <- approx(
      primary_evaluation$rho,
      primary_evaluation$corrected_ratio,
      xout = rho
    )$y
    inversion <- rho_fw_invert_ratio(ratio, primary)
    abs(inversion$rho_hat - rho)
  }, numeric(1L))

  map_monotone <- !is.unsorted(
    primary_evaluation$corrected_ratio,
    strictly = TRUE
  )
  maximum_z <- max(abs(z_difference))
  maximum_resolution_difference <- max(abs(
    fine$transfer - coarse$transfer
  ))
  small_tau_distance <- abs(primary_evaluation$transfer[1L] - 1)
  status <- if (
    map_monotone &&
      max(fixed_error) < 1e-12 &&
      maximum_z <= 3 &&
      maximum_resolution_difference <= 0.02 &&
      small_tau_distance <= 0.05
  ) {
    "pass"
  } else {
    "fail"
  }

  data.frame(
    n_steps = cell$n_steps,
    terminal = cell$terminal,
    q = unique(primary$q),
    h = unique(primary$h),
    correction_profile = cell$correction_profile,
    map_monotone = map_monotone,
    maximum_fixed_point_error = max(fixed_error),
    maximum_independent_z_difference = maximum_z,
    maximum_128_to_256_difference =
      maximum_resolution_difference,
    small_tau_transfer_distance_from_one =
      small_tau_distance,
    status = status,
    stringsAsFactors = FALSE
  )
}

validation <- do.call(
  rbind,
  lapply(cell_jobs, validate_profile)
)
write.csv(
  validation,
  "overarching_transfer_profile_validation.csv",
  row.names = FALSE
)
print(build_summary)
print(validation)
if (any(validation$status != "pass")) {
  stop("At least one observation-design transfer profile failed validation.")
}
