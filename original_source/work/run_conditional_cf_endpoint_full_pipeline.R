# Full observed-data validation for the conditional-CF endpoint HMM.

source("work/diagnose_exact_cf_endpoint_hmm.R")

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1L])
}

quick <- any(args == "--quick")
resolution_experiment <- any(
  args == "--validation-resolution-experiment"
)
legacy_oracle_tail <- any(args == "--oracle-tail")
tail_input_mode <- match.arg(option_value(
  "tail-input", if (legacy_oracle_tail) "oracle" else "estimated"
), c("estimated", "true-alpha", "true-sigma2", "oracle"))
case_index <- as.integer(option_value("case", 1L))
design_file <- option_value("design-file", "")
output_tag <- gsub("[^A-Za-z0-9_-]", "", option_value("output-tag", ""))
n_steps <- as.integer(option_value("n-steps",
                                    if (quick) 100000L else 10000000L))
terminal <- as.numeric(option_value("terminal", if (quick) 10 else 50))
if (!quick && n_steps < 10000000L && !resolution_experiment) {
  stop("Full-pipeline validation requires n_steps >= 10,000,000.")
}

cases <- data.frame(
  case = seq_len(9L),
  seed = c(2301L, 2311L, 2321L, 2401L, 2411L, 2421L, 2431L, 2441L, 2451L),
  alpha = c(1.25, 1.25, 1.25, 1.10, 1.50, 1.35, 1.15, 1.55, 1.40),
  sigma2 = c(1.00, 1.00, 1.00, 0.80, 1.20, 0.70, 1.30, 0.90, 1.10),
  sigma1 = c(1.20, 1.20, 1.20, 0.80, 1.60, 1.00, 1.50, 0.90, 1.40),
  rho = c(1.45, 1.45, 1.45, 1.80, 1.00, 1.70, 1.10, 1.90, 1.30)
)
panel_name <- "baseline"
run_id <- sprintf("B%02d", case_index)
if (nzchar(design_file)) {
  design <- utils::read.csv(design_file, stringsAsFactors = FALSE)
  required <- c("panel", "run_id", "case", "seed", "alpha", "sigma2",
                "sigma1", "rho", "n_steps", "terminal")
  if (!all(required %in% names(design))) {
    stop("Design file is missing required columns: ",
         paste(setdiff(required, names(design)), collapse = ", "))
  }
  cases <- design[, c("case", "seed", "alpha", "sigma2", "sigma1", "rho")]
  matched <- design[design$case == case_index, , drop = FALSE]
  if (nrow(matched) == 1L) {
    panel_name <- matched$panel
    run_id <- matched$run_id
    n_steps <- as.integer(matched$n_steps)
    terminal <- as.numeric(matched$terminal)
  }
}
if (!is.finite(case_index) || !case_index %in% cases$case) {
  stop("--case does not identify exactly one configured case.")
}
cfg <- cases[cases$case == case_index, , drop = FALSE]
if (nrow(cfg) != 1L) stop("Configured case identifiers must be unique.")
if (!quick && n_steps < 10000000L && !resolution_experiment) {
  stop("Full-pipeline validation requires n_steps >= 10,000,000.")
}

default_rho_grid <- paste(format(seq(0.7, 2.3, by = 0.1), nsmall = 1),
                          collapse = ",")
rho_grid <- as.numeric(strsplit(option_value("rho-grid", default_rho_grid),
                                ",", fixed = TRUE)[[1L]])
rho_grid <- sort(unique(rho_grid[is.finite(rho_grid) & rho_grid > 0]))
block_horizon <- as.numeric(option_value("block-horizon", 1 / 20))
n_states <- as.integer(option_value("n-states", if (quick) 9L else 21L))
transition_sim_blocks <- as.integer(option_value(
  "transition-sim-blocks", if (quick) 9000L else 630000L
))
substeps <- as.integer(option_value("substeps", if (quick) 8L else 32L))
max_path_nodes <- as.integer(option_value("max-path-nodes",
                                           if (quick) 5L else 31L))
kernel_seed <- as.integer(option_value("kernel-seed", 5000L))
beta_profile_grid_size <- as.integer(option_value("beta-profile-grid-size",
                                                   25L))
covariance_model <- match.arg(option_value(
  "covariance-model", "conditional_diag"
), c("conditional_diag", "conditional_full", "conditional_complex_full",
     "empirical_profiled"))
chunked_likelihood <- any(args == "--chunked-likelihood")
groups_per_chunk <- as.integer(option_value("groups-per-chunk", 16L))
if (chunked_likelihood) {
  if (!identical(covariance_model, "conditional_diag")) {
    stop("The exact-equivalent chunked evaluator currently supports conditional_diag only.")
  }
  conditional_exact_cf_forward_loglik_full_matrix <-
    conditional_exact_cf_forward_loglik
  source("work/chunked_exact_endpoint_likelihood.R")
  conditional_exact_cf_forward_loglik <- function(
      y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
      full_covariance = FALSE, return_filter = FALSE) {
    if (isTRUE(full_covariance)) {
      return(conditional_exact_cf_forward_loglik_full_matrix(
        y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
        full_covariance = TRUE, return_filter = return_filter
      ))
    }
    conditional_exact_cf_forward_loglik_chunked_library(
      y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
      groups_per_chunk = groups_per_chunk,
      return_filter = return_filter
    )
  }
}
kernel_cache_path <- option_value(
  "kernel-cache", "outputs/cache/conditional_cf_endpoint_kernel_production.rds"
)
kernel_models <- NULL
kernel_cache_used <- FALSE
if (nzchar(kernel_cache_path) && file.exists(kernel_cache_path)) {
  cache <- readRDS(kernel_cache_path)
  cache_config <- cache$config
  compatible <- isTRUE(all.equal(cache_config$block_horizon, block_horizon,
                                 tolerance = 1e-10)) &&
    identical(as.integer(cache_config$n_states), n_states) &&
    identical(as.integer(cache_config$transition_sim_blocks),
              transition_sim_blocks) &&
    identical(as.integer(cache_config$substeps), substeps) &&
    identical(as.integer(cache_config$max_path_nodes), max_path_nodes) &&
    identical(as.integer(cache_config$kernel_seed), kernel_seed)
  match_index <- if (compatible) {
    vapply(rho_grid, function(rho) {
      hit <- which(abs(cache_config$rho_grid - rho) < 1e-10)
      if (length(hit)) hit[1L] else NA_integer_
    }, integer(1L))
  } else rep(NA_integer_, length(rho_grid))
  if (compatible && all(is.finite(match_index))) {
    kernel_models <- cache$models[match_index]
    kernel_cache_used <- TRUE
  }
}

source_md5 <- unname(tools::md5sum(c(
  "work/diagnose_exact_cf_endpoint_hmm.R",
  "work/diagnose_markov_additive_endpoint_hmm.R",
  "work/run_conditional_cf_endpoint_full_pipeline.R"
)))
chunked_source_md5 <- if (chunked_likelihood) unname(tools::md5sum(
  "work/chunked_exact_endpoint_likelihood.R"
)) else NA_character_

cat(sprintf(
  paste0("CONDITIONAL-CF FULL PIPELINE case=%d seed=%d n=%d T=%.3f ",
         "true=(%.3f, %.3f, %.3f, %.3f)\n"),
  cfg$case, cfg$seed, n_steps, terminal, cfg$alpha, cfg$sigma2,
  cfg$sigma1, cfg$rho
))

simulation_time <- system.time({
  sim <- simulate_residual_euler(
    n_steps = n_steps, terminal = terminal,
    alpha = cfg$alpha, sigma2 = cfg$sigma2,
    sigma1 = cfg$sigma1, rho = cfg$rho,
    seed = cfg$seed, keep_v = FALSE
  )
})[["elapsed"]]

tail_time <- system.time({
  estimated_tail <- if (identical(tail_input_mode, "oracle")) {
    list(alpha_hat = cfg$alpha, sigma2_hat = cfg$sigma2,
         k_hat = NA_integer_, hill_score = NA_real_)
  } else {
    estimate_alpha_sigma2_hill(sim$R, sim$dt)
  }
  tail <- estimated_tail
  if (tail_input_mode %in% c("true-alpha", "oracle")) {
    tail$alpha_hat <- cfg$alpha
  }
  if (tail_input_mode %in% c("true-sigma2", "oracle")) {
    tail$sigma2_hat <- cfg$sigma2
  }
})[["elapsed"]]
cat(sprintf("Tail (%s): alpha=%.6f sigma2=%.6f k=%s\n",
            tail_input_mode, tail$alpha_hat, tail$sigma2_hat,
            as.character(tail$k_hat)))

checkpoint_root <- option_value("checkpoint-root", "outputs/checkpoints")
checkpoint_enabled <- !any(args == "--no-checkpoint") &&
  nzchar(checkpoint_root)
checkpoint_tag <- if (nzchar(output_tag)) output_tag else "default"
checkpoint_code <- paste(substr(source_md5, 1L, 6L), collapse = "-")
checkpoint_dir <- if (checkpoint_enabled) file.path(
  checkpoint_root,
  sprintf(
    "c%02d_s%d_%s_t%s_c%s_n%03d_src%s",
    cfg$case, cfg$seed, checkpoint_tag,
    substr(tail_input_mode, 1L, 1L),
    substr(covariance_model, 1L, 3L), max_path_nodes, checkpoint_code
  )
) else NULL
kernel_info <- if (kernel_cache_used) file.info(kernel_cache_path) else NULL
checkpoint_signature <- paste(
  "version=1",
  paste0("source_md5=", paste(source_md5, collapse = ",")),
  paste0("chunked_source_md5=", chunked_source_md5),
  paste0("case=", cfg$case), paste0("seed=", cfg$seed),
  paste0("n_steps=", n_steps), paste0("terminal=", format(terminal, digits = 17)),
  paste0("alpha_hat=", format(tail$alpha_hat, digits = 17)),
  paste0("sigma2_hat=", format(tail$sigma2_hat, digits = 17)),
  paste0("rho_grid=", paste(format(rho_grid, digits = 17), collapse = ",")),
  paste0("block_horizon=", format(block_horizon, digits = 17)),
  paste0("freq_mults=", paste(c(0.25, 0.50, 0.75), collapse = ",")),
  paste0("n_states=", n_states),
  paste0("transition_sim_blocks=", transition_sim_blocks),
  paste0("substeps=", substeps), paste0("max_path_nodes=", max_path_nodes),
  paste0("kernel_seed=", kernel_seed),
  paste0("beta_profile_grid_size=", beta_profile_grid_size),
  paste0("likelihood_evaluator=",
         if (chunked_likelihood) "chunked_exact" else "full_matrix_exact"),
  paste0("groups_per_chunk=", if (chunked_likelihood) groups_per_chunk else NA),
  paste0("covariance_model=", covariance_model),
  paste0("kernel_cache=", normalizePath(kernel_cache_path, mustWork = FALSE)),
  paste0("kernel_cache_size=", if (kernel_cache_used) kernel_info$size else NA),
  paste0("kernel_cache_mtime=", if (kernel_cache_used) {
    format(as.numeric(kernel_info$mtime), digits = 17)
  } else NA),
  sep = "|"
)

endpoint_time <- system.time({
  endpoint <- estimate_exact_cf_endpoint_hmm(
    R = sim$R, dt = sim$dt,
    alpha_hat = tail$alpha_hat, sigma2_hat = tail$sigma2_hat,
    block_horizon = block_horizon, rho_grid = rho_grid,
    freq_mults = c(0.25, 0.50, 0.75),
    n_states = n_states,
    transition_sim_blocks = transition_sim_blocks,
    substeps_per_block = substeps,
    max_path_nodes_per_pair = max_path_nodes,
    kernel_seed = kernel_seed,
    beta_profile_grid_size = beta_profile_grid_size,
    kernel_models = kernel_models,
    covariance_model = covariance_model,
    checkpoint_dir = checkpoint_dir,
    checkpoint_signature = checkpoint_signature,
    return_observations = FALSE,
    verbose = TRUE
  )
})[["elapsed"]]

profile_diag <- endpoint$profile_diagnostics
filter_diag <- endpoint$filter_diagnostics
alpha_error <- abs(tail$alpha_hat - cfg$alpha)
sigma2_error <- abs(tail$sigma2_hat - cfg$sigma2)
sigma1_error <- abs(endpoint$sigma1_hat - cfg$sigma1)
rho_error <- abs(endpoint$rho_hat - cfg$rho)
result <- data.frame(
  panel = panel_name,
  run_id = run_id,
  case = cfg$case,
  seed = cfg$seed,
  n_steps = n_steps,
  terminal = terminal,
  alpha_true = cfg$alpha,
  sigma2_true = cfg$sigma2,
  sigma1_true = cfg$sigma1,
  rho_true = cfg$rho,
  beta_true = cfg$sigma1 * cfg$rho,
  alpha_hat = tail$alpha_hat,
  sigma2_hat = tail$sigma2_hat,
  sigma1_hat = endpoint$sigma1_hat,
  rho_hat = endpoint$rho_hat,
  beta_hat = endpoint$beta_hat,
  variance_scale = endpoint$kappa_hat,
  alpha_abs_err = alpha_error,
  sigma2_abs_err = sigma2_error,
  sigma1_abs_err = sigma1_error,
  rho_abs_err = rho_error,
  beta_abs_err = abs(endpoint$beta_hat - cfg$sigma1 * cfg$rho),
  alpha_under_010 = alpha_error < 0.10,
  sigma2_under_010 = sigma2_error < 0.10,
  sigma1_under_010 = sigma1_error < 0.10,
  rho_under_010 = rho_error < 0.10,
  all_four_under_010 = all(c(alpha_error, sigma2_error,
                             sigma1_error, rho_error) < 0.10),
  k_hat = tail$k_hat,
  hill_score = tail$hill_score,
  block_size = endpoint$block_size,
  block_horizon = endpoint$block_horizon,
  n_blocks = endpoint$n_blocks,
  freq_mult_1 = endpoint$freq_mults[1L],
  freq_mult_2 = endpoint$freq_mults[2L],
  freq_mult_3 = endpoint$freq_mults[3L],
  rho_grid_n = length(rho_grid),
  rho_grid_spacing = if (length(rho_grid) > 1L) median(diff(rho_grid))
    else NA_real_,
  rho_support_lower = profile_diag$rho_support_lower,
  rho_support_upper = profile_diag$rho_support_upper,
  rho_profile_width = profile_diag$rho_profile_width,
  sigma1_support_lower = profile_diag$sigma1_support_lower,
  sigma1_support_upper = profile_diag$sigma1_support_upper,
  rho_profile_gap = profile_diag$rho_profile_gap,
  rho_at_grid_boundary = profile_diag$rho_at_grid_boundary,
  profile_quadratic_refined = endpoint$profile_quadratic_refined,
  profile_quadratic_curvature = endpoint$profile_quadratic_curvature,
  optimizer_convergence = endpoint$optimizer_convergence,
  optimizer_iterations = endpoint$optimizer_iterations,
  mean_filter_max_prob = if (is.null(filter_diag)) NA_real_
    else filter_diag$mean_filter_max_prob,
  mean_filter_entropy = if (is.null(filter_diag)) NA_real_
    else filter_diag$mean_filter_entropy,
  mean_filter_boundary_prob = if (is.null(filter_diag)) NA_real_
    else filter_diag$mean_filter_boundary_prob,
  warning = if (length(endpoint$warnings)) paste(endpoint$warnings,
                                                  collapse = ";") else "",
  covariance_model = endpoint$covariance_model,
  endpoint_state = "q",
  transition_sim_blocks = transition_sim_blocks,
  max_path_nodes = max_path_nodes,
  substeps = substeps,
  kernel_seed = kernel_seed,
  beta_profile_grid_size = beta_profile_grid_size,
  likelihood_evaluator = if (chunked_likelihood) {
    "chunked_exact_equivalent"
  } else {
    "full_matrix_exact"
  },
  groups_per_chunk = if (chunked_likelihood) groups_per_chunk else NA_integer_,
  kernel_cache_used = kernel_cache_used,
  kernel_cache_path = if (kernel_cache_used) kernel_cache_path else "",
  checkpoint_rows_reused = endpoint$checkpoint_rows_reused,
  checkpoint_dir = endpoint$checkpoint_dir,
  exact_cf_source_md5 = source_md5[1L],
  markov_core_source_md5 = source_md5[2L],
  runner_source_md5 = source_md5[3L],
  chunked_source_md5 = chunked_source_md5,
  tail_input_mode = tail_input_mode,
  simulation_seconds = simulation_time,
  tail_seconds = tail_time,
  endpoint_seconds = endpoint_time,
  total_seconds = simulation_time + tail_time + endpoint_time,
  peak_r_mb = max(gc(reset = FALSE)[, 6L]),
  numerical_convergence_status = "not_audited",
  profile_boundary_status = if (isTRUE(profile_diag$rho_at_grid_boundary))
    "rho_boundary" else "interior",
  failure_status = if (length(endpoint$warnings)) "warning" else "none",
  stringsAsFactors = FALSE
)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
prefix <- if (quick) "quick" else "full"
tail_tag <- if (identical(tail_input_mode, "estimated")) "" else
  paste0("_tail", gsub("-", "", tail_input_mode, fixed = TRUE))
output_tag <- if (nzchar(output_tag)) paste0("_", output_tag) else ""
csv_path <- file.path("outputs", sprintf(
  "conditional_cf_endpoint_%s_pipeline_case%02d_seed%d%s%s_%s.csv",
  prefix, cfg$case, cfg$seed, tail_tag, output_tag, stamp
))
rds_path <- sub("[.]csv$", ".rds", csv_path)
utils::write.csv(result, csv_path, row.names = FALSE)
saveRDS(list(result = result, profile = endpoint$profile,
             profile_diagnostics = profile_diag), rds_path)

print(result[, c("case", "alpha_hat", "sigma2_hat", "sigma1_hat", "rho_hat",
                 "alpha_abs_err", "sigma2_abs_err", "sigma1_abs_err",
                 "rho_abs_err", "all_four_under_010", "variance_scale",
                 "rho_profile_width", "warning")], row.names = FALSE)
cat("Wrote:", csv_path, "\n")
cat("Profile:", rds_path, "\n")
