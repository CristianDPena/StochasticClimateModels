# Build one reusable exact-CF endpoint kernel for a numerical-resolution audit.

source("work/diagnose_exact_cf_endpoint_hmm.R")

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1L])
}

rho <- as.numeric(option_value("rho", NA_real_))
max_path_nodes <- as.integer(option_value("max-path-nodes", NA_integer_))
if (!is.finite(rho) || rho <= 0) stop("--rho must be positive.")
if (!is.finite(max_path_nodes) || max_path_nodes < 1L) {
  stop("--max-path-nodes must be a positive integer.")
}

config <- list(
  block_horizon = as.numeric(option_value("block-horizon", 0.05)),
  n_states = as.integer(option_value("n-states", 21L)),
  transition_sim_blocks = as.integer(option_value(
    "transition-sim-blocks", 630000L
  )),
  substeps = as.integer(option_value("substeps", 32L)),
  max_path_nodes = max_path_nodes,
  kernel_seed = as.integer(option_value("kernel-seed", 5000L))
)
state_info <- make_q_endpoint_breaks(config$n_states)
elapsed <- system.time({
  model <- simulate_q_endpoint_cf_model(
    rho = rho,
    block_h = config$block_horizon,
    state_info = state_info,
    n_sim_blocks = config$transition_sim_blocks,
    substeps_per_block = config$substeps,
    max_path_nodes_per_pair = config$max_path_nodes,
    state_model = "conditional",
    seed = config$kernel_seed
  )
})[["elapsed"]]

dir.create("outputs/cache/components", recursive = TRUE, showWarnings = FALSE)
path <- file.path("outputs/cache/components", sprintf(
  "exact_cf_kernel_nodes%03d_rho%03d_seed%d_mc%d_sub%d_states%d.rds",
  max_path_nodes, as.integer(round(100 * rho)), config$kernel_seed,
  config$transition_sim_blocks, config$substeps, config$n_states
))
saveRDS(list(config = config, rho = rho, model = model,
             elapsed_seconds = elapsed), path)
cat(sprintf("Wrote %s in %.1f seconds\n", path, elapsed))
