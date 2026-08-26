# Assemble independently generated rho components into one runner-compatible cache.

args <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1L])
}

max_path_nodes <- as.integer(option_value("max-path-nodes", NA_integer_))
if (!is.finite(max_path_nodes)) stop("--max-path-nodes is required.")
kernel_seed <- as.integer(option_value("kernel-seed", 5000L))
transition_sim_blocks <- as.integer(option_value(
  "transition-sim-blocks", 630000L
))
substeps <- as.integer(option_value("substeps", 32L))
n_states <- as.integer(option_value("n-states", 21L))
rho_grid <- as.numeric(strsplit(
  option_value(
    "rho-grid",
    paste(format(seq(0.7, 2.3, by = 0.1), nsmall = 1), collapse = ",")
  ), ",", fixed = TRUE
)[[1L]])
rho_grid <- sort(unique(rho_grid[is.finite(rho_grid) & rho_grid > 0]))
if (!length(rho_grid)) stop("--rho-grid must contain positive values.")
paths <- file.path("outputs/cache/components", sprintf(
  "exact_cf_kernel_nodes%03d_rho%03d_seed%d_mc%d_sub%d_states%d.rds",
  max_path_nodes, as.integer(round(100 * rho_grid)), kernel_seed,
  transition_sim_blocks, substeps, n_states
))
missing <- paths[!file.exists(paths)]
if (length(missing)) stop("Missing components: ", paste(missing, collapse = ", "))
parts <- lapply(paths, readRDS)
reference <- parts[[1L]]$config
if (!all(vapply(parts, function(x) identical(x$config, reference), logical(1L)))) {
  stop("Component configurations differ.")
}
if (any(abs(vapply(parts, `[[`, numeric(1L), "rho") - rho_grid) > 1e-10)) {
  stop("Component rho values do not match the frozen grid.")
}
config <- c(list(rho_grid = rho_grid), reference)
cache <- list(config = config, models = lapply(parts, `[[`, "model"),
              component_seconds = vapply(parts, `[[`, numeric(1L),
                                           "elapsed_seconds"))
default_path <- file.path("outputs/cache", sprintf(
  "conditional_cf_endpoint_kernel_nodes%03d_seed%d_mc%d_sub%d_states%d.rds",
  max_path_nodes, kernel_seed, transition_sim_blocks, substeps, n_states
))
path <- option_value("output", default_path)
saveRDS(cache, path)
cat("Wrote", path, "\n")
