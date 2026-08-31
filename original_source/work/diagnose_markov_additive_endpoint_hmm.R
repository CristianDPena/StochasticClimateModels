# Markov-additive endpoint-state empirical-CF HMM diagnostics.
#
# This prototype corrects the failed Qbar-state HMM by using endpoint V bins
# as the Markov state and treating Qbar as a block additive functional.

source("work/occupation_transition_clock_prototype.R")

args <- commandArgs(trailingOnly = TRUE)
quick <- any(args == "--quick")
focused <- any(args == "--focused")
recommended <- any(args == "--recommended")
corrected <- any(args == "--corrected") || recommended
option_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1L])
}
model_label <- option_value(
  "label", if (recommended) "recommended_endpoint_hmm" else "conditional_joint"
)
state_model <- match.arg(option_value("state-model", "conditional"),
                         c("conditional", "collocation"))
endpoint_state <- match.arg(option_value(
  "endpoint-state", if (recommended) "q" else "v"
), c("v", "q"))
kernel_method <- match.arg(option_value(
  "kernel-method", "monte_carlo"
), c("monte_carlo", "feynman_kac", "hybrid"))
emission_dispersion <- as.numeric(option_value("emission-dispersion", 1.0))
if (!is.finite(emission_dispersion) || emission_dispersion <= 0) {
  stop("--emission-dispersion must be a positive finite number.")
}
continuous_dispersion <- any(args == "--continuous-dispersion") || recommended
continuous_beta <- any(args == "--continuous-beta") || recommended
profile_engine <- match.arg(option_value(
  "profile-engine", if (recommended) "em" else "optim"
), c("em", "optim"))
discrete_dispersion <- any(args == "--profile-dispersion") &&
  !continuous_dispersion
dispersion_grid <- if (discrete_dispersion) {
  c(0.75, 1.00, 1.50, 2.00, 3.00, 4.00, 6.00)
} else {
  emission_dispersion
}
dispersion_bounds <- as.numeric(strsplit(option_value(
  "dispersion-bounds", "0.25,Inf"
), ",", fixed = TRUE)[[1L]])
if (length(dispersion_bounds) != 2L ||
    any(is.na(dispersion_bounds)) ||
    dispersion_bounds[1L] <= 0 ||
    dispersion_bounds[2L] <= dispersion_bounds[1L] ||
    dispersion_bounds[2L] == -Inf) {
  stop(paste0("--dispersion-bounds must contain a positive lower bound and ",
              "a larger upper bound (Inf is allowed)."))
}
dispersion_tol <- as.numeric(option_value("dispersion-tol", 0.02))
if (!is.finite(dispersion_tol) || dispersion_tol <= 0) {
  stop("--dispersion-tol must be positive.")
}
beta_bounds <- as.numeric(strsplit(option_value(
  "beta-bounds", "0.20,5.00"
), ",", fixed = TRUE)[[1L]])
if (length(beta_bounds) != 2L || any(!is.finite(beta_bounds)) ||
    beta_bounds[1L] <= 0 || beta_bounds[2L] <= beta_bounds[1L]) {
  stop("--beta-bounds must be two positive increasing numbers.")
}
covariance_structure <- match.arg(option_value(
  "covariance", "diagonal"
), c("diagonal", "full"))
tail_input_mode <- match.arg(option_value(
  "tail-input", if (any(args == "--estimate-tail")) "estimated" else "oracle"
), c("oracle", "estimated", "true-alpha", "true-sigma2"))
observed_only <- any(args == "--observed-only")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

fixed <- list(
  alpha = as.numeric(option_value("alpha", 1.25)),
  sigma2 = as.numeric(option_value("sigma2", 1.0)),
  sigma1 = as.numeric(option_value("sigma1", 1.2)),
  rho = as.numeric(option_value("rho", 1.45))
)
beta_true <- fixed$sigma1 * fixed$rho

cases <- if (quick) {
  expand.grid(seed = c(1301L), n_steps = c(100000L), terminal = c(10),
              block_size = c(500L), KEEP.OUT.ATTRS = FALSE)
} else if (focused) {
  expand.grid(seed = c(1301L, 1311L), n_steps = c(300000L),
              terminal = c(20), block_size = c(500L),
              KEEP.OUT.ATTRS = FALSE)
} else {
  expand.grid(seed = c(1301L, 1311L), n_steps = c(300000L, 600000L),
              terminal = c(20), block_size = c(500L),
              KEEP.OUT.ATTRS = FALSE)
}

freq_mults <- if (recommended || corrected) {
  c(0.25, 0.50, 0.75)
} else {
  c(0.50, 0.75, 1.00)
}
sigma1_grid <- c(0.60, 0.80, 1.00, 1.20, 1.40, 1.70, 2.00, 2.40)
rho_grid <- c(0.80, 1.00, 1.20, 1.35, 1.45, 1.60, 1.80, 2.10, 2.40)
if (any(args == "--coarse-grid")) {
  sigma1_grid <- c(1.00, 1.20, 1.40, 1.70)
  rho_grid <- c(1.20, 1.35, 1.45, 1.60, 1.80)
}
if (any(args == "--local-grid")) {
  sigma1_grid <- c(1.00, 1.20, 1.40)
  rho_grid <- c(1.35, 1.45, 1.60)
}
if (any(args == "--diagnostic-grid")) {
  sigma1_grid <- c(0.80, 1.00, 1.20, 1.40, 1.70)
  rho_grid <- c(1.20, 1.35, 1.45, 1.60, 1.80)
}
rho_grid_option <- option_value("rho-grid", NULL)
if (!is.null(rho_grid_option)) {
  rho_grid <- sort(unique(as.numeric(strsplit(rho_grid_option, ",",
                                               fixed = TRUE)[[1L]])))
  rho_grid <- rho_grid[is.finite(rho_grid) & rho_grid > 0]
  if (length(rho_grid) < 3L) {
    stop("--rho-grid must contain at least three positive values.")
  }
}
n_states <- if (quick) 11L else 13L
transition_sim_blocks <- if (quick) 6000L else 12000L
substeps_per_block <- 8L
max_qbar_samples_per_pair <- if (quick) 60L else 80L
if (recommended) {
  n_states <- 21L
  transition_sim_blocks <- 210000L
  substeps_per_block <- 32L
  max_qbar_samples_per_pair <- 7L
}
n_states <- as.integer(option_value("n-states", n_states))
transition_sim_blocks <- as.integer(option_value(
  "transition-sim-blocks", transition_sim_blocks
))
substeps_per_block <- as.integer(option_value("substeps", substeps_per_block))
max_qbar_samples_per_pair <- as.integer(option_value(
  "max-q-samples", max_qbar_samples_per_pair
))
kernel_seed <- as.integer(option_value("kernel-seed", 5000L))
if (!is.finite(kernel_seed)) stop("--kernel-seed must be an integer.")
n_steps_override <- suppressWarnings(as.integer(option_value("n-steps", NA)))
if (is.finite(n_steps_override) && n_steps_override > 0L) {
  cases$n_steps <- n_steps_override
}
terminal_override <- suppressWarnings(as.numeric(option_value("terminal", NA)))
if (is.finite(terminal_override) && terminal_override > 0) {
  cases$terminal <- terminal_override
}
block_size_override <- suppressWarnings(as.integer(option_value(
  "block-size", NA
)))
block_horizon <- suppressWarnings(as.numeric(option_value(
  "block-horizon", if (recommended) 1 / 20 else NA
)))
if (is.finite(block_size_override) && block_size_override >= 10L) {
  cases$block_size <- block_size_override
} else if (is.finite(block_horizon) && block_horizon > 0) {
  cases$block_size <- pmax(10L, as.integer(round(
    block_horizon * cases$n_steps / cases$terminal
  )))
}
seed_override <- option_value("seeds", NULL)
if (!is.null(seed_override)) {
  seed_values <- suppressWarnings(as.integer(strsplit(seed_override, ",",
                                                       fixed = TRUE)[[1L]]))
  seed_values <- seed_values[is.finite(seed_values)]
  if (!length(seed_values)) stop("--seeds must contain at least one integer.")
  case_shapes <- unique(cases[, c("n_steps", "terminal", "block_size"),
                              drop = FALSE])
  cases <- merge(data.frame(seed = seed_values), case_shapes)
  cases <- cases[, c("seed", "n_steps", "terminal", "block_size"),
                 drop = FALSE]
}

tests <- c("log", "inv1", "inv4", "ratio1", "exp025")
instruments <- c(
  "a1", "center_log", "bounded_poly1", "bounded_poly2", "tanh_log",
  "logistic_low", "logistic_high", "low_bin", "mid_bin", "high_bin"
)
lags <- c(1L, 2L)

bind_rows_fill <- function(rows) {
  rows <- rows[lengths(rows) > 0L]
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  filled <- lapply(rows, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, filled)
}

make_endpoint_breaks <- function(n_states, seed = 881L) {
  v <- stationary_v_sample(300000L, seed = seed)
  probs <- seq(0, 1, length.out = n_states + 1L)
  breaks <- stats::quantile(v, probs = probs, names = FALSE, type = 8)
  breaks <- cummax(breaks)
  for (i in 2:length(breaks)) {
    if (breaks[i] <= breaks[i - 1L]) {
      breaks[i] <- breaks[i - 1L] + .Machine$double.eps
    }
  }
  cut_breaks <- breaks
  cut_breaks[1L] <- -Inf
  cut_breaks[length(cut_breaks)] <- Inf
  centers <- numeric(n_states)
  bin <- cut(v, breaks = cut_breaks, labels = FALSE, right = TRUE)
  for (s in seq_len(n_states)) {
    centers[s] <- stats::median(v[bin == s], na.rm = TRUE)
  }
  list(breaks = breaks, cut_breaks = cut_breaks, centers = centers)
}

make_q_endpoint_breaks <- function(n_states, seed = 881L) {
  q <- stationary_q_sample(300000L, seed = seed)
  probs <- seq(0, 1, length.out = n_states + 1L)
  breaks <- stats::quantile(q, probs = probs, names = FALSE, type = 8)
  breaks <- cummax(breaks)
  for (i in 2:length(breaks)) {
    if (breaks[i] <= breaks[i - 1L]) {
      breaks[i] <- breaks[i - 1L] + .Machine$double.eps
    }
  }
  cut_breaks <- breaks
  cut_breaks[1L] <- 1
  cut_breaks[length(cut_breaks)] <- Inf
  bin <- cut(q, breaks = cut_breaks, labels = FALSE, right = TRUE,
             include.lowest = TRUE)
  centers <- vapply(seq_len(n_states), function(s) {
    stats::median(q[bin == s], na.rm = TRUE)
  }, numeric(1L))
  list(breaks = breaks, cut_breaks = cut_breaks, centers = centers,
       endpoint_state = "q")
}

state_index <- function(v, cut_breaks) {
  as.integer(cut(v, breaks = cut_breaks, labels = FALSE, right = TRUE,
                 include.lowest = TRUE))
}

block_truth <- function(V, block_size, n_blocks, state_info,
                        endpoint_state = "v") {
  q_path <- 1 + V[seq_len(n_blocks * block_size)] ^ 2
  q_bar <- colMeans(matrix(q_path, nrow = block_size, ncol = n_blocks))
  endpoint_idx <- 1L + (0:n_blocks) * block_size
  endpoints <- if (identical(endpoint_state, "q")) {
    1 + V[endpoint_idx] ^ 2
  } else {
    V[endpoint_idx]
  }
  states <- state_index(endpoints, state_info$cut_breaks)
  list(
    q_bar = q_bar,
    endpoints = endpoints,
    start_state = states[-length(states)],
    end_state = states[-1L]
  )
}

simulate_endpoint_additive_model <- function(rho, block_h, state_info,
                                             n_sim_blocks,
                                             substeps_per_block = 8L,
                                             max_qbar_samples_per_pair = 80L,
                                             initial_v = 0,
                                             state_model = c("conditional",
                                                             "collocation"),
                                             seed = 1L) {
  set.seed(seed)
  state_model <- match.arg(state_model)
  n_states <- length(state_info$centers)
  dt_sub <- block_h / substeps_per_block
  rho2 <- rho ^ 2
  n_paths_per_state <- ceiling(n_sim_blocks / n_states)

  # This is a finite-state Feynman-Kac approximation.  Conditional paths are
  # sampled separately from every starting state so all rows of the endpoint
  # kernel and their additive-functional laws receive comparable Monte Carlo
  # resolution, including the tails.  The collocation version starts at a
  # definite grid value, whereas the conditional version averages over the
  # stationary law inside a cell and is retained for the approximation study.
  simulate_block_paths <- function(v_start) {
    v <- as.numeric(v_start)
    q_sum <- numeric(length(v))
    for (s in seq_len(substeps_per_block)) {
      q_now <- 1 + v ^ 2
      q_sum <- q_sum + q_now
      v <- v + (-rho2 * v * dt_sub) +
        rho * sqrt(pmax(q_now, 1e-12)) * sqrt(dt_sub) *
        stats::rnorm(length(v))
    }
    list(end_state = state_index(v, state_info$cut_breaks),
         q_bar = q_sum / substeps_per_block)
  }

  if (identical(state_model, "conditional")) {
    reference_v <- stationary_v_sample(max(300000L, 2000L * n_states),
                                       seed = seed + 17L)
    reference_state <- state_index(reference_v, state_info$cut_breaks)
  }
  counts <- matrix(0, n_states, n_states)
  q_samples <- vector("list", n_states * n_states)
  pair_index <- function(a, b) (a - 1L) * n_states + b
  for (a in seq_len(n_states)) {
    starts <- if (identical(state_model, "conditional")) {
      candidates <- reference_v[reference_state == a]
      if (!length(candidates)) stop("Endpoint grid has an empty starting bin.")
      sample(candidates, n_paths_per_state, replace = TRUE)
    } else {
      rep(state_info$centers[a], n_paths_per_state)
    }
    paths <- simulate_block_paths(starts)
    counts[a, ] <- tabulate(paths$end_state, nbins = n_states)
    for (b in which(counts[a, ] > 0)) {
      q_samples[[pair_index(a, b)]] <- paths$q_bar[paths$end_state == b]
    }
  }

  transition <- counts / n_paths_per_state
  active_end_states <- lapply(seq_len(n_states), function(a) {
    which(counts[a, ] > 0)
  })
  initial_paths <- simulate_block_paths(rep(initial_v, n_paths_per_state))
  initial_counts <- tabulate(initial_paths$end_state, nbins = n_states)
  initial_transition <- initial_counts / n_paths_per_state
  initial_q_samples <- vector("list", n_states)
  for (b in which(initial_counts > 0)) {
    initial_q_samples[[b]] <- initial_paths$q_bar[
      initial_paths$end_state == b
    ]
  }

  compress_samples <- function(qs) {
    qs <- qs[is.finite(qs) & qs > 0]
    if (length(qs) > max_qbar_samples_per_pair) {
      probs <- (seq_len(max_qbar_samples_per_pair) - 0.5) /
        max_qbar_samples_per_pair
      qs <- as.numeric(stats::quantile(qs, probs = probs, names = FALSE,
                                       type = 8))
    }
    qs
  }
  for (idx in seq_along(q_samples)) {
    q_samples[idx] <- list(compress_samples(q_samples[[idx]]))
  }
  for (b in seq_along(initial_q_samples)) {
    initial_q_samples[b] <- list(compress_samples(initial_q_samples[[b]]))
  }
  active_pairs <- which(counts > 0, arr.ind = TRUE)
  pair_node_count <- vapply(seq_len(nrow(active_pairs)), function(i) {
    length(q_samples[[pair_index(active_pairs[i, 1L], active_pairs[i, 2L])]])
  }, integer(1L))
  max_pair_nodes <- max(pair_node_count)
  pair_q_nodes <- matrix(NA_real_, nrow = nrow(active_pairs),
                         ncol = max_pair_nodes)
  for (i in seq_len(nrow(active_pairs))) {
    qs <- q_samples[[pair_index(active_pairs[i, 1L], active_pairs[i, 2L])]]
    pair_q_nodes[i, seq_along(qs)] <- qs
  }
  list(
    rho = rho,
    transition = transition,
    active_end_states = active_end_states,
    active_pairs = active_pairs,
    pair_q_nodes = pair_q_nodes,
    pair_node_count = pair_node_count,
    init_prob = rep(1 / n_states, n_states),
    initial_v = initial_v,
    state_model = state_model,
    initial_transition = initial_transition,
    initial_q_samples = initial_q_samples,
    q_samples = q_samples,
    counts = counts,
    state_info = state_info,
    n_states = n_states,
    block_h = block_h,
    n_sim_blocks = n_sim_blocks,
    n_paths_per_state = n_paths_per_state,
    substeps_per_block = substeps_per_block
  )
}

simulate_q_endpoint_additive_model <- function(
    rho, block_h, state_info, n_sim_blocks,
    substeps_per_block = 32L, max_qbar_samples_per_pair = 7L,
    initial_q = 1, state_model = c("conditional", "collocation"),
    seed = 1L) {
  set.seed(seed)
  state_model <- match.arg(state_model)
  n_states <- length(state_info$centers)
  dt_sub <- block_h / substeps_per_block
  rho2 <- rho ^ 2
  n_paths_per_state <- max(100L, as.integer(ceiling(n_sim_blocks / n_states)))
  stationary_pool <- stationary_q_sample(max(300000L, 20L * n_paths_per_state),
                                         seed = seed + 103L)
  stationary_bin <- state_index(stationary_pool, state_info$cut_breaks)
  counts <- matrix(0, n_states, n_states)
  q_samples <- vector("list", n_states * n_states)

  draw_start_q <- function(s, n) {
    if (identical(state_model, "collocation")) {
      return(rep(state_info$centers[s], n))
    }
    pool <- stationary_pool[stationary_bin == s]
    if (!length(pool)) return(rep(state_info$centers[s], n))
    sample(pool, n, replace = length(pool) < n)
  }
  simulate_paths <- function(q_start) {
    n <- length(q_start)
    signs <- sample(c(-1, 1), n, replace = TRUE)
    v <- signs * sqrt(pmax(q_start - 1, 0))
    q_acc <- numeric(n)
    for (step in seq_len(substeps_per_block)) {
      q_acc <- q_acc + (1 + v ^ 2)
      d_w <- stats::rnorm(n, sd = sqrt(dt_sub))
      v <- v - rho2 * v * dt_sub +
        rho * sqrt(pmax(1 + v ^ 2, 1e-12)) * d_w
    }
    list(q_end = 1 + v ^ 2, q_bar = q_acc / substeps_per_block)
  }
  for (a in seq_len(n_states)) {
    paths <- simulate_paths(draw_start_q(a, n_paths_per_state))
    end_state <- state_index(paths$q_end, state_info$cut_breaks)
    for (b in seq_len(n_states)) {
      q <- paths$q_bar[end_state == b]
      counts[a, b] <- length(q)
      if (length(q)) q_samples[[pair_index(a, b, n_states)]] <- q
    }
  }
  transition <- counts / rowSums(counts)

  initial_paths <- simulate_paths(rep(initial_q, n_paths_per_state))
  initial_end <- state_index(initial_paths$q_end, state_info$cut_breaks)
  initial_counts <- tabulate(initial_end, nbins = n_states)
  initial_transition <- initial_counts / sum(initial_counts)
  initial_q_samples <- lapply(seq_len(n_states), function(b) {
    initial_paths$q_bar[initial_end == b]
  })

  compress <- function(q) {
    q <- q[is.finite(q) & q >= 1]
    if (length(q) > max_qbar_samples_per_pair) {
      probs <- (seq_len(max_qbar_samples_per_pair) - 0.5) /
        max_qbar_samples_per_pair
      q <- as.numeric(stats::quantile(q, probs = probs, names = FALSE,
                                      type = 8))
    }
    q
  }
  q_samples <- lapply(q_samples, compress)
  initial_q_samples <- lapply(initial_q_samples, compress)
  active_pairs <- which(counts > 0, arr.ind = TRUE)
  pair_node_count <- vapply(seq_len(nrow(active_pairs)), function(i) {
    length(q_samples[[pair_index(active_pairs[i, 1L], active_pairs[i, 2L],
                                  n_states)]])
  }, integer(1L))
  pair_q_nodes <- matrix(NA_real_, nrow = nrow(active_pairs),
                         ncol = max(pair_node_count))
  for (i in seq_len(nrow(active_pairs))) {
    q <- q_samples[[pair_index(active_pairs[i, 1L], active_pairs[i, 2L],
                                n_states)]]
    pair_q_nodes[i, seq_along(q)] <- q
  }
  list(
    rho = rho,
    transition = transition,
    active_end_states = lapply(seq_len(n_states), function(a) {
      which(counts[a, ] > 0)
    }),
    active_pairs = active_pairs,
    pair_q_nodes = pair_q_nodes,
    pair_node_count = pair_node_count,
    init_prob = rep(1 / n_states, n_states),
    initial_v = 0,
    state_model = paste0("q_", state_model),
    initial_transition = initial_transition,
    initial_q_samples = initial_q_samples,
    q_samples = q_samples,
    counts = counts,
    state_info = state_info,
    n_states = n_states,
    block_h = block_h,
    n_sim_blocks = n_sim_blocks,
    n_paths_per_state = n_paths_per_state,
    substeps_per_block = substeps_per_block,
    endpoint_state = "q"
  )
}

activity_quadrature_nodes <- function(mean_q, var_q, n_nodes = 7L) {
  n_nodes <- max(1L, as.integer(n_nodes))
  if (!is.finite(mean_q) || mean_q < 1) mean_q <- 1
  if (!is.finite(var_q) || var_q <= 1e-12 || n_nodes == 1L) return(mean_q)
  z <- stats::qnorm((seq_len(n_nodes) - 0.5) / n_nodes)
  z <- z / sqrt(mean(z ^ 2))
  positive_mean <- mean_q - 1
  if (positive_mean > 1e-8) {
    log_var <- log1p(var_q / positive_mean ^ 2)
    log_mean <- log(positive_mean) - 0.5 * log_var
    return(1 + exp(log_mean + sqrt(log_var) * z))
  }
  pmax(1, mean_q + sqrt(var_q) * z)
}

build_feynman_kac_endpoint_model <- function(
    rho, block_h, state_info, substeps_per_block = 32L,
    q_nodes_per_pair = 7L) {
  n_states <- length(state_info$centers)
  n_substeps <- max(1L, as.integer(substeps_per_block))
  delta <- block_h / n_substeps
  z_centers <- asinh(state_info$centers)
  z_breaks <- asinh(state_info$cut_breaks)
  q_centers <- cosh(z_centers) ^ 2
  sd_step <- rho * sqrt(delta)
  step_transition <- matrix(0, n_states, n_states)
  for (a in seq_len(n_states)) {
    mean_step <- z_centers[a] - 1.5 * rho ^ 2 *
      tanh(z_centers[a]) * delta
    upper <- stats::pnorm((z_breaks[-1L] - mean_step) / sd_step)
    lower <- stats::pnorm((z_breaks[-length(z_breaks)] - mean_step) / sd_step)
    step_transition[a, ] <- pmax(upper - lower, .Machine$double.xmin)
    step_transition[a, ] <- step_transition[a, ] /
      sum(step_transition[a, ])
  }
  step_cost <- delta * outer(q_centers, q_centers, "+") / 2
  kernel0 <- step_transition
  kernel1 <- step_transition * step_cost
  kernel2 <- step_transition * step_cost ^ 2

  transition <- diag(n_states)
  additive1 <- matrix(0, n_states, n_states)
  additive2 <- matrix(0, n_states, n_states)
  for (step in seq_len(n_substeps)) {
    next_transition <- transition %*% kernel0
    next_additive1 <- additive1 %*% kernel0 + transition %*% kernel1
    next_additive2 <- additive2 %*% kernel0 +
      2 * additive1 %*% kernel1 + transition %*% kernel2
    transition <- next_transition
    additive1 <- next_additive1
    additive2 <- next_additive2
  }
  transition <- transition / rowSums(transition)
  moment_transition <- transition
  mean_additive <- additive1 / moment_transition
  var_additive <- pmax(additive2 / moment_transition - mean_additive ^ 2, 0)
  mean_qbar <- mean_additive / block_h
  var_qbar <- var_additive / block_h ^ 2

  # Local-linear Gaussian transition in Lamperti coordinates avoids the
  # artificial support gaps of a simulated histogram.  The substep chain
  # above is retained only for conditional additive-functional moments.
  transition <- matrix(0, n_states, n_states)
  for (a in seq_len(n_states)) {
    lambda <- 1.5 * rho ^ 2 / cosh(z_centers[a]) ^ 2
    drift <- -1.5 * rho ^ 2 * tanh(z_centers[a])
    mean_h <- z_centers[a] + drift * (-expm1(-lambda * block_h)) / lambda
    var_h <- rho ^ 2 * (-expm1(-2 * lambda * block_h)) / (2 * lambda)
    upper <- stats::pnorm((z_breaks[-1L] - mean_h) / sqrt(var_h))
    lower <- stats::pnorm((z_breaks[-length(z_breaks)] - mean_h) /
                            sqrt(var_h))
    transition[a, ] <- pmax(upper - lower, .Machine$double.xmin)
    transition[a, ] <- transition[a, ] / sum(transition[a, ])
  }

  initial_v <- 0
  initial_z <- asinh(initial_v)
  initial_q <- 1 + initial_v ^ 2
  initial_mean <- initial_z - 1.5 * rho ^ 2 * tanh(initial_z) * delta
  initial_transition <- pmax(
    stats::pnorm((z_breaks[-1L] - initial_mean) / sd_step) -
      stats::pnorm((z_breaks[-length(z_breaks)] - initial_mean) / sd_step),
    .Machine$double.xmin
  )
  initial_transition <- initial_transition / sum(initial_transition)
  initial_cost <- delta * (initial_q + q_centers) / 2
  initial_additive1 <- initial_transition * initial_cost
  initial_additive2 <- initial_transition * initial_cost ^ 2
  if (n_substeps >= 2L) {
    for (step in 2:n_substeps) {
      next_transition <- as.numeric(initial_transition %*% kernel0)
      next_additive1 <- as.numeric(initial_additive1 %*% kernel0 +
                                     initial_transition %*% kernel1)
      next_additive2 <- as.numeric(initial_additive2 %*% kernel0 +
                                     2 * initial_additive1 %*% kernel1 +
                                     initial_transition %*% kernel2)
      initial_transition <- next_transition
      initial_additive1 <- next_additive1
      initial_additive2 <- next_additive2
    }
  }
  initial_transition <- initial_transition / sum(initial_transition)
  initial_moment_transition <- initial_transition
  initial_mean_additive <- initial_additive1 / initial_moment_transition
  initial_var_additive <- pmax(
    initial_additive2 / initial_moment_transition - initial_mean_additive ^ 2, 0
  )
  initial_lambda <- 1.5 * rho ^ 2 / cosh(initial_z) ^ 2
  initial_drift <- -1.5 * rho ^ 2 * tanh(initial_z)
  initial_mean_h <- initial_z + initial_drift *
    (-expm1(-initial_lambda * block_h)) / initial_lambda
  initial_var_h <- rho ^ 2 * (-expm1(-2 * initial_lambda * block_h)) /
    (2 * initial_lambda)
  initial_transition <- pmax(
    stats::pnorm((z_breaks[-1L] - initial_mean_h) / sqrt(initial_var_h)) -
      stats::pnorm((z_breaks[-length(z_breaks)] - initial_mean_h) /
                     sqrt(initial_var_h)),
    .Machine$double.xmin
  )
  initial_transition <- initial_transition / sum(initial_transition)

  q_samples <- vector("list", n_states * n_states)
  for (a in seq_len(n_states)) {
    for (b in seq_len(n_states)) {
      idx <- pair_index(a, b, n_states)
      nodes <- activity_quadrature_nodes(mean_qbar[a, b], var_qbar[a, b],
                                         q_nodes_per_pair)
      if (length(nodes) < q_nodes_per_pair) nodes <- rep(nodes, q_nodes_per_pair)
      q_samples[[idx]] <- nodes
    }
  }
  initial_q_samples <- lapply(seq_len(n_states), function(b) {
    activity_quadrature_nodes(
      initial_mean_additive[b] / block_h,
      initial_var_additive[b] / block_h ^ 2,
      q_nodes_per_pair
    )
  })
  active_pairs <- which(transition > 0, arr.ind = TRUE)
  pair_q_nodes <- matrix(NA_real_, nrow = nrow(active_pairs),
                         ncol = q_nodes_per_pair)
  for (i in seq_len(nrow(active_pairs))) {
    pair_q_nodes[i, ] <- q_samples[[pair_index(active_pairs[i, 1L],
                                                active_pairs[i, 2L],
                                                n_states)]]
  }
  active_end_states <- lapply(seq_len(n_states), function(a) {
    which(transition[a, ] > 0)
  })
  list(
    rho = rho,
    transition = transition,
    active_end_states = active_end_states,
    active_pairs = active_pairs,
    pair_q_nodes = pair_q_nodes,
    pair_node_count = rep(as.integer(q_nodes_per_pair), nrow(active_pairs)),
    init_prob = rep(1 / n_states, n_states),
    initial_v = initial_v,
    state_model = "feynman_kac_moment",
    initial_transition = initial_transition,
    initial_q_samples = initial_q_samples,
    q_samples = q_samples,
    counts = transition,
    state_info = state_info,
    n_states = n_states,
    block_h = block_h,
    n_sim_blocks = NA_integer_,
    n_paths_per_state = NA_integer_,
    substeps_per_block = n_substeps,
    additive_mean = mean_qbar,
    additive_variance = var_qbar
  )
}

build_hybrid_endpoint_additive_model <- function(
    rho, block_h, state_info, n_sim_blocks = 210000L,
    substeps_per_block = 32L, q_nodes_per_pair = 7L,
    seed = 5000L) {
  analytic <- build_feynman_kac_endpoint_model(
    rho = rho,
    block_h = block_h,
    state_info = state_info,
    substeps_per_block = substeps_per_block,
    q_nodes_per_pair = q_nodes_per_pair
  )
  additive_mc <- simulate_endpoint_additive_model(
    rho = rho,
    block_h = block_h,
    state_info = state_info,
    n_sim_blocks = n_sim_blocks,
    substeps_per_block = substeps_per_block,
    max_qbar_samples_per_pair = q_nodes_per_pair,
    initial_v = 0,
    state_model = "collocation",
    seed = seed
  )
  n_states <- analytic$n_states
  q_samples <- vector("list", n_states * n_states)
  missing_pair <- matrix(FALSE, n_states, n_states)
  for (a in seq_len(n_states)) {
    for (b in seq_len(n_states)) {
      idx <- pair_index(a, b, n_states)
      q <- additive_mc$q_samples[[idx]]
      if (!length(q)) {
        q <- analytic$q_samples[[idx]]
        missing_pair[a, b] <- TRUE
      }
      if (length(q) < q_nodes_per_pair) {
        q <- q[as.integer(round(seq(1, length(q),
                                  length.out = q_nodes_per_pair)))]
      }
      q_samples[[idx]] <- q
    }
  }
  initial_q_samples <- lapply(seq_len(n_states), function(b) {
    q <- additive_mc$initial_q_samples[[b]]
    if (length(q)) q else analytic$initial_q_samples[[b]]
  })
  active_pairs <- which(analytic$transition > 0, arr.ind = TRUE)
  max_nodes <- max(vapply(q_samples, length, integer(1L)))
  pair_q_nodes <- matrix(NA_real_, nrow = nrow(active_pairs),
                         ncol = max_nodes)
  pair_node_count <- integer(nrow(active_pairs))
  for (i in seq_len(nrow(active_pairs))) {
    q <- q_samples[[pair_index(active_pairs[i, 1L], active_pairs[i, 2L],
                                n_states)]]
    pair_q_nodes[i, seq_along(q)] <- q
    pair_node_count[i] <- length(q)
  }
  analytic$q_samples <- q_samples
  analytic$initial_q_samples <- initial_q_samples
  analytic$active_pairs <- active_pairs
  analytic$active_end_states <- lapply(seq_len(n_states), function(a) {
    which(analytic$transition[a, ] > 0)
  })
  analytic$pair_q_nodes <- pair_q_nodes
  analytic$pair_node_count <- pair_node_count
  analytic$counts <- additive_mc$counts
  analytic$n_sim_blocks <- n_sim_blocks
  analytic$n_paths_per_state <- additive_mc$n_paths_per_state
  analytic$state_model <- "hybrid_analytic_transition_mc_additive"
  analytic$additive_fallback_pair_count <- sum(missing_pair)
  analytic$additive_fallback_transition_mass <-
    mean(rowSums(analytic$transition * missing_pair))
  analytic
}

simulate_y_from_qbar <- function(q_bar, ecf_obs, sigma1, rho, seed = 1L) {
  set.seed(seed)
  n_blocks <- min(length(q_bar), nrow(ecf_obs$y_dejumped))
  freqs <- ecf_obs$freqs
  y <- matrix(NA_real_, nrow = n_blocks, ncol = length(freqs))
  gamma <- (sigma1 * rho) ^ 2
  for (j in seq_len(n_blocks)) {
    mu <- -0.5 * freqs ^ 2 * gamma * q_bar[j]
    R <- regularized_chol(ecf_obs$obs_cov[j, , ])
    y[j, ] <- as.numeric(mu + t(R) %*% stats::rnorm(length(freqs)))
  }
  y
}

emission_terms <- function(y, S, freqs,
                           emission_mode = c("full_gaussian",
                                             "projected_quasi",
                                             "projected_diag")) {
  emission_mode <- match.arg(emission_mode)
  if (identical(emission_mode, "projected_diag")) {
    S <- diag(pmax(diag(as.matrix(S)), .Machine$double.eps), nrow(S))
  }
  R <- regularized_chol(S)
  inv_y <- backsolve(R, y, transpose = TRUE)
  c0 <- -0.5 * freqs ^ 2
  inv_c <- backsolve(R, c0, transpose = TRUE)
  if (identical(emission_mode, "projected_quasi") ||
      identical(emission_mode, "projected_diag")) {
    denom <- sum(inv_c ^ 2)
    if (!is.finite(denom) || denom <= 0) {
      return(list(mode = "projected_quasi", t = NA_real_, var = NA_real_))
    }
    return(list(
      mode = "projected_quasi",
      t = sum(inv_c * inv_y) / denom,
      var = 1 / denom
    ))
  }
  list(
    mode = "full_gaussian",
    ySiY = sum(inv_y ^ 2),
    cSiY = sum(inv_c * inv_y),
    cSiC = sum(inv_c ^ 2),
    logdet = 2 * sum(log(diag(R))),
    p = length(freqs)
  )
}

log_emission_q <- function(terms, gamma, q, dispersion = 1) {
  if (!is.finite(dispersion) || dispersion <= 0) {
    return(rep(-Inf, length(q)))
  }
  if (identical(terms$mode, "projected_quasi")) {
    if (!is.finite(terms$t) || !is.finite(terms$var) || terms$var <= 0) {
      return(rep(-Inf, length(q)))
    }
    var <- dispersion * terms$var
    return(-0.5 * (log(2 * pi * var) +
                     (terms$t - gamma * q) ^ 2 / var))
  }
  qf <- terms$ySiY - 2 * gamma * q * terms$cSiY +
    gamma ^ 2 * q ^ 2 * terms$cSiC
  -0.5 * (terms$p * log(2 * pi) + terms$logdet + qf)
}

logmeanexp <- function(x) log_sum_exp(x) - log(length(x))

pair_index <- function(a, b, n_states) (a - 1L) * n_states + b

initial_log_kernel <- function(terms, gamma, model, dispersion = 1) {
  out <- rep(-Inf, model$n_states)
  for (b in seq_len(model$n_states)) {
    p_b <- model$initial_transition[b]
    qs <- model$initial_q_samples[[b]]
    if (!is.finite(p_b) || p_b <= 0 || !length(qs)) next
    out[b] <- log(p_b) + logmeanexp(log_emission_q(
      terms, gamma, qs, dispersion = dispersion
    ))
  }
  out
}

pair_log_emissions <- function(terms, gamma, model, dispersion = 1) {
  if (!identical(terms$mode, "projected_quasi")) return(NULL)
  if (!is.finite(dispersion) || dispersion <= 0) return(NULL)
  var <- dispersion * terms$var
  q <- model$pair_q_nodes
  z <- -0.5 * (log(2 * pi * var) +
                 (terms$t - gamma * q) ^ 2 / var)
  row_max <- apply(z, 1L, max, na.rm = TRUE)
  row_max + log(rowSums(exp(z - row_max), na.rm = TRUE) /
                  model$pair_node_count)
}

known_endpoint_oracle_qbar_loglik <- function(y_mat, obs_cov, freqs,
                                              start_state, end_state,
                                              q_bar, model, sigma1, rho,
                                              emission_mode = "full_gaussian",
                                              emission_dispersion = 1,
                                              terms_list = NULL) {
  gamma <- (sigma1 * rho) ^ 2
  loglik <- 0
  for (j in seq_along(q_bar)) {
    a <- start_state[j]
    b <- end_state[j]
    if (!is.finite(a) || !is.finite(b) || a < 1L || b < 1L ||
        a > model$n_states || b > model$n_states) return(-Inf)
    p_ab <- if (j == 1L) model$initial_transition[b] else model$transition[a, b]
    if (!is.finite(p_ab) || p_ab <= 0) return(-Inf)
    terms <- if (is.null(terms_list)) {
      emission_terms(y_mat[j, ], obs_cov[j, , ], freqs,
                     emission_mode = emission_mode)
    } else {
      terms_list[[j]]
    }
    loglik <- loglik + log(p_ab) + log_emission_q(
      terms, gamma, q_bar[j], dispersion = emission_dispersion
    )
  }
  loglik
}

known_endpoint_transition_loglik <- function(start_state, end_state, model) {
  loglik <- 0
  for (j in seq_along(start_state)) {
    a <- start_state[j]
    b <- end_state[j]
    if (!is.finite(b) || b < 1L || b > model$n_states) return(-Inf)
    if (j > 1L && (!is.finite(a) || a < 1L || a > model$n_states)) {
      return(-Inf)
    }
    p_ab <- if (j == 1L) model$initial_transition[b] else model$transition[a, b]
    if (!is.finite(p_ab) || p_ab <= 0) return(-Inf)
    loglik <- loglik + log(p_ab)
  }
  loglik
}

known_endpoint_integrated_qbar_loglik <- function(y_mat, obs_cov, freqs,
                                                  start_state, end_state,
                                                  model, sigma1, rho,
                                                  emission_mode = "full_gaussian",
                                                  emission_dispersion = 1,
                                                  terms_list = NULL) {
  gamma <- (sigma1 * rho) ^ 2
  n_states <- model$n_states
  loglik <- 0
  for (j in seq_along(start_state)) {
    a <- start_state[j]
    b <- end_state[j]
    if (!is.finite(a) || !is.finite(b) || a < 1L || b < 1L ||
        a > n_states || b > n_states) return(-Inf)
    p_ab <- if (j == 1L) model$initial_transition[b] else model$transition[a, b]
    if (!is.finite(p_ab) || p_ab <= 0) return(-Inf)
    qs <- if (j == 1L) model$initial_q_samples[[b]] else {
      model$q_samples[[pair_index(a, b, n_states)]]
    }
    if (!length(qs)) return(-Inf)
    terms <- if (is.null(terms_list)) {
      emission_terms(y_mat[j, ], obs_cov[j, , ], freqs,
                     emission_mode = emission_mode)
    } else {
      terms_list[[j]]
    }
    loglik <- loglik + log(p_ab) +
      logmeanexp(log_emission_q(terms, gamma, qs,
                                dispersion = emission_dispersion))
  }
  loglik
}

latent_endpoint_hmm_loglik <- function(y_mat, obs_cov, freqs, model,
                                       sigma1, rho,
                                       emission_mode = "full_gaussian",
                                       emission_dispersion = 1,
                                       terms_list = NULL,
                                       return_filter = FALSE) {
  gamma <- (sigma1 * rho) ^ 2
  n_states <- model$n_states
  logP <- log(model$transition)
  logP[!is.finite(logP)] <- -Inf
  log_alpha <- NULL
  loglik <- 0
  filter_prob <- if (return_filter) {
    matrix(NA_real_, nrow = nrow(y_mat), ncol = n_states)
  } else {
    NULL
  }
  loglik_increment <- if (return_filter) numeric(nrow(y_mat)) else NULL
  for (j in seq_len(nrow(y_mat))) {
    terms <- if (is.null(terms_list)) {
      emission_terms(y_mat[j, ], obs_cov[j, , ], freqs,
                     emission_mode = emission_mode)
    } else {
      terms_list[[j]]
    }
    if (j == 1L) {
      next_alpha <- initial_log_kernel(terms, gamma, model,
                                       dispersion = emission_dispersion)
    } else {
      logE <- matrix(-Inf, n_states, n_states)
      pair_terms <- pair_log_emissions(terms, gamma, model,
                                       dispersion = emission_dispersion)
      if (!is.null(pair_terms)) {
        logE[model$active_pairs] <- pair_terms
      } else {
        for (a in seq_len(n_states)) {
          for (b in model$active_end_states[[a]]) {
            qs <- model$q_samples[[pair_index(a, b, n_states)]]
            if (!length(qs)) next
            logE[a, b] <- logmeanexp(log_emission_q(
              terms, gamma, qs, dispersion = emission_dispersion
            ))
          }
        }
      }
      next_alpha <- rep(-Inf, n_states)
      for (b in seq_len(n_states)) {
        next_alpha[b] <- log_sum_exp(log_alpha + logP[, b] + logE[, b])
      }
    }
    inc <- log_sum_exp(next_alpha)
    if (!is.finite(inc)) return(-Inf)
    loglik <- loglik + inc
    log_alpha <- next_alpha - inc
    if (return_filter) {
      filter_prob[j, ] <- exp(log_alpha)
      loglik_increment[j] <- inc
    }
  }
  if (!return_filter) return(loglik)
  entropy <- -rowSums(filter_prob * log(pmax(filter_prob,
                                             .Machine$double.xmin)))
  list(
    loglik = loglik,
    filter_prob = filter_prob,
    loglik_increment = loglik_increment,
    mean_filter_max_prob = mean(apply(filter_prob, 1L, max)),
    mean_filter_entropy = mean(entropy) / log(n_states),
    mean_filter_boundary_prob = mean(filter_prob[, 1L] +
                                       filter_prob[, n_states])
  )
}

pair_emission_moments <- function(terms, gamma, model, dispersion) {
  var <- dispersion * terms$var
  q <- model$pair_q_nodes
  z <- -0.5 * (log(2 * pi * var) + (terms$t - gamma * q) ^ 2 / var)
  z[!is.finite(q)] <- -Inf
  row_max <- apply(z, 1L, max)
  w <- exp(z - row_max)
  w[!is.finite(w)] <- 0
  sum_w <- rowSums(w)
  list(
    log_mean = row_max + log(sum_w / model$pair_node_count),
    q_mean = rowSums(w * q, na.rm = TRUE) / sum_w,
    q2_mean = rowSums(w * q ^ 2, na.rm = TRUE) / sum_w
  )
}

initial_emission_moments <- function(terms, gamma, model, dispersion) {
  out <- matrix(NA_real_, nrow = model$n_states, ncol = 3L,
                dimnames = list(NULL, c("log_mean", "q_mean", "q2_mean")))
  for (b in seq_len(model$n_states)) {
    q <- model$initial_q_samples[[b]]
    if (!length(q)) next
    z <- log_emission_q(terms, gamma, q, dispersion = dispersion)
    z_max <- max(z)
    w <- exp(z - z_max)
    sum_w <- sum(w)
    out[b, ] <- c(
      z_max + log(sum_w / length(q)),
      sum(w * q) / sum_w,
      sum(w * q ^ 2) / sum_w
    )
  }
  out
}

fit_latent_emission_em <- function(terms_list, model, beta_initial,
                                   dispersion_initial = 1,
                                   beta_bounds = c(0.20, 5.00),
                                   dispersion_bounds = c(0.25, Inf),
                                   max_iter = 100L, tol = 1e-4) {
  n_obs <- length(terms_list)
  n_states <- model$n_states
  gamma <- min(max(beta_initial ^ 2, beta_bounds[1L] ^ 2),
               beta_bounds[2L] ^ 2)
  dispersion <- min(max(dispersion_initial, dispersion_bounds[1L]),
                    dispersion_bounds[2L])
  logP <- log(model$transition)
  logP[!is.finite(logP)] <- -Inf

  e_step <- function(gamma, dispersion, retain = TRUE) {
    log_alpha <- matrix(NA_real_, nrow = n_obs, ncol = n_states)
    increments <- numeric(n_obs)
    logK <- vector("list", n_obs)
    q_mean <- vector("list", n_obs)
    q2_mean <- vector("list", n_obs)

    init <- initial_emission_moments(terms_list[[1L]], gamma, model,
                                     dispersion)
    logK[[1L]] <- log(model$initial_transition) + init[, "log_mean"]
    q_mean[[1L]] <- init[, "q_mean"]
    q2_mean[[1L]] <- init[, "q2_mean"]
    increments[1L] <- log_sum_exp(logK[[1L]])
    log_alpha[1L, ] <- logK[[1L]] - increments[1L]

    if (n_obs >= 2L) {
      for (j in 2:n_obs) {
        moments <- pair_emission_moments(terms_list[[j]], gamma, model,
                                         dispersion)
        edge_logK <- matrix(-Inf, n_states, n_states)
        edge_q <- matrix(NA_real_, n_states, n_states)
        edge_q2 <- matrix(NA_real_, n_states, n_states)
        edge_logK[model$active_pairs] <-
          logP[model$active_pairs] + moments$log_mean
        edge_q[model$active_pairs] <- moments$q_mean
        edge_q2[model$active_pairs] <- moments$q2_mean
        logK[[j]] <- edge_logK
        q_mean[[j]] <- edge_q
        q2_mean[[j]] <- edge_q2
        joint <- sweep(edge_logK, 1L, log_alpha[j - 1L, ], "+")
        next_alpha <- apply(joint, 2L, log_sum_exp)
        increments[j] <- log_sum_exp(next_alpha)
        log_alpha[j, ] <- next_alpha - increments[j]
      }
    }

    log_beta <- matrix(0, nrow = n_obs, ncol = n_states)
    if (n_obs >= 2L) {
      for (j in n_obs:2L) {
        backward_joint <- sweep(logK[[j]], 2L, log_beta[j, ], "+")
        log_beta[j - 1L, ] <- apply(backward_joint, 1L, log_sum_exp) -
          increments[j]
      }
    }

    expected_q <- expected_q2 <- numeric(n_obs)
    initial_log_weight <- logK[[1L]] + log_beta[1L, ] - increments[1L]
    initial_weight <- exp(initial_log_weight - log_sum_exp(initial_log_weight))
    expected_q[1L] <- sum(initial_weight * q_mean[[1L]], na.rm = TRUE)
    expected_q2[1L] <- sum(initial_weight * q2_mean[[1L]], na.rm = TRUE)
    if (n_obs >= 2L) {
      for (j in 2:n_obs) {
        edge_log_weight <- sweep(logK[[j]], 1L, log_alpha[j - 1L, ], "+")
        edge_log_weight <- sweep(edge_log_weight, 2L, log_beta[j, ], "+") -
          increments[j]
        normalizer <- log_sum_exp(as.numeric(edge_log_weight))
        edge_weight <- exp(edge_log_weight - normalizer)
        expected_q[j] <- sum(edge_weight * q_mean[[j]], na.rm = TRUE)
        expected_q2[j] <- sum(edge_weight * q2_mean[[j]], na.rm = TRUE)
      }
    }
    list(
      loglik = sum(increments),
      expected_q = expected_q,
      expected_q2 = expected_q2,
      log_alpha = if (retain) log_alpha else NULL,
      increments = if (retain) increments else NULL
    )
  }

  converged <- FALSE
  relative_change <- Inf
  loglik_path <- numeric(max_iter)
  for (iter in seq_len(max_iter)) {
    estep <- e_step(gamma, dispersion, retain = FALSE)
    loglik_path[iter] <- estep$loglik
    t_obs <- vapply(terms_list, `[[`, numeric(1L), "t")
    working_var <- vapply(terms_list, `[[`, numeric(1L), "var")
    gamma_new <- sum(t_obs * estep$expected_q / working_var) /
      sum(estep$expected_q2 / working_var)
    gamma_new <- min(max(gamma_new, beta_bounds[1L] ^ 2),
                     beta_bounds[2L] ^ 2)
    residual2 <- t_obs ^ 2 - 2 * gamma_new * t_obs * estep$expected_q +
      gamma_new ^ 2 * estep$expected_q2
    dispersion_new <- mean(residual2 / working_var)
    dispersion_new <- min(max(dispersion_new, dispersion_bounds[1L]),
                          dispersion_bounds[2L])
    relative_change <- max(abs(gamma_new - gamma) / max(gamma, 1e-12),
                           abs(dispersion_new - dispersion) /
                             max(dispersion, 1e-12))
    gamma <- gamma_new
    dispersion <- dispersion_new
    if (relative_change < tol) {
      converged <- TRUE
      break
    }
  }
  final <- e_step(gamma, dispersion, retain = TRUE)
  filter_prob <- exp(final$log_alpha)
  entropy <- -rowSums(filter_prob * log(pmax(filter_prob,
                                             .Machine$double.xmin)))
  list(
    beta_hat = sqrt(gamma),
    emission_dispersion_hat = dispersion,
    loglik = final$loglik,
    convergence = if (converged) 0L else 1L,
    n_iter = iter,
    relative_change = relative_change,
    loglik_path = loglik_path[seq_len(iter)],
    filter_diagnostics = list(
      filter_prob = filter_prob,
      loglik_increment = final$increments,
      mean_filter_max_prob = mean(apply(filter_prob, 1L, max)),
      mean_filter_entropy = mean(entropy) / log(n_states),
      mean_filter_boundary_prob = mean(filter_prob[, 1L] +
                                         filter_prob[, n_states])
    )
  )
}

refine_rho_profile_quadratic <- function(grid) {
  profile_rows <- do.call(rbind, lapply(split(grid, grid$rho), function(x) {
    x[which.max(x$loglik), , drop = FALSE]
  }))
  profile_rows <- profile_rows[order(profile_rows$rho), , drop = FALSE]
  best_i <- which.max(profile_rows$loglik)
  base <- profile_rows[best_i, , drop = FALSE]
  out <- list(
    rho = base$rho,
    beta = base$beta,
    sigma1 = base$beta / base$rho,
    dispersion = base$emission_dispersion,
    loglik = base$loglik,
    refined = FALSE,
    curvature = NA_real_
  )
  if (best_i <= 1L || best_i >= nrow(profile_rows)) return(out)
  local <- profile_rows[(best_i - 1L):(best_i + 1L), , drop = FALSE]
  x <- log(local$rho)
  fit <- stats::lm(local$loglik ~ x + I(x ^ 2))
  coef <- stats::coef(fit)
  if (length(coef) < 3L || any(!is.finite(coef)) || coef[3L] >= 0) return(out)
  x_star <- -coef[2L] / (2 * coef[3L])
  if (!is.finite(x_star) || x_star <= min(x) || x_star >= max(x)) return(out)
  rho_star <- exp(x_star)
  beta_star <- stats::approx(x, local$beta, xout = x_star,
                             method = "linear", rule = 2)$y
  dispersion_star <- stats::approx(x, local$emission_dispersion,
                                   xout = x_star, method = "linear",
                                   rule = 2)$y
  list(
    rho = rho_star,
    beta = beta_star,
    sigma1 = beta_star / rho_star,
    dispersion = dispersion_star,
    loglik = unname(coef[1L] + coef[2L] * x_star + coef[3L] * x_star ^ 2),
    refined = TRUE,
    curvature = -2 * unname(coef[3L])
  )
}

evaluate_endpoint_grid <- function(y_mat, obs_cov, freqs, start_state,
                                   end_state, q_bar, models, sigma1_grid,
                                   rho_grid, likelihood_type,
                                   emission_mode = "full_gaussian",
                                   emission_dispersion = 1,
                                   profile_beta = FALSE,
                                   beta_bounds = c(0.20, 5.00),
                                   profile_dispersion = FALSE,
                                   profile_engine = "optim",
                                   dispersion_bounds = c(0.25, Inf),
                                   dispersion_tol = 0.02) {
  n_dispersion <- if (profile_dispersion) 1L else length(emission_dispersion)
  n_sigma <- if (profile_beta) 1L else length(sigma1_grid)
  rows <- vector("list", n_sigma * length(rho_grid) * n_dispersion)
  k <- 0L
  terms_list <- lapply(seq_len(nrow(y_mat)), function(j) {
    emission_terms(y_mat[j, ], obs_cov[j, , ], freqs,
                   emission_mode = emission_mode)
  })
  projected_t <- vapply(terms_list, function(x) {
    if (!is.null(x$t) && is.finite(x$t)) x$t else NA_real_
  }, numeric(1L))
  q_reference_median <- stationary_q_quantile(0.5)
  beta_initial <- sqrt(max(stats::median(projected_t, na.rm = TRUE) /
                             q_reference_median,
                           beta_bounds[1L] ^ 2))
  beta_initial <- min(max(beta_initial, beta_bounds[1L]), beta_bounds[2L])
  em_fits <- vector("list", length(rho_grid))
  eval_loglik <- function(sigma1, rho, model, dispersion) {
    switch(
      likelihood_type,
      oracle_qbar_known_endpoints = known_endpoint_oracle_qbar_loglik(
        y_mat, obs_cov, freqs, start_state, end_state, q_bar, model,
        sigma1, rho, emission_mode = emission_mode,
        emission_dispersion = dispersion, terms_list = terms_list
      ),
      integrated_qbar_known_endpoints = known_endpoint_integrated_qbar_loglik(
        y_mat, obs_cov, freqs, start_state, end_state, model, sigma1, rho,
        emission_mode = emission_mode,
        emission_dispersion = dispersion, terms_list = terms_list
      ),
      latent_endpoint_hmm = latent_endpoint_hmm_loglik(
        y_mat, obs_cov, freqs, model, sigma1, rho,
        emission_mode = emission_mode,
        emission_dispersion = dispersion, terms_list = terms_list
      ),
      stop("Unknown likelihood_type.")
    )
  }
  for (i in seq_along(rho_grid)) {
    rho <- rho_grid[i]
    model <- models[[i]]
    if (profile_beta) {
      optimizer_iterations <- NA_integer_
      optimizer_relative_change <- NA_real_
      if (profile_dispersion && identical(profile_engine, "em") &&
          identical(likelihood_type, "latent_endpoint_hmm") &&
          emission_mode %in% c("projected_quasi", "projected_diag")) {
        em <- fit_latent_emission_em(
          terms_list = terms_list,
          model = model,
          beta_initial = beta_initial,
          dispersion_initial = max(1, emission_dispersion[1L]),
          beta_bounds = beta_bounds,
          dispersion_bounds = dispersion_bounds
        )
        beta_hat <- em$beta_hat
        dispersion_hat <- em$emission_dispersion_hat
        loglik_hat <- em$loglik
        convergence <- em$convergence
        optimizer_iterations <- em$n_iter
        optimizer_relative_change <- em$relative_change
        em_fits[[i]] <- em
      } else if (profile_dispersion) {
        opt <- stats::optim(
          par = log(c(beta_initial, max(1, emission_dispersion[1L]))),
          fn = function(par) {
            beta <- exp(par[1L])
            dispersion <- exp(par[2L])
            ll <- eval_loglik(beta / rho, rho, model, dispersion)
            if (is.finite(ll)) -ll else .Machine$double.xmax / 100
          },
          method = "L-BFGS-B",
          lower = log(c(beta_bounds[1L], dispersion_bounds[1L])),
          upper = log(c(beta_bounds[2L], dispersion_bounds[2L])),
          control = list(factr = 1e8, pgtol = 1e-6)
        )
        beta_hat <- exp(opt$par[1L])
        dispersion_hat <- exp(opt$par[2L])
        loglik_hat <- -opt$value
        convergence <- opt$convergence
      } else {
        opt <- stats::optimize(
          function(log_beta) {
            ll <- eval_loglik(exp(log_beta) / rho, rho, model,
                              emission_dispersion[1L])
            if (is.finite(ll)) -ll else .Machine$double.xmax / 100
          },
          interval = log(beta_bounds), tol = dispersion_tol
        )
        beta_hat <- exp(opt$minimum)
        dispersion_hat <- emission_dispersion[1L]
        loglik_hat <- -opt$objective
        convergence <- 0L
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        sigma1 = beta_hat / rho,
        rho = rho,
        beta = beta_hat,
        emission_dispersion = dispersion_hat,
        optimizer_convergence = convergence,
        optimizer_iterations = optimizer_iterations,
        optimizer_relative_change = optimizer_relative_change,
        loglik = loglik_hat
      )
      next
    }
    for (sigma1 in sigma1_grid) {
      if (profile_dispersion) {
        opt <- stats::optimize(
          function(log_dispersion) {
            ll <- eval_loglik(sigma1, rho, model, exp(log_dispersion))
            if (is.finite(ll)) -ll else .Machine$double.xmax / 100
          },
          interval = log(dispersion_bounds),
          tol = dispersion_tol
        )
        dispersion_values <- exp(opt$minimum)
        loglik_values <- -opt$objective
      } else {
        dispersion_values <- emission_dispersion
        loglik_values <- vapply(dispersion_values, function(dispersion) {
          eval_loglik(sigma1, rho, model, dispersion)
        }, numeric(1L))
      }
      for (d_i in seq_along(dispersion_values)) {
        k <- k + 1L
        rows[[k]] <- data.frame(
          sigma1 = sigma1,
          rho = rho,
          beta = sigma1 * rho,
          emission_dispersion = dispersion_values[d_i],
          optimizer_convergence = 0L,
          optimizer_iterations = NA_integer_,
          optimizer_relative_change = NA_real_,
          loglik = loglik_values[d_i]
        )
      }
    }
  }
  grid <- do.call(rbind, rows)
  grid <- grid[is.finite(grid$loglik), , drop = FALSE]
  if (!nrow(grid)) stop("No finite endpoint likelihood evaluations.")
  best <- grid[which.max(grid$loglik), , drop = FALSE]
  refined <- if (profile_beta) refine_rho_profile_quadratic(grid) else list(
    rho = best$rho,
    beta = best$beta,
    sigma1 = best$sigma1,
    dispersion = best$emission_dispersion,
    loglik = best$loglik,
    refined = FALSE,
    curvature = NA_real_
  )
  filter_diagnostics <- NULL
  if (identical(likelihood_type, "latent_endpoint_hmm")) {
    model_index <- which.min(abs(rho_grid - best$rho))
    if (profile_beta && profile_dispersion && identical(profile_engine, "em") &&
        !is.null(em_fits[[model_index]])) {
      filter_diagnostics <- em_fits[[model_index]]$filter_diagnostics
    } else {
      filter_diagnostics <- latent_endpoint_hmm_loglik(
        y_mat, obs_cov, freqs, models[[model_index]], best$sigma1, best$rho,
        emission_mode = emission_mode,
        emission_dispersion = best$emission_dispersion,
        terms_list = terms_list,
        return_filter = TRUE
      )
    }
  }
  list(
    sigma1_hat = refined$sigma1,
    rho_hat = refined$rho,
    beta_hat = refined$beta,
    emission_dispersion_hat = refined$dispersion,
    sigma1_hat_grid = best$sigma1,
    rho_hat_grid = best$rho,
    beta_hat_grid = best$beta,
    profile_quadratic_refined = refined$refined,
    profile_quadratic_curvature = refined$curvature,
    optimizer_convergence = best$optimizer_convergence,
    optimizer_iterations = best$optimizer_iterations,
    optimizer_relative_change = best$optimizer_relative_change,
    loglik = best$loglik,
    grid = grid,
    filter_diagnostics = filter_diagnostics,
    likelihood_type = likelihood_type
  )
}

evaluate_endpoint_transition_grid <- function(start_state, end_state, models,
                                              rho_grid) {
  rows <- vector("list", length(rho_grid))
  for (i in seq_along(rho_grid)) {
    rows[[i]] <- data.frame(
      sigma1 = NA_real_,
      rho = rho_grid[i],
      beta = NA_real_,
      emission_dispersion = NA_real_,
      loglik = known_endpoint_transition_loglik(start_state, end_state,
                                                models[[i]])
    )
  }
  grid <- do.call(rbind, rows)
  grid <- grid[is.finite(grid$loglik), , drop = FALSE]
  if (!nrow(grid)) stop("No finite endpoint transition likelihood values.")
  best <- grid[which.max(grid$loglik), , drop = FALSE]
  list(
    sigma1_hat = NA_real_,
    rho_hat = best$rho,
    beta_hat = NA_real_,
    emission_dispersion_hat = NA_real_,
    loglik = best$loglik,
    grid = grid,
    likelihood_type = "known_endpoint_transition_only"
  )
}

profile_diag <- function(grid, sigma1_true, rho_true) {
  max_ll <- max(grid$loglik, na.rm = TRUE)
  if (all(is.na(grid$sigma1))) {
    dist_true <- (log(grid$rho / rho_true)) ^ 2
    sig_width <- NA_real_
  } else {
    dist_true <- (log(grid$sigma1 / sigma1_true)) ^ 2 +
      (log(grid$rho / rho_true)) ^ 2
    prof_sig <- aggregate(loglik ~ sigma1, data = grid, FUN = max)
    prof_sig$deviance <- -2 * (prof_sig$loglik - max(prof_sig$loglik))
    sig_width <- diff(range(prof_sig$sigma1[
      prof_sig$deviance <= stats::qchisq(0.95, 1)
    ], finite = TRUE))
  }
  near_candidates <- which(dist_true == min(dist_true, na.rm = TRUE))
  near <- near_candidates[which.max(grid$loglik[near_candidates])]
  prof_rho <- aggregate(loglik ~ rho, data = grid, FUN = max)
  prof_rho$deviance <- -2 * (prof_rho$loglik - max(prof_rho$loglik))
  prof_rho <- prof_rho[order(prof_rho$rho), , drop = FALSE]
  best_rho_i <- which.max(prof_rho$loglik)
  local_i <- seq.int(max(1L, best_rho_i - 1L),
                     min(nrow(prof_rho), best_rho_i + 1L))
  rho_log_curvature <- NA_real_
  if (length(local_i) == 3L) {
    x <- log(prof_rho$rho[local_i])
    fit <- stats::lm(prof_rho$loglik[local_i] ~ x + I(x ^ 2))
    rho_log_curvature <- -2 * unname(stats::coef(fit)[3L])
  }
  sorted_ll <- sort(prof_rho$loglik, decreasing = TRUE)
  data.frame(
    likelihood_gap_true = max_ll - grid$loglik[near],
    sigma1_profile_width = sig_width,
    rho_profile_width = diff(range(prof_rho$rho[
      prof_rho$deviance <= stats::qchisq(0.95, 1)
    ], finite = TRUE)),
    rho_profile_gap = if (length(sorted_ll) >= 2L) sorted_ll[1L] - sorted_ll[2L]
      else NA_real_,
    rho_log_curvature = rho_log_curvature,
    rho_at_grid_boundary = best_rho_i %in% c(1L, nrow(prof_rho))
  )
}

observed_profile_diagnostics <- function(grid, confidence_level = 0.95) {
  profile_rows <- do.call(rbind, lapply(split(grid, grid$rho), function(x) {
    x[which.max(x$loglik), , drop = FALSE]
  }))
  profile_rows <- profile_rows[order(profile_rows$rho), , drop = FALSE]
  dense_rho <- seq(min(profile_rows$rho), max(profile_rows$rho),
                   length.out = 4001L)
  dense_loglik <- stats::approx(profile_rows$rho, profile_rows$loglik,
                                xout = dense_rho, rule = 2)$y
  dense_sigma1 <- stats::approx(profile_rows$rho, profile_rows$sigma1,
                                xout = dense_rho, rule = 2)$y
  dense_beta <- stats::approx(profile_rows$rho, profile_rows$beta,
                              xout = dense_rho, rule = 2)$y
  max_ll <- max(dense_loglik)
  dense_deviance <- -2 * (dense_loglik - max_ll)
  cutoff <- stats::qchisq(confidence_level, df = 1L)
  supported <- dense_deviance <= cutoff
  best_i <- which.max(profile_rows$loglik)
  sorted_ll <- sort(profile_rows$loglik, decreasing = TRUE)
  data.frame(
    rho_support_lower = min(dense_rho[supported]),
    rho_support_upper = max(dense_rho[supported]),
    rho_profile_width = diff(range(dense_rho[supported])),
    sigma1_support_lower = min(dense_sigma1[supported]),
    sigma1_support_upper = max(dense_sigma1[supported]),
    beta_support_lower = min(dense_beta[supported]),
    beta_support_upper = max(dense_beta[supported]),
    rho_profile_gap = if (length(sorted_ll) >= 2L) sorted_ll[1L] - sorted_ll[2L]
      else NA_real_,
    rho_at_grid_boundary = best_i %in% c(1L, nrow(profile_rows)),
    stringsAsFactors = FALSE
  )
}

estimate_markov_additive_endpoint_hmm <- function(
    R, dt, alpha_hat, sigma2_hat,
    block_horizon = 1 / 20,
    rho_grid = seq(0.7, 2.3, by = 0.2),
    beta_bounds = c(0.20, 5.00),
    dispersion_bounds = c(0.25, Inf),
    freq_mults = c(0.25, 0.50, 0.75),
    n_states = 21L,
    transition_sim_blocks = 210000L,
    substeps_per_block = 32L,
    max_qbar_samples_per_pair = 7L,
    kernel_seed = 5000L,
    kernel_method = "monte_carlo",
    endpoint_state = "q",
    confidence_level = 0.95,
    return_observations = FALSE,
    verbose = TRUE) {
  n <- min(length(R), length(dt))
  if (n < 1000L) stop("At least 1000 observed increments are required.")
  R <- as.numeric(R[seq_len(n)])
  dt <- as.numeric(dt[seq_len(n)])
  if (!is.finite(alpha_hat) || alpha_hat <= 0 || alpha_hat >= 2) {
    stop("alpha_hat must lie strictly between zero and two.")
  }
  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
    stop("sigma2_hat must be positive.")
  }
  dt_reference <- stats::median(dt[is.finite(dt) & dt > 0])
  if (!is.finite(dt_reference) || dt_reference <= 0) {
    stop("dt must contain positive finite sampling intervals.")
  }
  block_size <- max(10L, as.integer(round(block_horizon / dt_reference)))
  actual_block_h <- block_size * dt_reference
  if (is.null(freq_mults)) freq_mults <- c(0.25, 0.50, 0.75)
  rho_grid <- sort(unique(as.numeric(rho_grid)))
  rho_grid <- rho_grid[is.finite(rho_grid) & rho_grid > 0]
  if (length(rho_grid) < 3L) stop("rho_grid must contain at least three values.")

  endpoint_state <- match.arg(endpoint_state, c("v", "q"))
  state_info_local <- if (identical(endpoint_state, "q")) {
    make_q_endpoint_breaks(as.integer(n_states))
  } else {
    make_endpoint_breaks(as.integer(n_states))
  }
  if (verbose) {
    cat(sprintf("Building %d endpoint/additive kernels at H=%.6g...\n",
                length(rho_grid), actual_block_h))
  }
  kernel_method <- match.arg(kernel_method,
                             c("monte_carlo", "feynman_kac", "hybrid"))
  models_local <- lapply(rho_grid, function(rho) {
    if (identical(kernel_method, "hybrid")) {
      build_hybrid_endpoint_additive_model(
        rho = rho,
        block_h = actual_block_h,
        state_info = state_info_local,
        n_sim_blocks = as.integer(transition_sim_blocks),
        substeps_per_block = as.integer(substeps_per_block),
        q_nodes_per_pair = as.integer(max_qbar_samples_per_pair),
        seed = as.integer(kernel_seed)
      )
    } else if (identical(kernel_method, "feynman_kac")) {
      build_feynman_kac_endpoint_model(
        rho = rho,
        block_h = actual_block_h,
        state_info = state_info_local,
        substeps_per_block = as.integer(substeps_per_block),
        q_nodes_per_pair = as.integer(max_qbar_samples_per_pair)
      )
    } else {
      if (identical(endpoint_state, "q")) {
        simulate_q_endpoint_additive_model(
          rho = rho,
          block_h = actual_block_h,
          state_info = state_info_local,
          n_sim_blocks = as.integer(transition_sim_blocks),
          substeps_per_block = as.integer(substeps_per_block),
          max_qbar_samples_per_pair = as.integer(max_qbar_samples_per_pair),
          state_model = "conditional",
          seed = as.integer(kernel_seed)
        )
      } else {
        simulate_endpoint_additive_model(
          rho = rho,
          block_h = actual_block_h,
          state_info = state_info_local,
          n_sim_blocks = as.integer(transition_sim_blocks),
          substeps_per_block = as.integer(substeps_per_block),
          max_qbar_samples_per_pair = as.integer(max_qbar_samples_per_pair),
          state_model = "conditional",
          seed = as.integer(kernel_seed)
        )
      }
    }
  })
  observations <- block_cf_log_observations(
    R = R,
    dt = dt,
    alpha = alpha_hat,
    sigma2 = sigma2_hat,
    block_size = block_size,
    freq_mults = freq_mults
  )
  estimate <- evaluate_endpoint_grid(
    y_mat = observations$y_dejumped_bias_corrected,
    obs_cov = observations$obs_cov,
    freqs = observations$freqs,
    start_state = NULL,
    end_state = NULL,
    q_bar = NULL,
    models = models_local,
    sigma1_grid = 1,
    rho_grid = rho_grid,
    likelihood_type = "latent_endpoint_hmm",
    emission_mode = "projected_diag",
    emission_dispersion = 1,
    profile_beta = TRUE,
    beta_bounds = beta_bounds,
    profile_dispersion = TRUE,
    profile_engine = "em",
    dispersion_bounds = dispersion_bounds
  )
  profile <- observed_profile_diagnostics(estimate$grid,
                                          confidence_level = confidence_level)
  warnings <- character(0)
  if (isTRUE(profile$rho_at_grid_boundary)) {
    warnings <- c(warnings, "rho_profile_maximum_at_grid_boundary")
  }
  if (!identical(as.integer(estimate$optimizer_convergence), 0L)) {
    warnings <- c(warnings, "em_profile_not_converged")
  }
  if (estimate$beta_hat <= beta_bounds[1L] * 1.001 ||
      estimate$beta_hat >= beta_bounds[2L] / 1.001) {
    warnings <- c(warnings, "beta_at_parameter_boundary")
  }
  if (estimate$emission_dispersion_hat <= dispersion_bounds[1L] * 1.001 ||
      (is.finite(dispersion_bounds[2L]) &&
       estimate$emission_dispersion_hat >= dispersion_bounds[2L] / 1.001)) {
    warnings <- c(warnings, "dispersion_at_parameter_boundary")
  }
  output <- list(
    alpha_hat = alpha_hat,
    sigma2_hat = sigma2_hat,
    sigma1_hat = estimate$sigma1_hat,
    rho_hat = estimate$rho_hat,
    beta_hat = estimate$beta_hat,
    kappa_hat = estimate$emission_dispersion_hat,
    block_size = block_size,
    block_horizon = actual_block_h,
    n_blocks = observations$n_blocks,
    freq_mults = freq_mults,
    rho_grid = rho_grid,
    profile = estimate$grid,
    profile_diagnostics = profile,
    filter_diagnostics = estimate$filter_diagnostics,
    optimizer_convergence = estimate$optimizer_convergence,
    optimizer_iterations = estimate$optimizer_iterations,
    optimizer_relative_change = estimate$optimizer_relative_change,
    profile_quadratic_refined = estimate$profile_quadratic_refined,
    profile_quadratic_curvature = estimate$profile_quadratic_curvature,
    warnings = unique(warnings),
    method = "markov_additive_endpoint_ecf_hmm"
  )
  output$kernel_method <- kernel_method
  output$endpoint_state <- endpoint_state
  if (return_observations) output$observations <- observations
  if (verbose) {
    cat(sprintf(
      "Endpoint HMM: sigma1=%.6f rho=%.6f beta=%.6f kappa=%.4f blocks=%d\n",
      output$sigma1_hat, output$rho_hat, output$beta_hat,
      output$kappa_hat, output$n_blocks
    ))
    if (length(output$warnings)) {
      cat("Warnings:", paste(output$warnings, collapse = ", "), "\n")
    }
  }
  output
}

if (!identical(Sys.getenv("MARKOV_ENDPOINT_LIBRARY_ONLY", unset = "0"), "1")) {
state_info <- if (identical(endpoint_state, "q")) {
  make_q_endpoint_breaks(n_states)
} else {
  make_endpoint_breaks(n_states)
}
block_h <- cases$terminal[1L] / cases$n_steps[1L] * cases$block_size[1L]
cat("Precomputing endpoint Markov-additive transition models...\n")
models <- vector("list", length(rho_grid))
for (i in seq_along(rho_grid)) {
  models[[i]] <- if (identical(kernel_method, "hybrid")) {
    build_hybrid_endpoint_additive_model(
      rho = rho_grid[i],
      block_h = block_h,
      state_info = state_info,
      n_sim_blocks = transition_sim_blocks,
      substeps_per_block = substeps_per_block,
      q_nodes_per_pair = max_qbar_samples_per_pair,
      seed = kernel_seed
    )
  } else if (identical(kernel_method, "feynman_kac")) {
    build_feynman_kac_endpoint_model(
      rho = rho_grid[i],
      block_h = block_h,
      state_info = state_info,
      substeps_per_block = substeps_per_block,
      q_nodes_per_pair = max_qbar_samples_per_pair
    )
  } else {
    if (identical(endpoint_state, "q")) {
      simulate_q_endpoint_additive_model(
        rho = rho_grid[i], block_h = block_h, state_info = state_info,
        n_sim_blocks = transition_sim_blocks,
        substeps_per_block = substeps_per_block,
        max_qbar_samples_per_pair = max_qbar_samples_per_pair,
        state_model = state_model, seed = kernel_seed
      )
    } else {
      simulate_endpoint_additive_model(
        rho = rho_grid[i], block_h = block_h, state_info = state_info,
        n_sim_blocks = transition_sim_blocks,
        substeps_per_block = substeps_per_block,
        max_qbar_samples_per_pair = max_qbar_samples_per_pair,
        state_model = state_model, seed = kernel_seed
      )
    }
  }
}

result_rows <- list()
curve_rows <- list()
result_i <- 0L
curve_i <- 0L

for (case_i in seq_len(nrow(cases))) {
  cfg <- cases[case_i, ]
  case_tag <- sprintf("seed%d_n%d_T%s_b%d", cfg$seed, cfg$n_steps,
                      gsub("\\.", "p", as.character(cfg$terminal)),
                      cfg$block_size)
  cat(sprintf("CASE %02d/%02d %s\n", case_i, nrow(cases), case_tag))
  sim <- simulate_residual_euler(
    n_steps = cfg$n_steps,
    terminal = cfg$terminal,
    alpha = fixed$alpha,
    sigma2 = fixed$sigma2,
    sigma1 = fixed$sigma1,
    rho = fixed$rho,
    seed = cfg$seed,
    keep_v = !observed_only
  )
  tail_est <- if (identical(tail_input_mode, "oracle")) {
    list(alpha_hat = NA_real_, sigma2_hat = NA_real_, k_hat = NA_integer_,
         hill_score = NA_real_)
  } else {
    estimate_alpha_sigma2_hill(sim$R, sim$dt)
  }
  alpha_input <- switch(
    tail_input_mode,
    oracle = fixed$alpha,
    estimated = tail_est$alpha_hat,
    `true-alpha` = fixed$alpha,
    `true-sigma2` = tail_est$alpha_hat
  )
  sigma2_input <- switch(
    tail_input_mode,
    oracle = fixed$sigma2,
    estimated = tail_est$sigma2_hat,
    `true-alpha` = tail_est$sigma2_hat,
    `true-sigma2` = fixed$sigma2
  )
  cat(sprintf(
    "  tail=%s alpha_used=%.5f sigma2_used=%.5f k=%s\n",
    tail_input_mode, alpha_input, sigma2_input,
    if (is.finite(tail_est$k_hat)) as.character(tail_est$k_hat) else "oracle"
  ))
  ecf_obs <- block_cf_log_observations(
    R = sim$R,
    dt = sim$dt,
    alpha = alpha_input,
    sigma2 = sigma2_input,
    block_size = cfg$block_size,
    freq_mults = if (is.null(freq_mults)) {
      ecf_frequency_multipliers(cfg$block_size)
    } else freq_mults
  )
  truth <- if (observed_only) {
    list(q_bar = NULL, start_state = NULL, end_state = NULL)
  } else {
    block_truth(sim$V, cfg$block_size, ecf_obs$n_blocks, state_info,
                endpoint_state = endpoint_state)
  }
  model_y <- if (observed_only) NULL else simulate_y_from_qbar(
    q_bar = truth$q_bar,
    ecf_obs = ecf_obs,
    sigma1 = fixed$sigma1,
    rho = fixed$rho,
    seed = cfg$seed + 733L
  )
  y_sources <- list(
    model_y = model_y,
    empirical_y = ecf_obs$y_dejumped,
    empirical_y_bias_corrected = ecf_obs$y_dejumped_bias_corrected
  )
  y_source_names <- if (corrected) {
    c("model_y", "empirical_y_bias_corrected")
  } else {
    c("model_y", "empirical_y")
  }
  if (any(args == "--model-y-only")) y_source_names <- "model_y"
  if (any(args == "--empirical-y-only")) {
    y_source_names <- if (corrected) {
      "empirical_y_bias_corrected"
    } else {
      "empirical_y"
    }
  }
  if (observed_only) {
    y_source_names <- if (corrected) "empirical_y_bias_corrected" else "empirical_y"
  }
  emission_modes <- if (corrected) {
    if (identical(covariance_structure, "full")) {
      "projected_quasi"
    } else {
      "projected_diag"
    }
  } else {
    "full_gaussian"
  }

  if (!observed_only) {
  trans_est <- evaluate_endpoint_transition_grid(
    start_state = truth$start_state,
    end_state = truth$end_state,
    models = models,
    rho_grid = rho_grid
  )
  trans_diag <- profile_diag(trans_est$grid, fixed$sigma1, fixed$rho)
  result_i <- result_i + 1L
  result_rows[[result_i]] <- cbind(
    data.frame(
      case_tag = case_tag,
      seed = cfg$seed,
      n_steps = cfg$n_steps,
      terminal = cfg$terminal,
      block_size = cfg$block_size,
      block_h = block_h,
      n_blocks = ecf_obs$n_blocks,
      method = "known_endpoint_transition_only",
      y_source = "none",
      likelihood_type = "known_endpoint_transition_only",
      emission_mode = "none",
      alpha = fixed$alpha,
      sigma2 = fixed$sigma2,
      alpha_hat = tail_est$alpha_hat,
      sigma2_hat = tail_est$sigma2_hat,
      alpha_used = alpha_input,
      sigma2_used = sigma2_input,
      alpha_abs_err = if (is.finite(tail_est$alpha_hat))
        abs(tail_est$alpha_hat - fixed$alpha) else NA_real_,
      sigma2_abs_err = if (is.finite(tail_est$sigma2_hat))
        abs(tail_est$sigma2_hat - fixed$sigma2) else NA_real_,
      k_hat = tail_est$k_hat,
      sigma1_true = fixed$sigma1,
      rho_true = fixed$rho,
      beta_true = beta_true,
      sigma1_hat = NA_real_,
      rho_hat = trans_est$rho_hat,
      beta_hat = NA_real_,
      emission_dispersion_hat = NA_real_,
      sigma1_hat_grid = NA_real_,
      rho_hat_grid = trans_est$rho_hat,
      beta_hat_grid = NA_real_,
      profile_quadratic_refined = FALSE,
      profile_quadratic_curvature = NA_real_,
      optimizer_convergence = NA_integer_,
      optimizer_iterations = NA_integer_,
      optimizer_relative_change = NA_real_,
      beta_at_boundary = NA,
      dispersion_at_boundary = NA,
      mean_filter_max_prob = NA_real_,
      mean_filter_entropy = NA_real_,
      mean_filter_boundary_prob = NA_real_,
      sigma1_abs_err = NA_real_,
      rho_abs_err = abs(trans_est$rho_hat - fixed$rho),
      beta_abs_err = NA_real_,
      downstream_both_under_010 = NA,
      all_four_under_010 = NA,
      sigma1_rel_err = NA_real_,
      rho_rel_err = abs(trans_est$rho_hat - fixed$rho) / fixed$rho,
      beta_rel_err = NA_real_
    ),
    trans_diag
  )
  curve_i <- curve_i + 1L
  curve_rows[[curve_i]] <- cbind(
    data.frame(case_tag = case_tag,
               method = "known_endpoint_transition_only"),
    trans_est$grid
  )
  }

  likelihood_types <- c("oracle_qbar_known_endpoints",
                        "integrated_qbar_known_endpoints",
                        "latent_endpoint_hmm")
  if (observed_only) likelihood_types <- "latent_endpoint_hmm"
  if (any(args == "--latent-only")) likelihood_types <- "latent_endpoint_hmm"
  if (any(args == "--known-endpoints-only")) {
    likelihood_types <- c("oracle_qbar_known_endpoints",
                          "integrated_qbar_known_endpoints")
  }
  method_specs <- expand.grid(
    y_source = y_source_names,
    likelihood_type = likelihood_types,
    emission_mode = emission_modes,
    stringsAsFactors = FALSE
  )
  for (m_i in seq_len(nrow(method_specs))) {
    spec <- method_specs[m_i, ]
    est <- evaluate_endpoint_grid(
      y_mat = y_sources[[spec$y_source]],
      obs_cov = ecf_obs$obs_cov,
      freqs = ecf_obs$freqs,
      start_state = truth$start_state,
      end_state = truth$end_state,
      q_bar = truth$q_bar,
      models = models,
      sigma1_grid = sigma1_grid,
      rho_grid = rho_grid,
      likelihood_type = spec$likelihood_type,
      emission_mode = spec$emission_mode,
      emission_dispersion = dispersion_grid,
      profile_beta = continuous_beta,
      beta_bounds = beta_bounds,
      profile_dispersion = continuous_dispersion,
      profile_engine = profile_engine,
      dispersion_bounds = dispersion_bounds,
      dispersion_tol = dispersion_tol
    )
    diag <- profile_diag(est$grid, fixed$sigma1, fixed$rho)
    method_name <- paste(spec$likelihood_type, spec$y_source,
                         spec$emission_mode, sep = "_")
    result_i <- result_i + 1L
    result_rows[[result_i]] <- cbind(
      data.frame(
        case_tag = case_tag,
        seed = cfg$seed,
        n_steps = cfg$n_steps,
        terminal = cfg$terminal,
        block_size = cfg$block_size,
        block_h = block_h,
        n_blocks = ecf_obs$n_blocks,
        method = method_name,
        y_source = spec$y_source,
        likelihood_type = spec$likelihood_type,
        emission_mode = spec$emission_mode,
        alpha = fixed$alpha,
        sigma2 = fixed$sigma2,
        alpha_hat = tail_est$alpha_hat,
        sigma2_hat = tail_est$sigma2_hat,
        alpha_used = alpha_input,
        sigma2_used = sigma2_input,
        alpha_abs_err = if (is.finite(tail_est$alpha_hat))
          abs(tail_est$alpha_hat - fixed$alpha) else NA_real_,
        sigma2_abs_err = if (is.finite(tail_est$sigma2_hat))
          abs(tail_est$sigma2_hat - fixed$sigma2) else NA_real_,
        k_hat = tail_est$k_hat,
        sigma1_true = fixed$sigma1,
        rho_true = fixed$rho,
        beta_true = beta_true,
        sigma1_hat = est$sigma1_hat,
        rho_hat = est$rho_hat,
        beta_hat = est$beta_hat,
        emission_dispersion_hat = est$emission_dispersion_hat,
        sigma1_hat_grid = est$sigma1_hat_grid,
        rho_hat_grid = est$rho_hat_grid,
        beta_hat_grid = est$beta_hat_grid,
        profile_quadratic_refined = est$profile_quadratic_refined,
        profile_quadratic_curvature = est$profile_quadratic_curvature,
        optimizer_convergence = est$optimizer_convergence,
        optimizer_iterations = est$optimizer_iterations,
        optimizer_relative_change = est$optimizer_relative_change,
        beta_at_boundary = est$beta_hat <= beta_bounds[1L] * 1.001 ||
          est$beta_hat >= beta_bounds[2L] / 1.001,
        dispersion_at_boundary =
          est$emission_dispersion_hat <= dispersion_bounds[1L] * 1.001 ||
          est$emission_dispersion_hat >= dispersion_bounds[2L] / 1.001,
        mean_filter_max_prob = if (is.null(est$filter_diagnostics)) NA_real_
          else est$filter_diagnostics$mean_filter_max_prob,
        mean_filter_entropy = if (is.null(est$filter_diagnostics)) NA_real_
          else est$filter_diagnostics$mean_filter_entropy,
        mean_filter_boundary_prob = if (is.null(est$filter_diagnostics)) NA_real_
          else est$filter_diagnostics$mean_filter_boundary_prob,
        sigma1_abs_err = abs(est$sigma1_hat - fixed$sigma1),
        rho_abs_err = abs(est$rho_hat - fixed$rho),
        beta_abs_err = abs(est$beta_hat - beta_true),
        downstream_both_under_010 =
          abs(est$sigma1_hat - fixed$sigma1) < 0.10 &&
          abs(est$rho_hat - fixed$rho) < 0.10,
        all_four_under_010 =
          is.finite(tail_est$alpha_hat) && is.finite(tail_est$sigma2_hat) &&
          abs(tail_est$alpha_hat - fixed$alpha) < 0.10 &&
          abs(tail_est$sigma2_hat - fixed$sigma2) < 0.10 &&
          abs(est$sigma1_hat - fixed$sigma1) < 0.10 &&
          abs(est$rho_hat - fixed$rho) < 0.10,
        sigma1_rel_err = abs(est$sigma1_hat - fixed$sigma1) / fixed$sigma1,
        rho_rel_err = abs(est$rho_hat - fixed$rho) / fixed$rho,
        beta_rel_err = abs(est$beta_hat - beta_true) / beta_true
      ),
      diag
    )
    cat(sprintf(
      paste0("  endpoint HMM sigma1=%.5f (err %.5f) rho=%.5f (err %.5f) ",
             "beta=%.5f kappa=%.4f refined=%s\n"),
      est$sigma1_hat, abs(est$sigma1_hat - fixed$sigma1),
      est$rho_hat, abs(est$rho_hat - fixed$rho), est$beta_hat,
      est$emission_dispersion_hat, est$profile_quadratic_refined
    ))
    curve_i <- curve_i + 1L
    curve_rows[[curve_i]] <- cbind(
      data.frame(case_tag = case_tag, method = method_name),
      est$grid
    )
  }

  activity <- block_cf_activity(
    R = sim$R,
    dt = sim$dt,
    alpha = alpha_input,
    sigma2 = sigma2_input,
    block_size = cfg$block_size,
    freq_mults = if (is.null(freq_mults)) {
      ecf_frequency_multipliers(cfg$block_size)
    } else freq_mults
  )
  q_sample <- stationary_q_sample(120000L, seed = cfg$seed + 29L)
  q_thresholds <- stats::quantile(q_sample, probs = c(0.40, 0.75),
                                  names = FALSE, type = 8)
  scalar_beta <- estimate_beta_from_activity(activity$activity, q_sample,
                                             method = "wls")
  scalar_q <- activity$activity / scalar_beta$beta_hat ^ 2
  scalar_rho <- estimate_rho_generator_gmm(
    q_hat = scalar_q,
    block_h = activity$block_h,
    tests = tests,
    instruments = instruments,
    q_thresholds = q_thresholds,
    lags = lags,
    standardize = TRUE,
    rho_grid = rho_grid
  )
  result_i <- result_i + 1L
  result_rows[[result_i]] <- data.frame(
    case_tag = case_tag,
    seed = cfg$seed,
    n_steps = cfg$n_steps,
    terminal = cfg$terminal,
    block_size = cfg$block_size,
    block_h = block_h,
    n_blocks = ecf_obs$n_blocks,
    method = "scalar_point_generator_gmm_baseline",
    y_source = "scalar_activity",
    likelihood_type = "baseline",
    emission_mode = "none",
    alpha = fixed$alpha,
    sigma2 = fixed$sigma2,
    alpha_hat = tail_est$alpha_hat,
    sigma2_hat = tail_est$sigma2_hat,
    alpha_used = alpha_input,
    sigma2_used = sigma2_input,
    alpha_abs_err = if (is.finite(tail_est$alpha_hat))
      abs(tail_est$alpha_hat - fixed$alpha) else NA_real_,
    sigma2_abs_err = if (is.finite(tail_est$sigma2_hat))
      abs(tail_est$sigma2_hat - fixed$sigma2) else NA_real_,
    k_hat = tail_est$k_hat,
    sigma1_true = fixed$sigma1,
    rho_true = fixed$rho,
    beta_true = beta_true,
    sigma1_hat = scalar_beta$beta_hat / scalar_rho$rho_hat,
    rho_hat = scalar_rho$rho_hat,
    beta_hat = scalar_beta$beta_hat,
    emission_dispersion_hat = NA_real_,
    sigma1_hat_grid = scalar_beta$beta_hat / scalar_rho$rho_hat,
    rho_hat_grid = scalar_rho$rho_hat,
    beta_hat_grid = scalar_beta$beta_hat,
    profile_quadratic_refined = FALSE,
    profile_quadratic_curvature = NA_real_,
    optimizer_convergence = NA_integer_,
    optimizer_iterations = NA_integer_,
    optimizer_relative_change = NA_real_,
    beta_at_boundary = NA,
    dispersion_at_boundary = NA,
    mean_filter_max_prob = NA_real_,
    mean_filter_entropy = NA_real_,
    mean_filter_boundary_prob = NA_real_,
    sigma1_abs_err = abs(scalar_beta$beta_hat / scalar_rho$rho_hat -
                           fixed$sigma1),
    rho_abs_err = abs(scalar_rho$rho_hat - fixed$rho),
    beta_abs_err = abs(scalar_beta$beta_hat - beta_true),
    downstream_both_under_010 =
      abs(scalar_beta$beta_hat / scalar_rho$rho_hat - fixed$sigma1) < 0.10 &&
      abs(scalar_rho$rho_hat - fixed$rho) < 0.10,
    all_four_under_010 =
      is.finite(tail_est$alpha_hat) && is.finite(tail_est$sigma2_hat) &&
      abs(tail_est$alpha_hat - fixed$alpha) < 0.10 &&
      abs(tail_est$sigma2_hat - fixed$sigma2) < 0.10 &&
      abs(scalar_beta$beta_hat / scalar_rho$rho_hat - fixed$sigma1) < 0.10 &&
      abs(scalar_rho$rho_hat - fixed$rho) < 0.10,
    sigma1_rel_err = abs(scalar_beta$beta_hat / scalar_rho$rho_hat -
                           fixed$sigma1) / fixed$sigma1,
    rho_rel_err = abs(scalar_rho$rho_hat - fixed$rho) / fixed$rho,
    beta_rel_err = abs(scalar_beta$beta_hat - beta_true) / beta_true,
    likelihood_gap_true = NA_real_,
    sigma1_profile_width = NA_real_,
    rho_profile_width = NA_real_,
    rho_profile_gap = NA_real_,
    rho_log_curvature = NA_real_,
    rho_at_grid_boundary = NA
  )
}

results <- do.call(rbind, result_rows)
curves <- bind_rows_fill(curve_rows)
results$model_label <- model_label
results$endpoint_states <- n_states
results$transition_sim_blocks <- transition_sim_blocks
results$substeps_per_block <- substeps_per_block
results$initial_v_model <- 0
results$state_model <- state_model
results$kernel_seed <- kernel_seed
results$kernel_method <- kernel_method
results$endpoint_state <- endpoint_state
results$emission_dispersion <- emission_dispersion
results$emission_dispersion_grid <- paste(dispersion_grid, collapse = ";")
results$continuous_dispersion <- continuous_dispersion
results$continuous_beta <- continuous_beta
results$profile_engine <- profile_engine
results$beta_lower <- beta_bounds[1L]
results$beta_upper <- beta_bounds[2L]
results$dispersion_lower <- dispersion_bounds[1L]
results$dispersion_upper <- dispersion_bounds[2L]
results$dispersion_tol <- dispersion_tol
results$covariance_structure <- covariance_structure
results$tail_input_mode <- tail_input_mode
results$observed_only <- observed_only
curves$model_label <- model_label
curves$endpoint_states <- n_states
curves$transition_sim_blocks <- transition_sim_blocks
curves$substeps_per_block <- substeps_per_block
curves$state_model <- state_model
curves$kernel_seed <- kernel_seed
curves$kernel_method <- kernel_method
curves$endpoint_state <- endpoint_state
curves$emission_dispersion_grid <- paste(dispersion_grid, collapse = ";")
curves$continuous_dispersion <- continuous_dispersion
curves$continuous_beta <- continuous_beta
curves$profile_engine <- profile_engine
curves$beta_lower <- beta_bounds[1L]
curves$beta_upper <- beta_bounds[2L]
curves$dispersion_lower <- dispersion_bounds[1L]
curves$dispersion_upper <- dispersion_bounds[2L]
curves$dispersion_tol <- dispersion_tol
curves$covariance_structure <- covariance_structure

prefix <- paste0(
  if (quick) "quick_" else if (focused) "focused_" else "broad_",
  if (corrected) "corrected_" else "",
  gsub("[^A-Za-z0-9_-]", "_", model_label), "_"
)
results_path <- file.path(
  out_dir, sprintf("markov_additive_endpoint_hmm_diagnostics_%s%s.csv",
                   prefix, stamp)
)
summary_path <- file.path(
  out_dir, sprintf("markov_additive_endpoint_hmm_summary_%s%s.csv",
                   prefix, stamp)
)
curves_path <- file.path(
  out_dir, sprintf("markov_additive_endpoint_hmm_curves_%s%s.csv",
                   prefix, stamp)
)
utils::write.csv(results, results_path, row.names = FALSE)
utils::write.csv(curves, curves_path, row.names = FALSE)

summary <- aggregate(
  results[, c("sigma1_abs_err", "rho_abs_err", "beta_abs_err",
              "sigma1_rel_err", "rho_rel_err", "beta_rel_err",
              "likelihood_gap_true", "sigma1_profile_width",
              "rho_profile_width")],
  by = list(method = results$method, n_steps = results$n_steps,
            terminal = results$terminal),
  FUN = function(x) stats::median(x, na.rm = TRUE)
)
utils::write.csv(summary, summary_path, row.names = FALSE)

cat("\nWrote:\n")
cat("  ", results_path, "\n", sep = "")
cat("  ", summary_path, "\n", sep = "")
cat("  ", curves_path, "\n", sep = "")
}
