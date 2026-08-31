# Exact finite-frequency observation model for the Markov-additive endpoint HMM.
#
# Conditional on a within-block activity path Q_s, the jump-deconvolved
# empirical characteristic function has mean
#
#   log { H^{-1} integral exp(-u^2 beta^2 Q_s / 2) ds }.
#
# The linear-Qbar endpoint HMM replaces this expression by
# -u^2 beta^2 Qbar / 2.  This prototype retains representative activity paths
# in the endpoint kernel and evaluates the finite-frequency mean directly.

Sys.setenv(MARKOV_ENDPOINT_LIBRARY_ONLY = "1")
source("work/diagnose_markov_additive_endpoint_hmm.R")
Sys.unsetenv("MARKOV_ENDPOINT_LIBRARY_ONLY")

select_representative_activity_paths <- function(q_bar, q_path, max_nodes) {
  keep <- is.finite(q_bar) & q_bar >= 1 &
    apply(q_path, 1L, function(x) all(is.finite(x) & x >= 1))
  q_bar <- q_bar[keep]
  q_path <- q_path[keep, , drop = FALSE]
  if (!length(q_bar)) {
    return(list(q_bar = numeric(0), q_path = matrix(numeric(0), 0L,
                                                     ncol(q_path))))
  }
  order_index <- order(q_bar)
  if (length(order_index) > max_nodes) {
    probs <- (seq_len(max_nodes) - 0.5) / max_nodes
    ranks <- pmin(length(order_index), pmax(1L,
      as.integer(ceiling(probs * length(order_index)))))
    order_index <- order_index[ranks]
  }
  list(q_bar = q_bar[order_index],
       q_path = q_path[order_index, , drop = FALSE])
}

simulate_q_endpoint_cf_model <- function(
    rho, block_h, state_info, n_sim_blocks,
    substeps_per_block = 32L, max_path_nodes_per_pair = 7L,
    initial_q = 1, state_model = c("conditional", "collocation"),
    seed = 1L) {
  set.seed(seed)
  state_model <- match.arg(state_model)
  n_states <- length(state_info$centers)
  n_substeps <- as.integer(substeps_per_block)
  dt_sub <- block_h / n_substeps
  rho2 <- rho ^ 2
  n_paths_per_state <- max(100L, as.integer(ceiling(n_sim_blocks / n_states)))
  stationary_pool <- stationary_q_sample(max(300000L, 20L * n_paths_per_state),
                                         seed = seed + 103L)
  stationary_bin <- state_index(stationary_pool, state_info$cut_breaks)
  counts <- matrix(0, n_states, n_states)
  q_samples <- vector("list", n_states * n_states)
  q_path_samples <- vector("list", n_states * n_states)

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
    q_path <- matrix(NA_real_, nrow = n, ncol = n_substeps)
    for (step in seq_len(n_substeps)) {
      q_path[, step] <- 1 + v ^ 2
      v <- v - rho2 * v * dt_sub +
        rho * sqrt(pmax(1 + v ^ 2, 1e-12)) *
        stats::rnorm(n, sd = sqrt(dt_sub))
    }
    list(q_end = 1 + v ^ 2, q_bar = rowMeans(q_path), q_path = q_path)
  }

  for (a in seq_len(n_states)) {
    paths <- simulate_paths(draw_start_q(a, n_paths_per_state))
    end_state <- state_index(paths$q_end, state_info$cut_breaks)
    for (b in seq_len(n_states)) {
      hit <- which(end_state == b)
      counts[a, b] <- length(hit)
      if (!length(hit)) next
      selected <- select_representative_activity_paths(
        paths$q_bar[hit], paths$q_path[hit, , drop = FALSE],
        max_path_nodes_per_pair
      )
      idx <- pair_index(a, b, n_states)
      q_samples[[idx]] <- selected$q_bar
      q_path_samples[[idx]] <- selected$q_path
    }
  }
  transition <- counts / rowSums(counts)

  initial_paths_raw <- simulate_paths(rep(initial_q, n_paths_per_state))
  initial_end <- state_index(initial_paths_raw$q_end, state_info$cut_breaks)
  initial_counts <- tabulate(initial_end, nbins = n_states)
  initial_transition <- initial_counts / sum(initial_counts)
  initial_selected <- lapply(seq_len(n_states), function(b) {
    hit <- which(initial_end == b)
    if (!length(hit)) {
      return(list(q_bar = numeric(0),
                  q_path = matrix(numeric(0), 0L, n_substeps)))
    }
    select_representative_activity_paths(
      initial_paths_raw$q_bar[hit],
      initial_paths_raw$q_path[hit, , drop = FALSE],
      max_path_nodes_per_pair
    )
  })
  initial_q_samples <- lapply(initial_selected, `[[`, "q_bar")
  initial_q_path_samples <- lapply(initial_selected, `[[`, "q_path")

  active_pairs <- which(counts > 0, arr.ind = TRUE)
  pair_node_count <- vapply(seq_len(nrow(active_pairs)), function(i) {
    length(q_samples[[pair_index(active_pairs[i, 1L], active_pairs[i, 2L],
                                  n_states)]])
  }, integer(1L))
  max_nodes <- max(pair_node_count)
  pair_q_nodes <- matrix(NA_real_, nrow(active_pairs), max_nodes)
  pair_q_paths <- array(NA_real_,
                        dim = c(nrow(active_pairs), max_nodes, n_substeps))
  for (i in seq_len(nrow(active_pairs))) {
    idx <- pair_index(active_pairs[i, 1L], active_pairs[i, 2L], n_states)
    n_i <- pair_node_count[i]
    pair_q_nodes[i, seq_len(n_i)] <- q_samples[[idx]]
    pair_q_paths[i, seq_len(n_i), ] <- q_path_samples[[idx]]
  }

  initial_node_count <- vapply(initial_q_samples, length, integer(1L))
  initial_max_nodes <- max(initial_node_count)
  initial_q_paths <- array(NA_real_,
                           dim = c(n_states, initial_max_nodes, n_substeps))
  for (b in seq_len(n_states)) {
    if (initial_node_count[b] > 0L) {
      initial_q_paths[b, seq_len(initial_node_count[b]), ] <-
        initial_q_path_samples[[b]]
    }
  }

  list(
    rho = rho,
    transition = transition,
    active_pairs = active_pairs,
    pair_node_count = pair_node_count,
    pair_q_nodes = pair_q_nodes,
    pair_q_paths = pair_q_paths,
    initial_transition = initial_transition,
    initial_node_count = initial_node_count,
    initial_q_paths = initial_q_paths,
    counts = counts,
    state_info = state_info,
    n_states = n_states,
    block_h = block_h,
    n_sim_blocks = n_sim_blocks,
    n_paths_per_state = n_paths_per_state,
    substeps_per_block = n_substeps,
    endpoint_state = "q",
    observation_functional = "exact_block_cf"
  )
}

exact_cf_path_means <- function(q_paths, beta, freqs) {
  dims <- dim(q_paths)
  if (length(dims) != 3L) stop("q_paths must be a three-dimensional array.")
  n_group <- dims[1L]
  n_node <- dims[2L]
  q_flat <- matrix(q_paths, nrow = n_group * n_node,
                   ncol = dims[3L])
  means <- array(NA_real_, dim = c(n_group, n_node, length(freqs)))
  for (m in seq_along(freqs)) {
    z <- -0.5 * beta ^ 2 * freqs[m] ^ 2 * q_flat
    finite_row <- rowSums(is.finite(z)) > 0L
    row_max <- rep(NA_real_, nrow(z))
    row_max[finite_row] <- apply(z[finite_row, , drop = FALSE],
                                 1L, max, na.rm = TRUE)
    value <- rep(NA_real_, nrow(z))
    if (any(finite_row)) {
      centered <- z[finite_row, , drop = FALSE] - row_max[finite_row]
      value[finite_row] <- row_max[finite_row] +
        log(rowMeans(exp(centered), na.rm = TRUE))
    }
    means[, , m] <- matrix(value, nrow = n_group, ncol = n_node)
  }
  means
}

make_exact_cf_terms <- function(observations) {
  lapply(seq_len(observations$n_blocks), function(j) {
    variance <- pmax(diag(observations$obs_cov[j, , ]),
                     .Machine$double.eps)
    list(y = observations$y_dejumped_bias_corrected[j, ],
         variance = variance,
         log_variance = sum(log(variance)),
         p = length(variance))
  })
}

exact_cf_node_loglik <- function(terms, means, dispersion) {
  if (!is.finite(dispersion) || dispersion <= 0) {
    return(matrix(-Inf, nrow(means), ncol(means)))
  }
  out <- matrix(0, nrow(means), ncol(means))
  for (m in seq_len(dim(means)[3L])) {
    out <- out + (terms$y[m] - means[, , m]) ^ 2 / terms$variance[m]
  }
  -0.5 * (terms$p * log(2 * pi) + terms$log_variance +
            terms$p * log(dispersion) + out / dispersion)
}

exact_cf_emission_table <- function(terms_list, means, dispersion) {
  if (!is.finite(dispersion) || dispersion <= 0) {
    return(matrix(-Inf, nrow = length(terms_list),
                  ncol = dim(means)[1L] * dim(means)[2L]))
  }
  n_group <- dim(means)[1L]
  n_node <- dim(means)[2L]
  p <- dim(means)[3L]
  mu <- matrix(means, nrow = n_group * n_node, ncol = p)
  y <- do.call(rbind, lapply(terms_list, `[[`, "y"))
  variance <- do.call(rbind, lapply(terms_list, `[[`, "variance"))
  inverse_variance <- 1 / variance
  y_weighted <- y * inverse_variance
  base_qf <- rowSums(y ^ 2 * inverse_variance)
  qf <- base_qf - 2 * (y_weighted %*% t(mu)) +
    inverse_variance %*% t(mu ^ 2)
  log_normalizer <- -0.5 * (p * log(2 * pi) +
    rowSums(log(variance)) + p * log(dispersion))
  sweep(-0.5 * qf / dispersion, 1L, log_normalizer, "+")
}

aggregate_exact_cf_nodes <- function(node_loglik, n_groups, node_count,
                                     node_weights = NULL) {
  out <- matrix(-Inf, nrow(node_loglik), n_groups)
  for (i in seq_len(n_groups)) {
    n_i <- node_count[i]
    if (n_i <= 0L) next
    columns <- i + (seq_len(n_i) - 1L) * n_groups
    values <- node_loglik[, columns, drop = FALSE]
    weights <- if (is.null(node_weights)) rep(1 / n_i, n_i) else
      as.numeric(node_weights[i, seq_len(n_i)])
    valid <- is.finite(weights) & weights > 0
    if (!any(valid)) next
    weights <- weights[valid] / sum(weights[valid])
    values <- values[, valid, drop = FALSE] +
      rep(log(weights), each = nrow(values))
    row_max <- apply(values, 1L, max)
    out[, i] <- row_max +
      log(rowSums(exp(values - row_max)))
  }
  out
}

exact_cf_conditional_moments <- function(q_paths, beta, freqs, alpha,
                                         sigma2, dt_reference, block_size,
                                         full_covariance = FALSE) {
  dims <- dim(q_paths)
  n_group <- dims[1L]
  n_node <- dims[2L]
  p <- length(freqs)
  q_flat <- matrix(q_paths, nrow = n_group * n_node,
                   ncol = dims[3L])
  valid <- rowSums(is.finite(q_flat)) > 0L
  mean_matrix <- variance_matrix <- matrix(
    NA_real_, nrow(q_flat), p
  )
  covariance_array <- if (isTRUE(full_covariance)) {
    array(NA_real_, dim = c(nrow(q_flat), p, p))
  } else NULL
  mean_phi_matrix <- matrix(NA_real_, nrow(q_flat), p)
  q <- q_flat[valid, , drop = FALSE]
  use_phi_cache <- isTRUE(getOption("exact_cf.use_phi_cache", TRUE))
  phi_cache <- new.env(hash = TRUE, parent = emptyenv())
  phi_at <- function(u) {
    key <- if (use_phi_cache) sprintf("%a", u) else ""
    if (use_phi_cache && exists(key, envir = phi_cache, inherits = FALSE)) {
      return(get(key, envir = phi_cache, inherits = FALSE))
    }
    jump_u <- sigma2 ^ alpha * abs(u) ^ alpha *
      dt_reference ^ (1 - alpha / 2)
    value <- exp(-0.5 * beta ^ 2 * u ^ 2 * q - jump_u)
    if (use_phi_cache) assign(key, value, envir = phi_cache)
    value
  }
  phi_by_freq <- vector("list", p)
  for (m in seq_along(freqs)) {
    u <- freqs[m]
    jump_u <- sigma2 ^ alpha * abs(u) ^ alpha *
      dt_reference ^ (1 - alpha / 2)
    phi_u <- phi_at(u)
    phi_2u <- phi_at(2 * u)
    phi_by_freq[[m]] <- phi_u
    mean_phi <- rowMeans(phi_u)
    variance_real <- rowMeans(0.5 * (1 + phi_2u) - phi_u ^ 2) /
      block_size
    variance_imag <- rowMeans(0.5 * (1 - phi_2u)) / block_size
    variance_log <- variance_real /
      pmax(mean_phi ^ 2, .Machine$double.eps)
    log_modulus_bias <- 0.5 * (variance_imag - variance_real) /
      pmax(mean_phi ^ 2, .Machine$double.eps)
    mean_matrix[valid, m] <- log(mean_phi) + jump_u + log_modulus_bias
    variance_matrix[valid, m] <- pmax(variance_log,
                                      .Machine$double.eps)
    mean_phi_matrix[valid, m] <- mean_phi
  }
  if (isTRUE(full_covariance)) for (m in seq_len(p)) {
    for (k in seq_len(m)) {
      u <- freqs[m]
      v <- freqs[k]
      phi_diff <- if (abs(u - v) < .Machine$double.eps) {
        matrix(1, nrow(q), ncol(q))
      } else {
        phi_at(abs(u - v))
      }
      phi_sum <- phi_at(u + v)
      covariance_real <- rowMeans(
        0.5 * (phi_diff + phi_sum) -
          phi_by_freq[[m]] * phi_by_freq[[k]]
      ) / block_size
      covariance_log <- covariance_real /
        pmax(mean_phi_matrix[valid, m] * mean_phi_matrix[valid, k],
             .Machine$double.eps)
      covariance_array[valid, m, k] <- covariance_log
      covariance_array[valid, k, m] <- covariance_log
    }
  }
  list(
    mean = array(mean_matrix, dim = c(n_group, n_node, p)),
    variance = array(variance_matrix, dim = c(n_group, n_node, p)),
    covariance = if (isTRUE(full_covariance)) {
      array(covariance_array, dim = c(n_group, n_node, p, p))
    } else NULL
  )
}

conditional_exact_cf_emission_table <- function(y_mat, moments) {
  dims <- dim(moments$mean)
  n_group <- dims[1L]
  n_node <- dims[2L]
  p <- dims[3L]
  mu <- matrix(moments$mean, nrow = n_group * n_node, ncol = p)
  variance <- matrix(moments$variance,
                     nrow = n_group * n_node, ncol = p)
  inverse_variance <- 1 / variance
  qf <- y_mat ^ 2 %*% t(inverse_variance) -
    2 * y_mat %*% t(mu * inverse_variance) +
    matrix(rowSums(mu ^ 2 * inverse_variance), nrow(y_mat),
           nrow(mu), byrow = TRUE)
  log_normalizer <- -0.5 * (p * log(2 * pi) + rowSums(log(variance)))
  -0.5 * qf + rep(log_normalizer, each = nrow(qf))
}

conditional_exact_cf_emission_table_full <- function(y_mat, moments) {
  dims <- dim(moments$mean)
  n_group <- dims[1L]
  n_node <- dims[2L]
  p <- dims[3L]
  if (p != 3L) {
    stop("The vectorized full-covariance emission requires 3 frequencies.")
  }
  mu <- matrix(moments$mean, nrow = n_group * n_node, ncol = p)
  covariance <- matrix(moments$covariance,
                       nrow = n_group * n_node, ncol = p * p)
  scale <- pmax(covariance[, 1L], covariance[, 5L], covariance[, 9L],
                .Machine$double.eps)
  jitter <- scale * 1e-10
  a <- covariance[, 1L] + jitter
  b <- 0.5 * (covariance[, 4L] + covariance[, 2L])
  c <- 0.5 * (covariance[, 7L] + covariance[, 3L])
  d <- covariance[, 5L] + jitter
  e <- 0.5 * (covariance[, 8L] + covariance[, 6L])
  f <- covariance[, 9L] + jitter
  determinant <- a * d * f + 2 * b * c * e - a * e ^ 2 -
    d * c ^ 2 - f * b ^ 2
  invalid <- !is.finite(determinant) | determinant <=
    .Machine$double.eps * scale ^ 3
  if (any(invalid)) {
    extra <- scale[invalid] * 1e-6
    a[invalid] <- a[invalid] + extra
    d[invalid] <- d[invalid] + extra
    f[invalid] <- f[invalid] + extra
    determinant[invalid] <- a[invalid] * d[invalid] * f[invalid] +
      2 * b[invalid] * c[invalid] * e[invalid] -
      a[invalid] * e[invalid] ^ 2 - d[invalid] * c[invalid] ^ 2 -
      f[invalid] * b[invalid] ^ 2
  }
  determinant <- pmax(determinant, .Machine$double.xmin)
  inverse_a <- (d * f - e ^ 2) / determinant
  inverse_b <- (c * e - b * f) / determinant
  inverse_c <- (b * e - c * d) / determinant
  inverse_d <- (a * f - c ^ 2) / determinant
  inverse_e <- (b * c - a * e) / determinant
  inverse_f <- (a * d - b ^ 2) / determinant

  precision_mu_1 <- inverse_a * mu[, 1L] + inverse_b * mu[, 2L] +
    inverse_c * mu[, 3L]
  precision_mu_2 <- inverse_b * mu[, 1L] + inverse_d * mu[, 2L] +
    inverse_e * mu[, 3L]
  precision_mu_3 <- inverse_c * mu[, 1L] + inverse_e * mu[, 2L] +
    inverse_f * mu[, 3L]
  mu_qf <- mu[, 1L] * precision_mu_1 + mu[, 2L] * precision_mu_2 +
    mu[, 3L] * precision_mu_3

  y1 <- y_mat[, 1L]
  y2 <- y_mat[, 2L]
  y3 <- y_mat[, 3L]
  qf <- tcrossprod(y1 ^ 2, inverse_a)
  qf <- qf + tcrossprod(y2 ^ 2, inverse_d)
  qf <- qf + tcrossprod(y3 ^ 2, inverse_f)
  qf <- qf + 2 * tcrossprod(y1 * y2, inverse_b)
  qf <- qf + 2 * tcrossprod(y1 * y3, inverse_c)
  qf <- qf + 2 * tcrossprod(y2 * y3, inverse_e)
  qf <- qf - 2 * tcrossprod(y1, precision_mu_1)
  qf <- qf - 2 * tcrossprod(y2, precision_mu_2)
  qf <- qf - 2 * tcrossprod(y3, precision_mu_3)
  qf <- sweep(qf, 2L, mu_qf, "+")
  log_normalizer <- -0.5 * (p * log(2 * pi) + log(determinant))
  sweep(-0.5 * qf, 2L, log_normalizer, "+")
}

stable_gaussian_three_emission_table <- function(y_mat, mean_array,
                                                  covariance_array) {
  dims <- dim(mean_array)
  n_group <- dims[1L]
  n_node <- dims[2L]
  mu <- matrix(mean_array, nrow = n_group * n_node, ncol = 3L)
  covariance <- matrix(covariance_array,
                       nrow = n_group * n_node, ncol = 9L)
  a <- covariance[, 1L]
  b <- 0.5 * (covariance[, 4L] + covariance[, 2L])
  c <- 0.5 * (covariance[, 7L] + covariance[, 3L])
  d <- covariance[, 5L]
  e <- 0.5 * (covariance[, 8L] + covariance[, 6L])
  f <- covariance[, 9L]
  scale <- pmax(a, d, f, .Machine$double.eps)
  numerical_floor <- sqrt(.Machine$double.eps) * scale
  a <- a + numerical_floor
  d <- d + numerical_floor
  f <- f + numerical_floor
  determinant <- a * d * f + 2 * b * c * e - a * e ^ 2 -
    d * c ^ 2 - f * b ^ 2
  second_minor <- a * d - b ^ 2
  invalid <- is.finite(scale) & (!is.finite(determinant) |
    determinant <= 0 | !is.finite(second_minor) | second_minor <= 0)
  if (any(invalid)) {
    for (index in which(invalid)) {
      covariance_i <- matrix(c(
        a[index], b[index], c[index],
        b[index], d[index], e[index],
        c[index], e[index], f[index]
      ), 3L, 3L, byrow = TRUE)
      eigenvalues <- eigen(covariance_i, symmetric = TRUE,
                           only.values = TRUE)$values
      loading <- max(0, -min(eigenvalues)) +
        sqrt(.Machine$double.eps) * max(abs(eigenvalues),
                                        .Machine$double.eps)
      a[index] <- a[index] + loading
      d[index] <- d[index] + loading
      f[index] <- f[index] + loading
    }
    determinant <- a * d * f + 2 * b * c * e - a * e ^ 2 -
      d * c ^ 2 - f * b ^ 2
  }
  determinant[is.finite(determinant)] <- pmax(
    determinant[is.finite(determinant)], .Machine$double.xmin
  )
  inverse_a <- (d * f - e ^ 2) / determinant
  inverse_b <- (c * e - b * f) / determinant
  inverse_c <- (b * e - c * d) / determinant
  inverse_d <- (a * f - c ^ 2) / determinant
  inverse_e <- (b * c - a * e) / determinant
  inverse_f <- (a * d - b ^ 2) / determinant
  precision_mu_1 <- inverse_a * mu[, 1L] + inverse_b * mu[, 2L] +
    inverse_c * mu[, 3L]
  precision_mu_2 <- inverse_b * mu[, 1L] + inverse_d * mu[, 2L] +
    inverse_e * mu[, 3L]
  precision_mu_3 <- inverse_c * mu[, 1L] + inverse_e * mu[, 2L] +
    inverse_f * mu[, 3L]
  mu_qf <- mu[, 1L] * precision_mu_1 + mu[, 2L] * precision_mu_2 +
    mu[, 3L] * precision_mu_3
  y1 <- y_mat[, 1L]
  y2 <- y_mat[, 2L]
  y3 <- y_mat[, 3L]
  qf <- tcrossprod(y1 ^ 2, inverse_a) +
    tcrossprod(y2 ^ 2, inverse_d) + tcrossprod(y3 ^ 2, inverse_f) +
    2 * tcrossprod(y1 * y2, inverse_b) +
    2 * tcrossprod(y1 * y3, inverse_c) +
    2 * tcrossprod(y2 * y3, inverse_e) -
    2 * tcrossprod(y1, precision_mu_1) -
    2 * tcrossprod(y2, precision_mu_2) -
    2 * tcrossprod(y3, precision_mu_3)
  qf <- sweep(qf, 2L, mu_qf, "+")
  log_normalizer <- -0.5 * (3L * log(2 * pi) + log(determinant))
  out <- sweep(-0.5 * qf, 2L, log_normalizer, "+")
  out[!is.finite(out)] <- -Inf
  out
}

complex_ecf_conditional_moments <- function(
    q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size) {
  dims <- dim(q_paths)
  n_group <- dims[1L]
  n_node <- dims[2L]
  p <- length(freqs)
  if (p != 3L) stop("The complex-CF model requires three frequencies.")
  q_flat <- matrix(q_paths, nrow = n_group * n_node,
                   ncol = dims[3L])
  valid <- rowSums(is.finite(q_flat)) > 0L
  q <- q_flat[valid, , drop = FALSE]
  phi_at <- function(u) {
    jump_u <- sigma2 ^ alpha * abs(u) ^ alpha *
      dt_reference ^ (1 - alpha / 2)
    exp(-0.5 * beta ^ 2 * u ^ 2 * q - jump_u)
  }
  phi <- lapply(freqs, phi_at)
  mean_real <- matrix(NA_real_, nrow(q_flat), p)
  mean_real[valid, ] <- vapply(phi, rowMeans, numeric(nrow(q)))
  covariance_real <- covariance_imag <- array(
    NA_real_, dim = c(nrow(q_flat), p, p)
  )
  for (m in seq_len(p)) {
    for (k in seq_len(m)) {
      phi_diff <- if (abs(freqs[m] - freqs[k]) < .Machine$double.eps) {
        matrix(1, nrow(q), ncol(q))
      } else phi_at(abs(freqs[m] - freqs[k]))
      phi_sum <- phi_at(freqs[m] + freqs[k])
      cov_re <- rowMeans(0.5 * (phi_diff + phi_sum) -
                           phi[[m]] * phi[[k]]) / block_size
      cov_im <- rowMeans(0.5 * (phi_diff - phi_sum)) / block_size
      covariance_real[valid, m, k] <- cov_re
      covariance_real[valid, k, m] <- cov_re
      covariance_imag[valid, m, k] <- cov_im
      covariance_imag[valid, k, m] <- cov_im
    }
  }
  list(
    mean_real = array(mean_real, dim = c(n_group, n_node, p)),
    covariance_real = array(covariance_real,
                            dim = c(n_group, n_node, p, p)),
    covariance_imag = array(covariance_imag,
                            dim = c(n_group, n_node, p, p))
  )
}

complex_ecf_emission_table <- function(phi_hat, moments) {
  real_loglik <- stable_gaussian_three_emission_table(
    Re(phi_hat), moments$mean_real, moments$covariance_real
  )
  zero_mean <- moments$mean_real * 0
  imaginary_loglik <- stable_gaussian_three_emission_table(
    Im(phi_hat), zero_mean, moments$covariance_imag
  )
  real_loglik + imaginary_loglik
}

complex_ecf_forward_loglik <- function(
    phi_hat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    return_filter = FALSE) {
  pair_moments <- complex_ecf_conditional_moments(
    model$pair_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size
  )
  initial_moments <- complex_ecf_conditional_moments(
    model$initial_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size
  )
  pair_node_ll <- complex_ecf_emission_table(phi_hat, pair_moments)
  pair_ll <- aggregate_exact_cf_nodes(
    pair_node_ll, nrow(model$active_pairs), model$pair_node_count,
    model$pair_node_weights
  )
  initial_node_ll <- complex_ecf_emission_table(
    phi_hat[1L, , drop = FALSE], initial_moments
  )
  initial_ll <- as.numeric(aggregate_exact_cf_nodes(
    initial_node_ll, model$n_states, model$initial_node_count,
    model$initial_node_weights
  ))
  n_obs <- nrow(phi_hat)
  log_transition <- log(model$transition)
  log_transition[!is.finite(log_transition)] <- -Inf
  first <- log(model$initial_transition) + initial_ll
  first[!is.finite(first)] <- -Inf
  increments <- numeric(n_obs)
  increments[1L] <- log_sum_exp(first)
  alpha_filter <- first - increments[1L]
  log_alpha <- if (return_filter) {
    matrix(NA_real_, n_obs, model$n_states)
  } else NULL
  if (return_filter) log_alpha[1L, ] <- alpha_filter
  if (n_obs >= 2L) for (j in 2:n_obs) {
    edge <- matrix(-Inf, model$n_states, model$n_states)
    edge[model$active_pairs] <- log_transition[model$active_pairs] +
      pair_ll[j, ]
    joint <- edge + alpha_filter
    next_alpha <- apply(joint, 2L, log_sum_exp)
    increments[j] <- log_sum_exp(next_alpha)
    alpha_filter <- next_alpha - increments[j]
    if (return_filter) log_alpha[j, ] <- alpha_filter
  }
  result <- list(loglik = sum(increments))
  if (return_filter) {
    filter_prob <- exp(log_alpha)
    filter_prob[!is.finite(filter_prob)] <- 0
    filter_prob <- filter_prob / pmax(rowSums(filter_prob),
                                      .Machine$double.eps)
    entropy <- -rowSums(filter_prob * log(pmax(filter_prob,
                                               .Machine$double.xmin)))
    result$filter_diagnostics <- list(
      filter_prob = filter_prob,
      loglik_increment = increments,
      mean_filter_max_prob = mean(apply(filter_prob, 1L, max)),
      mean_filter_entropy = mean(entropy) / log(model$n_states),
      mean_filter_boundary_prob = mean(filter_prob[, 1L] +
                                         filter_prob[, model$n_states])
    )
  }
  result
}

complex_ecf_known_endpoint_loglik <- function(
    phi_hat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    start_state, end_state) {
  pair_moments <- complex_ecf_conditional_moments(
    model$pair_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size
  )
  initial_moments <- complex_ecf_conditional_moments(
    model$initial_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size
  )
  pair_node_ll <- complex_ecf_emission_table(phi_hat, pair_moments)
  pair_ll <- aggregate_exact_cf_nodes(
    pair_node_ll, nrow(model$active_pairs), model$pair_node_count,
    model$pair_node_weights
  )
  initial_node_ll <- complex_ecf_emission_table(
    phi_hat[1L, , drop = FALSE], initial_moments
  )
  initial_ll <- as.numeric(aggregate_exact_cf_nodes(
    initial_node_ll, model$n_states, model$initial_node_count,
    model$initial_node_weights
  ))
  first_end <- end_state[1L]
  if (!is.finite(first_end) || model$initial_transition[first_end] <= 0) {
    return(-Inf)
  }
  loglik <- log(model$initial_transition[first_end]) + initial_ll[first_end]
  pair_lookup <- matrix(NA_integer_, model$n_states, model$n_states)
  pair_lookup[model$active_pairs] <- seq_len(nrow(model$active_pairs))
  if (nrow(phi_hat) >= 2L) for (j in 2:nrow(phi_hat)) {
    a <- start_state[j]
    b <- end_state[j]
    pair_row <- pair_lookup[a, b]
    probability <- model$transition[a, b]
    if (!is.finite(pair_row) || !is.finite(probability) || probability <= 0) {
      return(-Inf)
    }
    loglik <- loglik + log(probability) + pair_ll[j, pair_row]
  }
  loglik
}

conditional_exact_cf_forward_loglik <- function(
    y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    return_filter = FALSE, full_covariance = FALSE) {
  pair_moments <- exact_cf_conditional_moments(
    model$pair_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size,
    full_covariance = full_covariance
  )
  initial_moments <- exact_cf_conditional_moments(
    model$initial_q_paths, beta, freqs, alpha, sigma2, dt_reference,
    block_size, full_covariance = full_covariance
  )
  emission_function <- if (isTRUE(full_covariance)) {
    conditional_exact_cf_emission_table_full
  } else {
    conditional_exact_cf_emission_table
  }
  pair_node_ll <- emission_function(y_mat, pair_moments)
  pair_ll <- aggregate_exact_cf_nodes(
    pair_node_ll, nrow(model$active_pairs), model$pair_node_count,
    model$pair_node_weights
  )
  initial_node_ll <- emission_function(
    y_mat[1L, , drop = FALSE], initial_moments
  )
  initial_ll <- as.numeric(aggregate_exact_cf_nodes(
    initial_node_ll, model$n_states, model$initial_node_count,
    model$initial_node_weights
  ))

  n_obs <- nrow(y_mat)
  n_states <- model$n_states
  log_transition <- log(model$transition)
  log_transition[!is.finite(log_transition)] <- -Inf
  first <- log(model$initial_transition) + initial_ll
  first[!is.finite(first)] <- -Inf
  increments <- numeric(n_obs)
  increments[1L] <- log_sum_exp(first)
  alpha_filter <- first - increments[1L]
  log_alpha <- if (return_filter) matrix(NA_real_, n_obs, n_states) else NULL
  if (return_filter) log_alpha[1L, ] <- alpha_filter
  if (n_obs >= 2L) {
    for (j in 2:n_obs) {
      edge <- matrix(-Inf, n_states, n_states)
      edge[model$active_pairs] <- log_transition[model$active_pairs] +
        pair_ll[j, ]
      joint <- edge + alpha_filter
      next_alpha <- apply(joint, 2L, log_sum_exp)
      increments[j] <- log_sum_exp(next_alpha)
      alpha_filter <- next_alpha - increments[j]
      if (return_filter) log_alpha[j, ] <- alpha_filter
    }
  }
  result <- list(loglik = sum(increments))
  if (return_filter) {
    filter_prob <- exp(log_alpha)
    filter_prob[!is.finite(filter_prob)] <- 0
    row_total <- rowSums(filter_prob)
    valid <- row_total > 0
    filter_prob[valid, ] <- filter_prob[valid, , drop = FALSE] /
      row_total[valid]
    entropy <- -rowSums(filter_prob * log(pmax(filter_prob,
                                               .Machine$double.xmin)))
    result$filter_diagnostics <- list(
      filter_prob = filter_prob,
      loglik_increment = increments,
      mean_filter_max_prob = mean(apply(filter_prob[valid, , drop = FALSE],
                                        1L, max)),
      mean_filter_entropy = mean(entropy[valid]) / log(n_states),
      mean_filter_boundary_prob = mean(filter_prob[valid, 1L] +
                                         filter_prob[valid, n_states])
    )
  }
  result
}

global_profile_optimize <- function(fn, interval, n_grid = 25L,
                                    tol = 1e-4) {
  n_grid <- max(9L, as.integer(n_grid))
  grid <- seq(interval[1L], interval[2L], length.out = n_grid)
  values <- vapply(grid, fn, numeric(1L))
  finite <- is.finite(values)
  if (!any(finite)) stop("No finite values in one-dimensional profile.")
  best_grid <- which.min(ifelse(finite, values, Inf))
  candidates <- data.frame(
    parameter = grid[best_grid], objective = values[best_grid]
  )
  evaluations <- n_grid
  if (best_grid > 1L && best_grid < n_grid) {
    counted_fn <- function(x) {
      evaluations <<- evaluations + 1L
      fn(x)
    }
    local <- stats::optimize(
      counted_fn, interval = grid[c(best_grid - 1L, best_grid + 1L)],
      tol = tol
    )
    candidates <- rbind(
      candidates,
      data.frame(parameter = local$minimum, objective = local$objective)
    )
  }
  best <- candidates[which.min(candidates$objective), , drop = FALSE]
  list(
    minimum = best$parameter,
    objective = best$objective,
    grid = grid,
    grid_objective = values,
    evaluations = evaluations
  )
}

exact_cf_forward_loglik <- function(terms_list, model, beta, freqs,
                                    dispersion, return_filter = FALSE) {
  pair_means <- exact_cf_path_means(model$pair_q_paths, beta, freqs)
  initial_means <- exact_cf_path_means(model$initial_q_paths, beta, freqs)
  n_obs <- length(terms_list)
  n_states <- model$n_states
  log_transition <- log(model$transition)
  log_transition[!is.finite(log_transition)] <- -Inf
  log_alpha <- if (return_filter) matrix(NA_real_, n_obs, n_states) else NULL
  increments <- numeric(n_obs)

  initial_node_ll <- exact_cf_emission_table(terms_list[1L], initial_means,
                                              dispersion)
  initial_ll <- as.numeric(aggregate_exact_cf_nodes(
    initial_node_ll, n_states, model$initial_node_count,
    model$initial_node_weights
  ))
  first <- log(model$initial_transition) + initial_ll
  first[!is.finite(first)] <- -Inf
  increments[1L] <- log_sum_exp(first)
  alpha <- first - increments[1L]
  if (return_filter) log_alpha[1L, ] <- alpha

  if (n_obs >= 2L) {
    pair_node_ll <- exact_cf_emission_table(terms_list, pair_means,
                                            dispersion)
    pair_ll <- aggregate_exact_cf_nodes(
      pair_node_ll, nrow(model$active_pairs), model$pair_node_count,
      model$pair_node_weights
    )
    for (j in 2:n_obs) {
      edge <- matrix(-Inf, n_states, n_states)
      edge[model$active_pairs] <- log_transition[model$active_pairs] +
        pair_ll[j, ]
      joint <- edge + alpha
      next_alpha <- apply(joint, 2L, log_sum_exp)
      increments[j] <- log_sum_exp(next_alpha)
      alpha <- next_alpha - increments[j]
      if (return_filter) log_alpha[j, ] <- alpha
    }
  }
  result <- list(loglik = sum(increments))
  if (return_filter) {
    filter_prob <- exp(log_alpha)
    filter_prob[!is.finite(filter_prob)] <- 0
    row_total <- rowSums(filter_prob)
    valid <- row_total > 0
    filter_prob[valid, ] <- filter_prob[valid, , drop = FALSE] / row_total[valid]
    entropy <- -rowSums(filter_prob * log(pmax(filter_prob,
                                               .Machine$double.xmin)))
    result$filter_diagnostics <- list(
      filter_prob = filter_prob,
      loglik_increment = increments,
      mean_filter_max_prob = mean(apply(filter_prob[valid, , drop = FALSE],
                                        1L, max)),
      mean_filter_entropy = mean(entropy[valid]) / log(n_states),
      mean_filter_boundary_prob = mean(filter_prob[valid, 1L] +
                                         filter_prob[valid, n_states])
    )
  }
  result
}

exact_cf_known_endpoint_loglik <- function(terms_list, model, beta, freqs,
                                           dispersion, start_state,
                                           end_state) {
  pair_means <- exact_cf_path_means(model$pair_q_paths, beta, freqs)
  initial_means <- exact_cf_path_means(model$initial_q_paths, beta, freqs)
  initial_node_ll <- exact_cf_emission_table(terms_list[1L], initial_means,
                                              dispersion)
  initial_ll <- as.numeric(aggregate_exact_cf_nodes(
    initial_node_ll, model$n_states, model$initial_node_count,
    model$initial_node_weights
  ))
  first_end <- end_state[1L]
  if (!is.finite(first_end) || model$initial_transition[first_end] <= 0) {
    return(-Inf)
  }
  loglik <- log(model$initial_transition[first_end]) + initial_ll[first_end]
  if (length(terms_list) == 1L) return(loglik)

  pair_node_ll <- exact_cf_emission_table(terms_list, pair_means, dispersion)
  pair_ll <- aggregate_exact_cf_nodes(
    pair_node_ll, nrow(model$active_pairs), model$pair_node_count,
    model$pair_node_weights
  )
  pair_lookup <- matrix(NA_integer_, model$n_states, model$n_states)
  pair_lookup[model$active_pairs] <- seq_len(nrow(model$active_pairs))
  for (j in 2:length(terms_list)) {
    a <- start_state[j]
    b <- end_state[j]
    pair_row <- pair_lookup[a, b]
    probability <- model$transition[a, b]
    if (!is.finite(pair_row) || !is.finite(probability) || probability <= 0) {
      return(-Inf)
    }
    loglik <- loglik + log(probability) + pair_ll[j, pair_row]
  }
  loglik
}

conditional_exact_cf_known_endpoint_loglik <- function(
    y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    start_state, end_state, full_covariance = FALSE) {
  pair_moments <- exact_cf_conditional_moments(
    model$pair_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size,
    full_covariance = full_covariance
  )
  initial_moments <- exact_cf_conditional_moments(
    model$initial_q_paths, beta, freqs, alpha, sigma2, dt_reference,
    block_size, full_covariance = full_covariance
  )
  emission_function <- if (isTRUE(full_covariance)) {
    conditional_exact_cf_emission_table_full
  } else {
    conditional_exact_cf_emission_table
  }
  initial_node_ll <- emission_function(
    y_mat[1L, , drop = FALSE], initial_moments
  )
  initial_ll <- as.numeric(aggregate_exact_cf_nodes(
    initial_node_ll, model$n_states, model$initial_node_count,
    model$initial_node_weights
  ))
  first_end <- end_state[1L]
  if (!is.finite(first_end) || model$initial_transition[first_end] <= 0) {
    return(-Inf)
  }
  loglik <- log(model$initial_transition[first_end]) + initial_ll[first_end]
  if (nrow(y_mat) == 1L) return(loglik)

  pair_node_ll <- emission_function(y_mat, pair_moments)
  pair_ll <- aggregate_exact_cf_nodes(
    pair_node_ll, nrow(model$active_pairs), model$pair_node_count,
    model$pair_node_weights
  )
  pair_lookup <- matrix(NA_integer_, model$n_states, model$n_states)
  pair_lookup[model$active_pairs] <- seq_len(nrow(model$active_pairs))
  for (j in 2:nrow(y_mat)) {
    a <- start_state[j]
    b <- end_state[j]
    pair_row <- pair_lookup[a, b]
    probability <- model$transition[a, b]
    if (!is.finite(pair_row) || !is.finite(probability) || probability <= 0) {
      return(-Inf)
    }
    loglik <- loglik + log(probability) + pair_ll[j, pair_row]
  }
  loglik
}

estimate_exact_cf_endpoint_hmm <- function(
    R, dt, alpha_hat, sigma2_hat,
    block_horizon = 1 / 20,
    rho_grid = seq(0.7, 2.3, by = 0.1),
    beta_bounds = c(0.20, 5.00),
    dispersion_bounds = c(0.25, Inf),
    freq_mults = c(0.25, 0.50, 0.75),
    n_states = 21L,
    transition_sim_blocks = 630000L,
    substeps_per_block = 32L,
    max_path_nodes_per_pair = 31L,
    kernel_seed = 5000L,
    beta_profile_grid_size = 25L,
    kernel_models = NULL,
    covariance_model = c("conditional_diag", "conditional_full",
                         "conditional_complex_full", "empirical_profiled"),
    confidence_level = 0.95,
    checkpoint_dir = NULL,
    checkpoint_signature = NULL,
    return_observations = FALSE,
    verbose = TRUE) {
  covariance_model <- match.arg(covariance_model)
  if (covariance_model %in% c("conditional_full",
                              "conditional_complex_full") &&
      length(freq_mults) != 3L) {
    stop("Full conditional covariance models require exactly three frequencies.")
  }
  n <- min(length(R), length(dt))
  R <- as.numeric(R[seq_len(n)])
  dt <- as.numeric(dt[seq_len(n)])
  dt_reference <- stats::median(dt[is.finite(dt) & dt > 0])
  block_size <- max(10L, as.integer(round(block_horizon / dt_reference)))
  actual_h <- block_size * dt_reference
  state_info <- make_q_endpoint_breaks(as.integer(n_states))
  if (is.null(kernel_models)) {
    if (verbose) {
      cat(sprintf("Building %d exact-CF endpoint kernels at H=%.6g...\n",
                  length(rho_grid), actual_h))
    }
    models <- lapply(rho_grid, function(rho) {
      simulate_q_endpoint_cf_model(
        rho = rho,
        block_h = actual_h,
        state_info = state_info,
        n_sim_blocks = as.integer(transition_sim_blocks),
        substeps_per_block = as.integer(substeps_per_block),
        max_path_nodes_per_pair = as.integer(max_path_nodes_per_pair),
        state_model = "conditional",
        seed = as.integer(kernel_seed)
      )
    })
    kernel_source <- "built"
  } else {
    models <- kernel_models
    model_rho <- vapply(models, `[[`, numeric(1L), "rho")
    if (length(models) != length(rho_grid) ||
        any(abs(model_rho - rho_grid) > 1e-10)) {
      stop("kernel_models do not match rho_grid.")
    }
    model_h <- vapply(models, `[[`, numeric(1L), "block_h")
    if (any(abs(model_h - actual_h) > 1e-10)) {
      stop("kernel_models do not match the requested block horizon.")
    }
    kernel_source <- "prebuilt"
    if (verbose) {
      cat(sprintf("Using %d prebuilt exact-CF endpoint kernels at H=%.6g...\n",
                  length(models), actual_h))
    }
  }
  observations <- block_cf_log_observations(
    R = R, dt = dt, alpha = alpha_hat, sigma2 = sigma2_hat,
    block_size = block_size, freq_mults = freq_mults
  )
  terms_list <- make_exact_cf_terms(observations)
  projected <- lapply(seq_len(observations$n_blocks), function(j) {
    emission_terms(observations$y_dejumped_bias_corrected[j, ],
                   observations$obs_cov[j, , ], observations$freqs,
                   emission_mode = "projected_diag")
  })
  projected_t <- vapply(projected, `[[`, numeric(1L), "t")
  beta_initial <- sqrt(max(stats::median(projected_t, na.rm = TRUE) /
                             stationary_q_quantile(0.5),
                           beta_bounds[1L] ^ 2))
  beta_initial <- min(max(beta_initial, beta_bounds[1L]), beta_bounds[2L])

  rows <- vector("list", length(rho_grid))
  fits <- vector("list", length(rho_grid))
  checkpoint_reused <- 0L
  checkpoint_enabled <- !is.null(checkpoint_dir) &&
    length(checkpoint_dir) == 1L && nzchar(checkpoint_dir)
  if (checkpoint_enabled) {
    if (is.null(checkpoint_signature) ||
        length(checkpoint_signature) != 1L ||
        !nzchar(checkpoint_signature)) {
      stop("checkpoint_signature is required when checkpoint_dir is used.")
    }
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(checkpoint_dir)) {
      stop("Could not create checkpoint directory: ", checkpoint_dir)
    }
  }
  for (i in seq_along(rho_grid)) {
    checkpoint_path <- if (checkpoint_enabled) {
      file.path(checkpoint_dir, sprintf(
        "rho_%03d.rds", as.integer(round(100 * rho_grid[i]))
      ))
    } else NULL
    if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
      checkpoint <- tryCatch(
        readRDS(checkpoint_path),
        error = function(e) {
          stop("Cannot read rho checkpoint ", checkpoint_path, ": ",
               conditionMessage(e))
        }
      )
      valid_checkpoint <- is.list(checkpoint) &&
        identical(checkpoint$version, 1L) &&
        identical(checkpoint$signature, checkpoint_signature) &&
        isTRUE(all.equal(checkpoint$rho, rho_grid[i], tolerance = 0)) &&
        is.data.frame(checkpoint$row) && nrow(checkpoint$row) == 1L &&
        is.list(checkpoint$fit) && is.finite(checkpoint$fit$loglik)
      if (!valid_checkpoint) {
        stop("Rho checkpoint does not match this exact evaluation: ",
             checkpoint_path)
      }
      rows[[i]] <- checkpoint$row
      fits[[i]] <- checkpoint$fit
      checkpoint_reused <- checkpoint_reused + 1L
      if (verbose) {
        cat(sprintf(
          "  rho[%02d/%02d]=%.3f beta=%.6f loglik=%.6f [checkpoint]\n",
          i, length(rho_grid), rho_grid[i], checkpoint$row$beta,
          checkpoint$row$loglik
        ))
        flush.console()
      }
      next
    }
    model <- models[[i]]
    if (covariance_model %in% c("conditional_diag", "conditional_full",
                                "conditional_complex_full")) {
      full_covariance <- identical(covariance_model, "conditional_full")
      opt <- global_profile_optimize(
        fn = function(log_beta) {
          value <- if (identical(covariance_model,
                                 "conditional_complex_full")) {
            complex_ecf_forward_loglik(
              observations$phi_hat, model, exp(log_beta),
              observations$freqs, alpha_hat, sigma2_hat, dt_reference,
              block_size
            )$loglik
          } else {
            conditional_exact_cf_forward_loglik(
              observations$y_dejumped, model, exp(log_beta),
              observations$freqs, alpha_hat, sigma2_hat, dt_reference,
              block_size, full_covariance = full_covariance
            )$loglik
          }
          if (is.finite(value)) -value else .Machine$double.xmax / 100
        },
        interval = log(beta_bounds), n_grid = beta_profile_grid_size,
        tol = 1e-4
      )
      beta_hat <- exp(opt$minimum)
      kappa_hat <- 1
      convergence <- 0L
      optimizer_iterations <- opt$evaluations
      fit <- if (identical(covariance_model, "conditional_complex_full")) {
        complex_ecf_forward_loglik(
          observations$phi_hat, model, beta_hat, observations$freqs,
          alpha_hat, sigma2_hat, dt_reference, block_size,
          return_filter = TRUE
        )
      } else {
        conditional_exact_cf_forward_loglik(
          observations$y_dejumped, model, beta_hat, observations$freqs,
          alpha_hat, sigma2_hat, dt_reference, block_size,
          return_filter = TRUE, full_covariance = full_covariance
        )
      }
    } else {
      opt <- stats::optim(
        par = log(c(beta_initial, 3)),
        fn = function(par) {
          beta <- exp(par[1L])
          dispersion <- exp(par[2L])
          value <- exact_cf_forward_loglik(
            terms_list, model, beta, observations$freqs, dispersion
          )$loglik
          if (is.finite(value)) -value else .Machine$double.xmax / 100
        },
        method = "L-BFGS-B",
        lower = log(c(beta_bounds[1L], dispersion_bounds[1L])),
        upper = log(c(beta_bounds[2L], dispersion_bounds[2L])),
        control = list(factr = 1e8, pgtol = 1e-6, maxit = 80L)
      )
      beta_hat <- exp(opt$par[1L])
      kappa_hat <- exp(opt$par[2L])
      convergence <- opt$convergence
      optimizer_iterations <- opt$counts[["function"]]
      fit <- exact_cf_forward_loglik(
        terms_list, model, beta_hat, observations$freqs, kappa_hat,
        return_filter = TRUE
      )
    }
    fits[[i]] <- fit
    rows[[i]] <- data.frame(
      sigma1 = beta_hat / rho_grid[i],
      rho = rho_grid[i],
      beta = beta_hat,
      emission_dispersion = kappa_hat,
      optimizer_convergence = convergence,
      optimizer_iterations = optimizer_iterations,
      optimizer_relative_change = NA_real_,
      loglik = fit$loglik
    )
    if (!is.null(checkpoint_path)) {
      temporary_path <- tempfile(
        pattern = paste0(basename(checkpoint_path), ".tmp-"),
        tmpdir = dirname(checkpoint_path)
      )
      saveRDS(list(
        version = 1L,
        signature = checkpoint_signature,
        rho = rho_grid[i],
        row = rows[[i]],
        fit = fits[[i]]
      ), temporary_path)
      if (!file.rename(temporary_path, checkpoint_path)) {
        unlink(temporary_path)
        stop("Could not commit rho checkpoint: ", checkpoint_path)
      }
    }
    if (verbose) {
      cat(sprintf("  rho[%02d/%02d]=%.3f beta=%.6f loglik=%.6f\n",
                  i, length(rho_grid), rho_grid[i], beta_hat, fit$loglik))
      flush.console()
    }
  }
  grid <- do.call(rbind, rows)
  best_index <- which.max(grid$loglik)
  refined <- refine_rho_profile_quadratic(grid)
  profile_diagnostics <- observed_profile_diagnostics(
    grid, confidence_level = confidence_level
  )
  warnings <- character(0)
  if (profile_diagnostics$rho_at_grid_boundary) {
    warnings <- c(warnings, "rho_profile_maximum_at_grid_boundary")
  }
  if (grid$optimizer_convergence[best_index] != 0L) {
    warnings <- c(warnings, "profile_optimizer_not_converged")
  }
  result <- list(
    alpha_hat = alpha_hat,
    sigma2_hat = sigma2_hat,
    sigma1_hat = refined$sigma1,
    rho_hat = refined$rho,
    beta_hat = refined$beta,
    kappa_hat = refined$dispersion,
    block_size = block_size,
    block_horizon = actual_h,
    n_blocks = observations$n_blocks,
    freq_mults = freq_mults,
    rho_grid = rho_grid,
    profile = grid,
    profile_diagnostics = profile_diagnostics,
    filter_diagnostics = fits[[best_index]]$filter_diagnostics,
    optimizer_convergence = grid$optimizer_convergence[best_index],
    optimizer_iterations = grid$optimizer_iterations[best_index],
    profile_quadratic_refined = refined$refined,
    profile_quadratic_curvature = refined$curvature,
    warnings = warnings,
    method = paste0("markov_additive_endpoint_exact_cf_hmm_",
                    covariance_model),
    covariance_model = covariance_model
  )
  result$beta_profile_grid_size <- as.integer(beta_profile_grid_size)
  result$kernel_source <- kernel_source
  result$checkpoint_rows_reused <- checkpoint_reused
  result$checkpoint_dir <- if (checkpoint_enabled) checkpoint_dir else ""
  if (return_observations) result$observations <- observations
  if (verbose) {
    cat(sprintf(
      "Exact-CF endpoint HMM: sigma1=%.6f rho=%.6f beta=%.6f kappa=%.4f\n",
      result$sigma1_hat, result$rho_hat, result$beta_hat, result$kappa_hat
    ))
  }
  result
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  option_value <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^--", name, "="), "", hit[1L])
  }
  seed <- as.integer(option_value("seed", 1301L))
  n_steps <- as.integer(option_value("n-steps", 600000L))
  terminal <- as.numeric(option_value("terminal", 20))
  fixed <- list(
    alpha = as.numeric(option_value("alpha", 1.25)),
    sigma2 = as.numeric(option_value("sigma2", 1.0)),
    sigma1 = as.numeric(option_value("sigma1", 1.2)),
    rho = as.numeric(option_value("rho", 1.45))
  )
  rho_grid <- as.numeric(strsplit(option_value(
    "rho-grid", "0.9,1.1,1.3,1.45,1.6,1.8,2.1"
  ), ",", fixed = TRUE)[[1L]])
  sim <- simulate_residual_euler(
    n_steps = n_steps, terminal = terminal,
    alpha = fixed$alpha, sigma2 = fixed$sigma2,
    sigma1 = fixed$sigma1, rho = fixed$rho,
    seed = seed, keep_v = FALSE
  )
  estimate <- estimate_exact_cf_endpoint_hmm(
    sim$R, sim$dt, fixed$alpha, fixed$sigma2,
    block_horizon = as.numeric(option_value("block-horizon", 1 / 20)),
    rho_grid = rho_grid,
    n_states = as.integer(option_value("n-states", 21L)),
    transition_sim_blocks = as.integer(option_value(
      "transition-sim-blocks", 630000L
    )),
    substeps_per_block = as.integer(option_value("substeps", 32L)),
    max_path_nodes_per_pair = as.integer(option_value("max-path-nodes", 31L)),
    kernel_seed = as.integer(option_value("kernel-seed", 5000L)),
    beta_profile_grid_size = as.integer(option_value(
      "beta-profile-grid-size", 25L
    )),
    covariance_model = if (any(args == "--empirical-covariance")) {
      "empirical_profiled"
    } else "conditional_diag"
  )
  result <- data.frame(
    seed = seed, n_steps = n_steps, terminal = terminal,
    alpha_true = fixed$alpha, sigma2_true = fixed$sigma2,
    sigma1_true = fixed$sigma1, rho_true = fixed$rho,
    sigma1_hat = estimate$sigma1_hat, rho_hat = estimate$rho_hat,
    beta_hat = estimate$beta_hat, kappa_hat = estimate$kappa_hat,
    sigma1_abs_err = abs(estimate$sigma1_hat - fixed$sigma1),
    rho_abs_err = abs(estimate$rho_hat - fixed$rho),
    beta_abs_err = abs(estimate$beta_hat - fixed$sigma1 * fixed$rho),
    profile_width = estimate$profile_diagnostics$rho_profile_width,
    covariance_model = estimate$covariance_model,
    kernel_seed = as.integer(option_value("kernel-seed", 5000L)),
    warning = paste(estimate$warnings, collapse = ";"),
    stringsAsFactors = FALSE
  )
  dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  path <- file.path("outputs", sprintf(
    "exact_cf_endpoint_hmm_seed%d_%s.rds", seed, stamp
  ))
  saveRDS(list(result = result, estimate = estimate), path)
  print(result, row.names = FALSE)
  cat("Wrote:", path, "\n")
}
