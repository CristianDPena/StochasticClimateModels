# Memory-bounded evaluation of the unchanged conditional-diagonal likelihood.

forward_from_emission_mixtures_chunked <- function(
    pair_ll, initial_ll, model, return_filter = FALSE) {
  n_obs <- nrow(pair_ll)
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
  if (n_obs >= 2L) for (j in 2:n_obs) {
    edge <- matrix(-Inf, n_states, n_states)
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
    row_total <- rowSums(filter_prob)
    valid <- row_total > 0
    filter_prob[valid, ] <- filter_prob[valid, , drop = FALSE] / row_total[valid]
    entropy <- -rowSums(filter_prob * log(pmax(
      filter_prob, .Machine$double.xmin
    )))
    result$filter_diagnostics <- list(
      filter_prob = filter_prob,
      loglik_increment = increments,
      mean_filter_max_prob = mean(apply(
        filter_prob[valid, , drop = FALSE], 1L, max
      )),
      mean_filter_entropy = mean(entropy[valid]) / log(n_states),
      mean_filter_boundary_prob = mean(filter_prob[valid, 1L] +
                                         filter_prob[valid, n_states])
    )
  }
  result
}

conditional_exact_cf_emission_aggregate_chunked_library <- function(
    y_mat, moments, node_count, node_weights = NULL,
    groups_per_chunk = 16L) {
  dims <- dim(moments$mean)
  n_groups <- dims[1L]
  n_nodes <- dims[2L]
  p <- dims[3L]
  groups_per_chunk <- max(1L, as.integer(groups_per_chunk))
  if (length(node_count) != n_groups) {
    stop("node_count does not match the moment groups.")
  }
  out <- matrix(-Inf, nrow(y_mat), n_groups)
  for (first in seq.int(1L, n_groups, by = groups_per_chunk)) {
    group <- seq.int(first, min(n_groups, first + groups_per_chunk - 1L))
    chunk_moments <- list(
      mean = array(
        moments$mean[group, , , drop = FALSE],
        dim = c(length(group), n_nodes, p)
      ),
      variance = array(
        moments$variance[group, , , drop = FALSE],
        dim = c(length(group), n_nodes, p)
      )
    )
    node_loglik <- conditional_exact_cf_emission_table(y_mat, chunk_moments)
    chunk_weights <- if (is.null(node_weights)) NULL else
      node_weights[group, , drop = FALSE]
    out[, group] <- aggregate_exact_cf_nodes(
      node_loglik, length(group), node_count[group], chunk_weights
    )
  }
  out
}

conditional_exact_cf_forward_loglik_chunked_library <- function(
    y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    groups_per_chunk = 16L, return_filter = FALSE) {
  pair_moments <- exact_cf_conditional_moments(
    model$pair_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size
  )
  initial_moments <- exact_cf_conditional_moments(
    model$initial_q_paths, beta, freqs, alpha, sigma2, dt_reference, block_size
  )
  pair_ll <- conditional_exact_cf_emission_aggregate_chunked_library(
    y_mat, pair_moments, model$pair_node_count, model$pair_node_weights,
    groups_per_chunk = groups_per_chunk
  )
  initial_ll <- as.numeric(
    conditional_exact_cf_emission_aggregate_chunked_library(
      y_mat[1L, , drop = FALSE], initial_moments,
      model$initial_node_count, model$initial_node_weights,
      groups_per_chunk = groups_per_chunk
    )
  )
  forward <- forward_from_emission_mixtures_chunked(
    pair_ll, initial_ll, model, return_filter = return_filter
  )
  forward$pair_ll <- pair_ll
  forward$initial_ll <- initial_ll
  forward
}
