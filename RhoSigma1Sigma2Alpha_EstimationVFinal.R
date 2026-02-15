# packages
library(yuima)
library(openxlsx)


#--------------------- Simulation parameters
# Terminal   <- 1
# n_steps    <- 200000
# alpha_true <- 1.75
# sigma2_true<- 1.2
# sigma1_true<- 0.8
# rho_true   <- 1.90
Uprime <- function(x) 16*x^3 - 16*x
timestep <- 0.1

#------------ Yuima model
# mod2D <- setModel(
#   drift = c("-(16*x^3 - 16*x) - sigma1*(rho^2)*v",
#             "-rho^2*v"),
#   diffusion = matrix(c("sigma1*rho*sqrt(1+v^2)",
#                        "rho*sqrt(1+v^2)"),
#                      nrow = 2, byrow = TRUE),
#   jump.coeff = matrix(c("sigma2","0"), nrow = 2),
#   measure.type = "code",
#   measure = list(df = "rstable(z, alpha, 0, 1, 0)"),
#   state.variable = c("x","v"),
#   solve.variable = c("x","v"),
#   xinit = c("x"=0,"v"=0)
# )
# 
# samp <- setSampling(Terminal = Terminal, n = n_steps)
# yu   <- setYuima(model = mod2D, sampling = samp)
# sim  <- simulate(yu, true.par = list(
#   sigma1 = sigma1_true, rho = rho_true,
#   sigma2 = sigma2_true, alpha = alpha_true
# ))
# z      <- yuima::get.zoo.data(sim)
# time   <- as.numeric(index(z[[1]]))
# X      <- as.numeric(z[[1]])
# V      <- as.numeric(z[[2]])

preX <- openxlsx::read.xlsx(
  "C:/Users/crist/OneDrive/School/Research/Long Research/GRIP-Data-Excel.xlsx",
  sheet = "GRIP_impurity",
  rows  = 20:216480,
  cols  = 3
)

X <- preX[[1]]

dt <- rep(timestep, length(X))
dX <- diff(X)
# dV <- diff(V)
R  <- dX + Uprime(X[-length(X)]) * dt
A  <- sort(R, decreasing = TRUE)
m  <- length(A)

# --------------------- Hill estimator
hill_alpha <- function(y_sorted_desc) {
  m <- length(y_sorted_desc)
  logY <- log(y_sorted_desc)
  csum_logY <- cumsum(logY[1:(m-1)])
  logY_kp1  <- logY[2:m] 
  k_vec     <- 1:(m-1) 
  gamma_k <- (csum_logY - k_vec * logY_kp1) / k_vec
  1 / gamma_k
}

runmed1 <- function(x, k = 5) {
  n <- length(x); if (n == 0) return(x)
  k <- max(1, as.integer(k)); if (k %% 2 == 0) k <- k + 1
  half <- k %/% 2
  xm <- x
  for (i in seq_len(n)) {
    lo <- max(1, i - half); hi <- min(n, i + half)
    xm[i] <- median(x[lo:hi], na.rm = TRUE)
  }
  xm
}


choose_k_via_plateau <- function(hill_vals, k_can, alpha_cap = 1.99) {
  m1 <- length(hill_vals)

  k_min <- 20
  k_max <- min(as.integer(m1 / 100), 20000, m1 - 10)

    #precompute for mean and variation
  cs  <- c(0, cumsum(hill_vals))
  cs2 <- c(0, cumsum(hill_vals^2))
  
  # create log-spaced candidate k-grid
  k_grid <- unique(as.integer(round(exp(seq(log(k_min), log(k_max), length.out = k_can)))))
  k_grid <- k_grid[k_grid >= k_min & k_grid <= k_max]
  
  best_k <- NA_integer_
  best_score <- Inf
  
  for (k in k_grid) {
    a <- as.integer(k / 2)
    if (a < 1) a <- 1
    len <- k - a + 1
    sum_seg   <- cs[k + 1]  - cs[a]
    sumsq_seg <- cs2[k + 1] - cs2[a]
    mean_seg  <- sum_seg / len
    var_seg   <- max(sumsq_seg / len - mean_seg * mean_seg, 0)
    sd_seg    <- sqrt(var_seg)
    
    if (!is.finite(sd_seg)) next
    if (hill_vals[k] > alpha_cap) next
    
    if (sd_seg < best_score) {
      best_score <- sd_seg
      best_k     <- k
    }
  }
  
  if (is.na(best_k)) best_k <- k_min
  list(k = best_k, score = best_score, range = range(k_grid))
}

# ----- Alpha estimation
hill_vec_raw <- hill_alpha(A)
hill_vec     <- runmed1(hill_vec_raw, k = 5)
k_can <- 100
sel      <- choose_k_via_plateau(hill_vec, k_can)
k_hat    <- sel$k
alpha_hat<- hill_vec[k_hat]

q_idx <- k_hat + 1
q_val <- if (q_idx <= m) A[q_idx] else A[m]
q_prob<- (q_idx)/(m+1)
cat(sprintf("\nEstimated alpha  = %.4f\n", alpha_hat))
#cat(sprintf("Tail threshold |R| at k: %.4e  (empirical exceed prob = %.4g)\n", q_val, q_prob))
#cat(sprintf("True alpha = %.4f\n", alpha_true))

# ===== sigma2 estimator
estimate_sigma2_tail_robust <- function(alpha_hat, A_sorted_desc, k_hat, dt_vec,
                                        w_half = 20, k_bounds = NULL) {
  m <- length(A_sorted_desc)

  dtrob <- median(dt_vec)
  c_alpha <- gamma(alpha_hat) * sin(pi * alpha_hat / 2) / pi
  
  if (is.null(k_bounds)) {
    k_lo <- max(2, k_hat - w_half)
    k_hi <- min(m-1, k_hat + w_half)
  } else {
    k_lo <- max(2, k_bounds[1])
    k_hi <- min(m-1, k_bounds[2])
  }
  
  k_grid <- k_lo:k_hi
  s2h <- numeric(length(k_grid))
  idx <- 0
  for (kj in k_grid) {
    idx <- idx + 1
    u_j   <- A_sorted_desc[kj]
    p_j   <- (kj - 0.5) / m
    num   <- p_j * u_j^alpha_hat
    den   <- 2 * c_alpha * dtrob
    s2a   <- num / den
    s2h[idx] <- if (is.finite(s2a) && s2a > 0) s2a^(1 / alpha_hat) else NA_real_
  }
  sigma2_hat_med <- median(s2h[is.finite(s2h)], na.rm = TRUE)
  
  list(
    sigma2_hat = sigma2_hat_med,
    details = list(k_hat = k_hat, k_lo = k_lo, k_hi = k_hi,
                   dtrob = dtrob, c_alpha = c_alpha, window_size = length(k_grid))
  )
}

sig2_out <- estimate_sigma2_tail_robust(alpha_hat, A, k_hat, dt)
cat(sprintf("Estimated sigma2  = %.4f\n", sig2_out$sigma2_hat))
#cat(sprintf("  window k in [%d,%d] (size=%d) | c(alpha)=%.5f | median(dt)=%.3e\n",
            # sig2_out$details$k_lo,
            # sig2_out$details$k_hi,
            # sig2_out$details$window_size,
            # sig2_out$details$c_alpha,
            # sig2_out$details$dtrob))
#cat(sprintf("True sigma2 = %.4f\n", sigma2_true))

