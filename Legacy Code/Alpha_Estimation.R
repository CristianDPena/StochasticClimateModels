

# ----- Simulation parameters -----
Terminal   <- 1
n_steps    <- 1000000
alpha_true <- 1.50
sigma2_true<- 1.20
sigma1_true<- 1.00
rho_true   <- 1.60

# ----- U'(x) model
Uprime <- function(x) 16*x^3 - 16*x

# ----- Simulate via Yumia
mod2D <- setModel(
  drift = c("-(16*x^3 - 16*x) - sigma1*(rho^2)*v",
            "-rho^2*v"),
  diffusion = matrix(c("sigma1*rho*sqrt(1+v^2)",
                       "rho*sqrt(1+v^2)"),
                     nrow = 2, byrow = TRUE),
  jump.coeff = matrix(c("sigma2","0"), nrow = 2),
  measure.type = "code",
  measure = list(df = "rstable(z, alpha, 0, 1, 0)"),
  state.variable = c("x","v"),
  solve.variable = c("x","v"),
  xinit = c("x"=0,"v"=0)
)

samp <- setSampling(Terminal = Terminal, n = n_steps)
yu   <- setYuima(model = mod2D, sampling = samp)
sim  <- simulate(yu, true.par = list(
  sigma1 = sigma1_true, rho = rho_true,
  sigma2 = sigma2_true, alpha = alpha_true
))
z      <- yuima::get.zoo.data(sim)
time   <- as.numeric(index(z[[1]]))
X      <- as.numeric(z[[1]])

# ----- Hill estimator
hill_alpha <- function(y_sorted_desc) {
  m <- length(y_sorted_desc)
  if (m < 5) stop("Not enough tail points for Hill.")
  logY <- log(y_sorted_desc)
  csum_logY <- cumsum(logY[1:(m-1)])
  logY_kp1  <- logY[2:m]
  k_vec     <- 1:(m-1)
  gamma_k <- (csum_logY - k_vec * logY_kp1) / k_vec
  1 / gamma_k
}

# ----- Strict extreme-tail k selector
choose_k_via_stability <- function(hill_vals,
                                   k_min_frac = 2e-4,   # 0.02%
                                   k_max_frac = 5e-3,   # 0.5%
                                   winsize    = 40,
                                   curvature_penalty = 0.15,
                                   alpha_cap  = 1.98) {
  m1 <- length(hill_vals)
  k_min <- max(20,  floor(k_min_frac * m1))
  k_max <- min(2000, ceiling(k_max_frac * m1), m1 - 10)
  if (k_max <= k_min + winsize) k_max <- min(m1 - 10, k_min + winsize + 20)
  if (k_min >= k_max - 10) stop("Tail range too narrow for stability selection.")
  
  d2 <- rep(NA_real_, m1)
  d2[2:(m1-1)] <- hill_vals[3:m1] - 2*hill_vals[2:(m1-1)] + hill_vals[1:(m1-2)]
  
  best_k <- NA_integer_; best_sc <- Inf
  for (k0 in seq(k_min, k_max - winsize)) {
    idx <- k0:(k0 + winsize)
    seg <- hill_vals[idx]
    if (!any(is.finite(seg))) next
    if (median(seg, na.rm = TRUE) > alpha_cap) next
    sd_seg <- sd(seg, na.rm = TRUE)
    curv   <- sqrt(mean(d2[idx]^2, na.rm = TRUE))
    score  <- sd_seg + curvature_penalty * curv
    if (is.finite(score) && score < best_sc) { best_sc <- score; best_k <- k0 + floor(winsize/2) }
  }
  if (is.na(best_k)) best_k <- max(20, floor(k_min + winsize/2))
  list(k = best_k, score = best_sc, range = c(k_min, k_max))
}

# ----- Alpha estimation -----
dt <- diff(time)
dX <- diff(X)
R  <- dX + Uprime(X[-length(X)]) * dt
A  <- sort(abs(R[is.finite(R) & R != 0]), decreasing = TRUE)
m  <- length(A)

hill_vec <- hill_alpha(A)
sel      <- choose_k_via_stability(hill_vec)
k_hat    <- sel$k
alpha_hat<- hill_vec[k_hat]

# ----- Diagnostics -----
# df_hill <- data.frame(k = 1:length(hill_vec), alpha_hat = hill_vec)
# p1 <- ggplot(df_hill, aes(k, alpha_hat)) +
#   geom_line() +
#   geom_vline(xintercept = k_hat, linetype = 2) +
#   labs(title = "Hill plot (extreme-tail selection)",
#        subtitle = paste0("k = ", k_hat, ",  alpha = ", sprintf("%.3f", alpha_hat)),
#        x = "k (top order stats)", y = expression(hat(alpha))) +
#   theme_minimal()
# 
# df_ll <- data.frame(
#   logy    = log(A),
#   logtail = log(1 - (1:m)/(m+1))
# )
# p2 <- ggplot(df_ll, aes(logy, logtail)) +
#   geom_point(alpha = 0.35, size = 0.6) +
#   labs(title = "Log–log tail plot", x = "log |R|", y = "log P(|R|>y)") +
#   theme_minimal()
# 
# print(p1); print(p2)

# ----- Console output -----
q_idx <- k_hat + 1L
q_val <- if (q_idx <= m) A[q_idx] else A[m]
q_prob<- (q_idx)/(m+1)

cat(sprintf("\nEstimated alpha (Hill) = %.4f  using k = %d top order stats\n", alpha_hat, k_hat))
cat(sprintf("Tail threshold |R| at k: %.4e  (empirical exceed prob = %.4g)\n", q_val, q_prob))
cat(sprintf("True alpha = %.4f\n", alpha_true))
