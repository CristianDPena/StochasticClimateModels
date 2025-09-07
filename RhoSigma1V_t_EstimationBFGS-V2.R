rm(list = ls())
library(yuima)
library(zoo)
library(ggplot2)

true_sigma1 <- 0.5
true_rho    <- 1.5

# Set the number of loop iterations
n_loops <- 1
mod2D <- setModel(
  drift = c(
    "-(16*x^3 - 16*x) - sigma1*(rho^2)*v",  # dX_t drift
    "-rho^2*v"                              # dV_t drift
  ),
  diffusion = matrix(
    c("sigma1*rho*sqrt(1+v^2)", "rho*sqrt(1+v^2)"),
    nrow = 2, ncol = 1, byrow = TRUE
  ),
  jump.coeff = matrix(c("sigma2", "0"), nrow = 2, ncol = 1),
  measure.type = "code",
  measure = list(df = "rstable(z, alpha, 0, 1, 0)"),
  state.variable = c("x", "v"),
  solve.variable = c("x", "v"),
  xinit = c("x" = 0, "v" = 0)
)

# Sampling
samp <- setSampling(Terminal = 1, n = 100000)

# store the parameter estimate
sigma_estimates <- numeric(n_loops)
rho_estimates <- numeric(n_loops)

start.time <- proc.time()

###1Simulation
for (i in 1:n_loops) {
  mod2D_yuima <- setYuima(model = mod2D, sampling = samp)
  res <- simulate(mod2D_yuima,
                  true.par = list(sigma1 = true_sigma1,
                                  sigma2 = 0,
                                  rho = true_rho,
                                  alpha = 1.5))
  sim_data <- get.zoo.data(res)
  
  #extract  X_t
  X_zoo <- sim_data[[1]]
  time_index <- index(X_zoo)
  X_obs <- as.numeric(coredata(X_zoo))
  time <- as.numeric(time_index)
  
###2 Parameter esitmation
  # Reconstruct the proxy for V_t using U'(x)
  Uprime <- function(x) { 16 * x^3 - 16 * x }
  dt <- diff(time)         # Time increments
  dX <- diff(X_obs)        # X increments
  
  # Reconstruct the proxy increments dV_dagger
  dV_dagger <- dX + Uprime(X_obs[-length(X_obs)]) * dt
  V_dagger <- c(0, cumsum(dV_dagger))   # Proxy for V_t (note: V_dagger = sigma1 * V_t)
  
  # Define the quasi–likelihood contrast function
  Q <- function(sigma_candidate, rho_candidate) {
    V_left <- V_dagger[1:(length(V_dagger) - 1)]  # left endpoints
    dV_current <- dV_dagger                       # observed increments
    var_term <- rho_candidate^2 * (sigma_candidate^2 + V_left^2) * dt
    err <- dV_current + rho_candidate^2 * V_left * dt
    sum( log(var_term) + (err^2) / var_term, na.rm = TRUE )
  }
  
  # objective function for optimizer
  objective <- function(params) {
    sigma_candidate <- params[1]
    rho_candidate <- params[2]
    Q(sigma_candidate, rho_candidate)
  }
  
  # initial guess
  init_params <- c(1, 1)
  
  #optimizer
  optim_result <- optim(
    par = init_params,
    fn = objective,
    method = "L-BFGS-B",
    lower = c(0.1, 1),
    upper = c(5, 5)
  )
  
  sigma_estimates[i] <- optim_result$par[1]
  rho_estimates[i] <- optim_result$par[2]
  
  cat("Iteration", i, "\n")
}

#3 output
sigma_avg <- mean(sigma_estimates)
rho_avg <- mean(rho_estimates)
sigma_avg_error <- mean(abs(sigma_estimates - true_sigma1))
rho_avg_error <- mean(abs(rho_estimates - true_rho))
sigma_max <- max(sigma_estimates)
sigma_min <- min(sigma_estimates)
rho_max <- max(rho_estimates)
rho_min <- min(rho_estimates)

cat("\nSummary for sigma1:\n")
cat("Average:      ", sigma_avg, "\n")
cat("Average Error:", sigma_avg_error, "\n")
cat("Max:          ", sigma_max, "\n")
cat("Min:          ", sigma_min, "\n")

cat("\nSummary for rho:\n")
cat("Average:      ", rho_avg, "\n")
cat("Average Error:", rho_avg_error, "\n")
cat("Max:          ", rho_max, "\n")
cat("Min:          ", rho_min, "\n")

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("\nTotal Time Passed:", time.taken, "\n")
