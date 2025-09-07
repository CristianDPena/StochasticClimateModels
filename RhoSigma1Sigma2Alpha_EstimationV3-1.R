rm(list = ls())
library(yuima)
library(zoo)
library(ggplot2)

true_sigma1 <- 0.2
true_sigma2 <- 0.9
true_rho    <- 2.2
true_alpha  <- 1.1

# Set the number of loop iterations
n_loops <- 10
mod2D <- setModel(
  drift         = c(
    "- (16*x^3 - 16*x) - sigma1*(rho^2)*v",  # dX_t drift
    "- rho^2 * v"                            # dV_t drift
  ),
  diffusion     = matrix(
    c("sigma1*rho*sqrt(1+v^2)",
      "rho*sqrt(1+v^2)"),
    nrow = 2, byrow = TRUE
  ),
  jump.coeff    = matrix(c("sigma2","0"), nrow=2),  
  measure.type  = "code",
  measure       = list(df = "rstable(z, alpha, 0, 1, 0)"),
  state.variable = c("x","v"),      
  solve.variable = c("x","v"),
  xinit         = c(x = 0, v = 0)   # start both processes at zero
)

# Sampling
samp <- setSampling(Terminal = 1, n = 1e6)

# store the parameter estimate
sigma_estimates <- numeric(n_loops)  
rho_estimates   <- numeric(n_loops)  
sig2_estimates  <- numeric(n_loops)  
alpha_estimates <- numeric(n_loops)

###1Koutrouvelis Regression
#define a helper function 
#fits log(-log|phi(u)|^2 against log|u| to get alpha and sig2
estimate_sig2_alpha <- function(residuals, dt,
                                J = 10, max_iter = 10, tol = 1e-4) {
  alpha_curr <- 1.5
  for (it in 1:max_iter) {
    # rescale residuals
    z       <- residuals / (dt^(1/alpha_curr))
    u_grid  <- seq(0.5, 1.5, length.out = J)
    #empirical characteristic function 
    phi_hat <- sapply(u_grid, function(u) mean(exp(1i * u * z)))
    phi_sq  <- Mod(phi_hat)^2
    
    good <- which(phi_sq > 0 & phi_sq < 1) #Keep only 0 < |phi|^2 < 1
    if (length(good) < 3) stop("Too few valid u’s; widen your u-range")
    
    y         <- log(-log(phi_sq[good]))
    v         <- log(abs(u_grid[good]))
    fit       <- lm(y ~ v)
    
    # Update alpha
    alpha_new <- coef(fit)[2]
    if (abs(alpha_new - alpha_curr) < tol) break
    alpha_curr <- alpha_new
  }
  intercept <- coef(fit)[1]
  # get sig2 from the intercept
  sig2      <- exp((intercept - log(2)) / alpha_curr)
  list(alpha = alpha_curr, sig2 = sig2)
}

###2Simulation
for (i in 1:n_loops) {
  # Simulate one path of (X, V) using the known “true” parameters
  mod2D_yuima <- setYuima(model = mod2D, sampling = samp)
  res <- simulate(mod2D_yuima,
                  true.par = list(
                    sigma1 = true_sigma1,
                    sigma2 = true_sigma2,
                    rho    = true_rho,
                    alpha  = true_alpha
                  ))
  sim_data <- get.zoo.data(res)
  X_zoo    <- sim_data[[1]]
  V_zoo    <- sim_data[[2]]
  
  #extract  X_t and V_t
  X_obs <- as.numeric(coredata(X_zoo))
  V_obs <- as.numeric(coredata(V_zoo))
  time  <- as.numeric(index(X_zoo))
  
  ###2 Parameter esitmation
  dt     <- diff(time)
  dV     <- diff(V_obs)
  V_left <- V_obs[-length(V_obs)]
  
  # Define the quasi–likelihood contrast function
  Q_rho <- function(rho_cand) {
    var_term <- rho_cand^2 * (1 + V_left^2) * dt
    err      <- dV + rho_cand^2 * V_left * dt
    sum(log(var_term) + (err^2) / var_term, na.rm = TRUE)
  }
  #optimizer
  rho_hat <- optimize(Q_rho, c(0.1, 5))$minimum
  
  #Estimate sig1 via simple regression
  Uprime      <- function(x) 16*x^3 - 16*x
  dX          <- diff(X_obs)
  resid_reg   <- dX + Uprime(X_obs[-length(X_obs)]) * dt
  sigma1_hat  <- sum(resid_reg * dV) / sum(dV^2)
  
  sigma_estimates[i] <- sigma1_hat
  rho_estimates[i]   <- rho_hat
  
  #get sig1*V_t to isolate the jump part
  residuals          <- resid_reg - sigma1_hat * dV
  #run Koutrouvelis
  kou                <- estimate_sig2_alpha(residuals, dt[1])
  sig2_estimates[i]  <- kou$sig2
  alpha_estimates[i] <- kou$alpha
  
  cat(sprintf(
    "Iteration %d sig1 = %.4f rho = %.4f sig2 = %.4f alpha = %.4f\n",
    i, sigma1_hat, rho_hat, kou$sig2, kou$alpha
  ))
}

#3 output
summ <- function(est, true) {
  avg <- mean(est)
  ae  <- mean(abs(est - true))
  mx  <- max(est)
  mn  <- min(est)
  sprintf(
    "Average:       %.7f \nAverage Error: %.9f \nMax:        %.7f \nMin:           %.7f",
    avg, ae, mx, mn
  )
}

cat("\nSummary for sigma1:\n")
cat(summ(sigma_estimates, true_sigma1), "\n\n")

cat("Summary for rho:\n")
cat(summ(rho_estimates, true_rho), "\n\n")

cat("Summary for sigma2:\n")
cat(summ(sig2_estimates, true_sigma2), "\n\n")

cat("Summary for alpha:\n")
cat(summ(alpha_estimates, true_alpha), "\n")
