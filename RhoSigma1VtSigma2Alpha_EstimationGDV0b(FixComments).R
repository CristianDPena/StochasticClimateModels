################################################################################
#  Two-step estimator for sigma1, rho, sigma2, alpha
#  • yuima simulation n_grid = 1e6
#  • Phase-1 BFGS for sigma1 and rho
#  • Phase-2 CF regression with analytic Brownian removal
#  • adaptive coarse residuals (k = 1000, 2000, 4000, …)
#  • robust log transform to avoid underflow
################################################################################
rm(list = ls())
library(yuima)
library(zoo)

# hidden true parameters (for testing only)
true_sig1  <- 0.5
true_sig2  <- 0.1
true_rho   <- 1.25
true_alpha <- 1.50

n_loops   <- 3
n_grid    <- 1e5
need_pts  <- 30
u_max_cap <- 4000

Uprime <- function(x) 16*x^3 - 16*x

# storage vectors
v_s1    <- numeric(n_loops)
v_rho   <- numeric(n_loops)
v_s2    <- numeric(n_loops)
v_alpha <- numeric(n_loops)

for (iter in seq_len(n_loops)) {
  
  # 0. simulate (X,v) path with yuima
  mod <- setModel(
    drift = c("-(16*x^3 - 16*x) - sigma1*(rho^2)*v",
              "-rho^2*v"),
    diffusion = matrix(c("sigma1*rho*sqrt(1+v^2)",
                         "rho*sqrt(1+v^2)"), 2, 1, byrow = TRUE),
    jump.coeff   = matrix(c("sigma2","0"), 2, 1),
    measure.type = "code",
    measure      = list(df = "rstable(z, alpha, 0, 1, 0)"),
    state.variable = c("x","v"),
    solve.variable = c("x","v"),
    xinit = c(x = 0, v = 0)
  )
  samp <- setSampling(Terminal = 1, n = n_grid)
  yu   <- setYuima(model = mod, sampling = samp)
  sim  <- simulate(yu,
                   true.par = list(
                     sigma1 = true_sig1,
                     sigma2 = true_sig2,
                     rho    = true_rho,
                     alpha  = true_alpha
                   ))
  zoo_data <- get.zoo.data(sim)
  X        <- coredata(zoo_data[[1]])
  time     <- index(zoo_data[[1]])
  
  # 1. Phase-1 BFGS: estimate sigma1 and rho on fine grid
  dt     <- diff(time)
  dX     <- diff(X)
  dV_dag <- dX + Uprime(X[-length(X)]) * dt
  V_dag  <- c(0, cumsum(dV_dag))
  
  obj <- function(p) {
    s1 <- p[1]; rh <- p[2]
    Vl <- V_dag[-length(V_dag)]
    var <- rh^2 * (s1^2 + Vl^2) * dt
    err <- dV_dag + rh^2 * Vl * dt
    sum(log(var) + (err^2)/var)
  }
  fit1 <- optim(c(1,1), obj, method = "L-BFGS-B",
                lower = c(0.1,1), upper = c(5,5))
  sigma1_hat <- fit1$par[1]
  rho_hat    <- fit1$par[2]
  v_s1[iter] <- sigma1_hat
  v_rho[iter]<- rho_hat
  
  V_hat  <- V_dag / sigma1_hat
  dV_hat <- diff(V_hat)
  
  # 2. adaptive coarse-step loop
  k <- 1000
  repeat {
    # ensure k yields at least need_pts+1 coarse points
    max_k <- floor(length(dX) / (need_pts + 1))
    if (k > max_k) k <- max_k
    if (k < 1) stop("Data too short for coarse residuals")
    
    idx <- seq(1, length(dX) - k + 1, by = k)
    # compute coarse increments
    dX_c  <- X[idx + k] - X[idx]
    dt_c  <- time[idx + k] - time[idx]
    dV_c  <- V_hat[idx + k] - V_hat[idx]
    resid <- dX_c + Uprime(X[idx]) * dt_c - sigma1_hat * dV_c
    
    # compute Brownian variance on THIS coarse grid (moved inside loop)
    var_B <- (rho_hat * sigma1_hat)^2 * mean(dt_c)
    
    # 3. adaptive u_max loop for CF regression
    u_max <- 25
    repeat {
      u_vals <- seq(0.1, u_max, length.out = 2000)
      good   <- logical(length(u_vals))
      modphi <- numeric(length(u_vals))
      
      for (j in seq_along(u_vals)) {
        u     <- u_vals[j]
        phi_r <- mean(exp(1i * u * resid))
        # analytic de-Gaussianisation
        phi   <- phi_r * exp(0.5 * var_B * u^2)
        m     <- Mod(phi)
        if (m <= 0 || m >= 1) next
        delta <- 1 - m
        # robust log transform
        y     <- if (delta < 1e-8) 2 * delta else -log(m^2)
        if (!is.finite(y) || y <= 0) next
        good[j]   <- TRUE
        modphi[j] <- m
      }
      
      if (sum(good) >= need_pts) break
      if (u_max >= u_max_cap) { u_max <- NA; break }
      u_max <- u_max * 1.5
    }
    
    if (!is.na(u_max)) {
      # 4. Koutrouvelis regression on good CF points
      u_use <- u_vals[good]
      m_use <- modphi[good]
      y_vec <- log(-log(m_use^2))
      v_vec <- log(u_use)
      fit2  <- lm(y_vec ~ v_vec)
      alpha_hat  <- coef(fit2)[2]
      c0         <- coef(fit2)[1]
      sigma2_hat <- exp((c0 - log(2) - log(mean(dt_c))) / alpha_hat)
      break  # exit coarse-step loop
    }
    
    k <- k * 2
    if (k > n_grid/10)
      stop("Even very coarse k gives no usable φ; sigma2 below numeric range")
  }
  
  v_s2[iter]    <- sigma2_hat
  v_alpha[iter] <- alpha_hat
  
  cat(sprintf(
    "iter %d | sigma1=%.4f rho=%.4f || sigma2=%.4f alpha=%.4f (k=%d)\n",
    iter, sigma1_hat, rho_hat, sigma2_hat, alpha_hat, k
  ))
}

# Summary
cat("\nSummary over", n_loops, "runs:\n")
summ <- function(true, est, name)
  cat(sprintf("%-7s true=%.4f mean=%.4f MAE=%.4f\n",
              name, true, mean(est), mean(abs(est - true))))
summ(true_sig1,  v_s1,    "sigma1")
summ(true_rho,   v_rho,   "rho")
summ(true_sig2,  v_s2,    "sigma2")
summ(true_alpha, v_alpha, "alpha")
