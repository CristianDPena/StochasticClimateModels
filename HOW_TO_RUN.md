# How To Run

## Environment

- Working directory:
  project root containing `work/`, `outputs/cache/`, and the frozen CSV files
- Validated R executable:
  `C:\Program Files\R\R-4.5.2\bin\Rscript.exe`
- R version:
  `4.5.2`
- R packages:
  no non-base packages are required by the authoritative chain
- Base/recommended namespaces used:
  `stats`, `utils`, `tools`; `parallel` is used only to generate transfer
  profiles
- Yuima:
  not used by the authoritative estimator
- Environment variables:
  none required; the source chain temporarily manages
  `MARKOV_ENDPOINT_LIBRARY_ONLY`
- R option:
  `exact_cf.use_phi_cache` defaults to `TRUE` and was not overridden

All commands below assume PowerShell in the project root.

## 1. Existing Lightweight Smoke Test

This test exercises transfer construction, Richardson evaluation, monotone
inversion, boundary failure behavior, and sigma1 recovery. It does not run the
full alpha/beta pipeline and does not write an output file.

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' `
  'test_finite_window_corrected_rho_estimator.R'
```

Expected final line:

```text
All finite-window correction implementation tests passed.
```

Audit verification: this existing test completed successfully in 3 seconds on
the current host on 2026-07-25. It uses modest memory and writes no files.

## 2. One Production-Resolution Validation Path

Do **not** run this command as part of the source audit. It is the direct
equivalent of one frozen controller invocation for case 10001:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' `
  'work/run_finite_window_corrected_full_pipeline.R' `
  '--design-file=frozen_run_design.csv' `
  '--case=10001' `
  '--tail-input=estimated' `
  '--output-tag=manual_10001' `
  '--chunked-likelihood' `
  '--groups-per-chunk=16' `
  '--kernel-cache=outputs/cache/conditional_cf_endpoint_kernel_production.rds' `
  '--checkpoint-root=D:\CodexOverarchingValidationCheckpoints' `
  '--correction-profile=finite_window_correction_numerical_profiles.csv' `
  '--validation-output=outputs/manual_10001_result.csv'
```

Required inputs:

- `frozen_run_design.csv`
- `outputs/cache/conditional_cf_endpoint_kernel_production.rds`
- `finite_window_correction_numerical_profiles.csv`
- the eight frozen statistical source files

Expected outputs:

- intermediate timestamped endpoint CSV and RDS under `outputs/`;
- final timestamped corrected CSV under `outputs/`;
- `outputs/manual_10001_result.csv`;
- optional per-rho checkpoint RDS files.

Observed validation resource use at `n=10,000,000,T=50`:

- median total runtime about 1,403 seconds;
- maximum reported path runtime about 1,661 seconds;
- median reported peak R memory about 993 MB;
- endpoint beta stage median about 1,332 seconds;
- corrected rho stage median about 6.3 seconds.

The checkpoint path used in the final validation was on `D:`. For a different
machine, `--checkpoint-root=outputs/checkpoints` is computationally
acceptable, but it is not the literal frozen controller path. Checkpoints
change restart behavior only when their exact signature matches.

## 3. User-Supplied Observed X Path

### Frozen-interface limitation

There is no authoritative command-line runner for a user-supplied `X`
vector. The frozen runner always simulates `R`. Creating such a wrapper was
outside this source-recovery task and no wrapper has been added.

The exact residual-level functions can nevertheless be invoked interactively
without altering them. The observed input must be an RDS list:

```text
X          numeric vector of length n+1
time       strictly increasing numeric vector of length n+1
drift_left numeric vector of U'(X_k), length n
```

For the frozen rho implementation, sampling must be equally spaced. A
matching frozen correction profile must exist for the resulting `(n,T)`.
Available cells are listed in `frozen_method_manifest.md:188-195`.

Start an interactive session:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\R.exe' --vanilla
```

Then execute the frozen functions in this order:

```r
source("rho_level4_oracle_residual_estimator.R")
source("finite_window_corrected_rho_estimator.R")
source("work/diagnose_exact_cf_endpoint_hmm.R")

conditional_exact_cf_forward_loglik_full_matrix <-
  conditional_exact_cf_forward_loglik
source("work/chunked_exact_endpoint_likelihood.R")
groups_per_chunk <- 16L
conditional_exact_cf_forward_loglik <- function(
    y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    full_covariance = FALSE, return_filter = FALSE) {
  if (isTRUE(full_covariance)) {
    return(conditional_exact_cf_forward_loglik_full_matrix(
      y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
      full_covariance = TRUE, return_filter = return_filter
    ))
  }
  conditional_exact_cf_forward_loglik_chunked_library(
    y_mat, model, beta, freqs, alpha, sigma2, dt_reference, block_size,
    groups_per_chunk = groups_per_chunk,
    return_filter = return_filter
  )
}

input <- readRDS("observed_path.rds")
X <- as.numeric(input$X)
time <- as.numeric(input$time)
dt <- diff(time)
drift_left <- as.numeric(input$drift_left)

stopifnot(
  length(X) == length(time),
  length(drift_left) == length(dt),
  all(is.finite(X)),
  all(is.finite(dt) & dt > 0),
  max(abs(dt - dt[1L])) <= 1e-12 * max(1, abs(dt[1L]))
)

R <- diff(X) + drift_left * dt
n_steps <- length(R)
terminal <- sum(dt)

tail <- estimate_alpha_sigma2_hill(R, dt)

cache_path <-
  "outputs/cache/conditional_cf_endpoint_kernel_production.rds"
cache <- readRDS(cache_path)
rho_grid <- seq(0.7, 2.3, by = 0.1)
stopifnot(
  isTRUE(all.equal(cache$config$block_horizon, 0.05, tolerance = 1e-10)),
  identical(as.integer(cache$config$n_states), 21L),
  identical(as.integer(cache$config$transition_sim_blocks), 630000L),
  identical(as.integer(cache$config$substeps), 32L),
  identical(as.integer(cache$config$max_path_nodes), 31L),
  identical(as.integer(cache$config$kernel_seed), 5000L),
  isTRUE(all.equal(cache$config$rho_grid, rho_grid, tolerance = 1e-10))
)

endpoint <- estimate_exact_cf_endpoint_hmm(
  R = R,
  dt = dt,
  alpha_hat = tail$alpha_hat,
  sigma2_hat = tail$sigma2_hat,
  block_horizon = 0.05,
  rho_grid = rho_grid,
  freq_mults = c(0.25, 0.50, 0.75),
  n_states = 21L,
  transition_sim_blocks = 630000L,
  substeps_per_block = 32L,
  max_path_nodes_per_pair = 31L,
  kernel_seed = 5000L,
  beta_profile_grid_size = 25L,
  kernel_models = cache$models,
  covariance_model = "conditional_diag",
  checkpoint_dir = NULL,
  return_observations = FALSE,
  verbose = TRUE
)

level4 <- rho_level4_estimate(
  dR = R,
  dt = dt[1L],
  alpha = tail$alpha_hat,
  sigma2 = tail$sigma2_hat,
  beta = endpoint$beta_hat,
  k = floor(sqrt(n_steps)),
  u = 1 / (
    endpoint$beta_hat * sqrt(log(exp(1) + n_steps))
  ),
  chunk_size = 100000L,
  return_observed_cache = FALSE
)

correction_profile_path <-
  "finite_window_correction_numerical_profiles.csv"
profile_rows <- read.csv(
  correction_profile_path,
  stringsAsFactors = FALSE
)
profile <- profile_rows[
  profile_rows$replication == "primary",
  ,
  drop = FALSE
]
expected_q <- 1 / log(exp(1) + n_steps)
expected_h <- floor(sqrt(n_steps)) * terminal / n_steps
stopifnot(
  nrow(profile) > 0,
  max(abs(profile$q - expected_q)) <= 1e-12,
  max(abs(profile$h - expected_h)) <= 1e-12
)

corrected <- rho_fw_correct_level4(level4, profile)

estimates <- list(
  alpha_hat = tail$alpha_hat,
  sigma2_hat = tail$sigma2_hat,
  beta_hat = endpoint$beta_hat,
  rho_hat = corrected$rho_hat,
  sigma1_hat = corrected$sigma1_hat,
  tail_diagnostics = tail,
  beta_profile = endpoint$profile,
  beta_profile_diagnostics = endpoint$profile_diagnostics,
  rho_uncorrected = level4,
  rho_corrected = corrected
)
print(estimates[c(
  "alpha_hat", "sigma2_hat", "beta_hat", "rho_hat", "sigma1_hat"
)])
```

This interactive sequence uses the exact frozen statistical functions and
settings but is **not** itself a frozen, validated entry point. For an
arbitrary observation design, the method cannot be run exactly until a
matching semigroup transfer profile exists.

## Cache Details

The large cache was intentionally not duplicated in the review directory:

```text
Path:
outputs/cache/conditional_cf_endpoint_kernel_production.rds

Type:
R serialized list

Size:
61,305,934 bytes

SHA-256:
3fb889a85e93af3a6c2782f02b06afae05a2ccf975d6ab430f078c94b5af032f
```

It contains `config` and 17 candidate-rho models. Each model stores Q-endpoint
transition probabilities and representative within-block activity paths.
The current component generator and assembler are included under
`original_source/work/`, but the cache itself is the frozen numerical
authority and should not be regenerated for reproduction.
