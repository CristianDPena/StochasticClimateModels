# Pipeline Execution Order

## Scope

The authoritative frozen validation entry point is:

`work/run_finite_window_corrected_full_pipeline.R`

The 254-path controller invokes it from the project root. The exact
statistical arguments are constructed in
`work/run_overarching_validation_batch.ps1:185-205`:

```text
Rscript work/run_finite_window_corrected_full_pipeline.R
  --design-file=<absolute frozen_run_design.csv>
  --case=<frozen case number>
  --tail-input=estimated
  --output-tag=ov_<case number>
  --chunked-likelihood
  --groups-per-chunk=16
  --kernel-cache=outputs/cache/conditional_cf_endpoint_kernel_production.rds
  --checkpoint-root=D:\CodexOverarchingValidationCheckpoints
  --correction-profile=<matching profile from frozen_seed_manifest.csv>
  --validation-output=<canonical result CSV>
```

For `n < 10,000,000`, the controller also adds
`--validation-resolution-experiment`. The final production and untouched
paths have `n = 10,000,000`, so that flag is absent.

## Important Input Boundary

The authoritative runner does **not** accept an observed `X` vector. It
simulates residual increments `R` directly at
`work/run_conditional_cf_endpoint_full_pipeline.R:155-162`. Therefore, the
frozen executable chain begins statistically at `(R, dt)`, not at `X`.

The model-level preprocessing expected for observed data is

```text
R_k = X_{k+1} - X_k + U'(X_k) dt_k.
```

That calculation appears in the older, unused OneDrive script at
`external_legacy/RhoSigma1Sigma2Alpha_EstimationVFinal.R:26,57-58`, but no
frozen production function exposes it. This is an interface gap, not a
missing dependency inside the validated residual-level estimator.

## Chronological Trace

### 1. Parse outer-runner arguments

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `8-27`
- Inputs: command-line arguments
- Outputs:
  - `correction_profile_path`
  - `validation_output_path`
- Disk reads: none
- Next: source the statistical pipeline.

### 2. Load the full statistical source chain

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `29-31`
- Direct sources:
  - `rho_level4_oracle_residual_estimator.R`
  - `finite_window_corrected_rho_estimator.R`
  - `work/run_conditional_cf_endpoint_full_pipeline.R`
- Transitive sources:
  - `work/diagnose_exact_cf_endpoint_hmm.R`
  - `work/diagnose_markov_additive_endpoint_hmm.R`
  - `work/occupation_transition_clock_prototype.R`
  - conditionally, `work/chunked_exact_endpoint_likelihood.R`
- Environment behavior:
  - `work/diagnose_exact_cf_endpoint_hmm.R:12-14` temporarily sets
    `MARKOV_ENDPOINT_LIBRARY_ONLY=1` while sourcing the mixed
    Markov-additive file, preventing its diagnostic main program from running.

Sourcing `work/run_conditional_cf_endpoint_full_pipeline.R` executes that
runner immediately before control returns to the outer script.

### 3. Resolve the frozen simulation design

- File: `work/run_conditional_cf_endpoint_full_pipeline.R`
- Lines: `5-64`
- Function: top-level `option_value`
- Inputs:
  - `--design-file=frozen_run_design.csv`
  - `--case=<case>`
- Disk read: `frozen_run_design.csv` at lines `40-55`
- Outputs:
  - `cfg` with `seed`, true simulation parameters, and case
  - `n_steps`
  - `terminal`
  - `panel_name`
  - `run_id`
- Boundary rule:
  - lines `26-28` and `62-64` reject non-quick runs below ten million unless
    `--validation-resolution-experiment` is present.
- Next: configure beta-stage numerics.

### 4. Configure and load the endpoint likelihood

- File: `work/run_conditional_cf_endpoint_full_pipeline.R`
- Lines: `66-146`
- Inputs/settings:
  - rho grid `0.7:2.3` by `0.1`
  - block horizon `0.05`
  - 21 endpoint states
  - 630,000 transition simulation blocks
  - 32 substeps
  - 31 path nodes per endpoint pair
  - kernel seed `5000`
  - 25 beta profile points
  - covariance model `conditional_diag`
  - chunked evaluator enabled, 16 groups per chunk
- Conditional source:
  - `work/chunked_exact_endpoint_likelihood.R:94`
  - lines `92-109` replace the in-memory full-matrix likelihood function with
    the exact-equivalent chunked implementation.
- Disk read:
  - `outputs/cache/conditional_cf_endpoint_kernel_production.rds`,
    lines `111-137`
- Cache checks:
  - block horizon, state count, transition count, substeps, node count, seed,
    and every requested rho must match.
  - If the cache is missing or incompatible, the code silently leaves
    `kernel_models=NULL`; the endpoint function then rebuilds candidate
    kernels in memory. The frozen controller supplied a compatible cache.
- Outputs:
  - `kernel_models`
  - `kernel_cache_used`
  - source MD5 values
- Next: simulate the observed residual path in validation.

### 5. Produce the validation path

- File: `work/run_conditional_cf_endpoint_full_pipeline.R`
- Lines: `155-162`
- Function: `simulate_residual_euler`
- Definition:
  - `work/occupation_transition_clock_prototype.R:21-44`
- Inputs:
  - `n_steps`, `terminal`
  - simulation-only `cfg$alpha`, `cfg$sigma2`, `cfg$sigma1`, `cfg$rho`
  - `cfg$seed`
  - `keep_v=FALSE`
- Internal simulator:
  - generates latent `v` only to construct `R`
  - returns `sim$R` and `sim$dt`
  - does not return the realized latent path
- Outputs: `sim$R`, `sim$dt`
- Next: tail estimation.

For real data, this stage must be replaced by an externally supplied
residual vector and matching positive time increments. The frozen source has
no command-line path for doing that.

### 6. Estimate alpha and sigma2

- File: `work/run_conditional_cf_endpoint_full_pipeline.R`
- Lines: `164-178`
- Production branch: line `169`
- Function: `estimate_alpha_sigma2_hill`
- Definition: `work/occupation_transition_clock_prototype.R:116-152`
- Inputs:
  - `sim$R`
  - `sim$dt`
- Internal calls:
  - `hill_alpha_path`, lines `46-52`
  - `choose_hill_plateau`, lines `54-87`
  - `estimate_sigma2_tail_robust`, lines `89-114`
- Outputs:
  - `tail$alpha_hat`
  - `tail$sigma2_hat`
  - `tail$k_hat`
  - Hill and sigma2-window diagnostics
- Production setting:
  - controller sets `--tail-input=estimated`
- Next: beta endpoint likelihood.

### 7. Estimate beta by the conditional empirical-CF endpoint likelihood

- File: `work/run_conditional_cf_endpoint_full_pipeline.R`
- Lines: `226-245`
- Function: `estimate_exact_cf_endpoint_hmm`
- Definition: `work/diagnose_exact_cf_endpoint_hmm.R:922-1210`
- Inputs:
  - `R=sim$R`, `dt=sim$dt`
  - `alpha_hat=tail$alpha_hat`
  - `sigma2_hat=tail$sigma2_hat`
  - frozen beta and kernel numerical settings
  - prebuilt `kernel_models`
- Observation construction:
  - `block_cf_log_observations`
  - `work/occupation_transition_clock_prototype.R:264-370`
- Candidate transition/additive-functional model:
  - cached objects originally produced by
    `simulate_q_endpoint_cf_model`
  - definition: `work/diagnose_exact_cf_endpoint_hmm.R:36-161`
- Conditional likelihood:
  - theoretical moments:
    `work/diagnose_exact_cf_endpoint_hmm.R:254-334`
  - chunked emission aggregation and forward recursion:
    `work/chunked_exact_endpoint_likelihood.R:3-110`
- Profile/optimization:
  - for each candidate rho, profile beta over `[0.20,5.00]`
  - 25-point global log-beta grid plus local `optimize`
  - code: `work/diagnose_exact_cf_endpoint_hmm.R:1023-1160`
  - local quadratic rho-profile refinement:
    `work/diagnose_markov_additive_endpoint_hmm.R:1155-1194`
- Outputs used later:
  - `endpoint$beta_hat`
  - beta-stage profile and diagnostics
- Outputs not used as the final split:
  - `endpoint$rho_hat`
  - `endpoint$sigma1_hat`
- Disk I/O:
  - optional exact-signature rho checkpoints at
    `work/diagnose_exact_cf_endpoint_hmm.R:1010-1153`
  - intermediate endpoint CSV/RDS written at
    `work/run_conditional_cf_endpoint_full_pipeline.R:346-359`
- Next: return to the outer runner and load the finite-window profile.

### 8. Load and validate the finite-window correction profile

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `33-52`
- Disk read: `correction_profile_path`
- Selection: rows with `replication == "primary"`
- Required design match:
  - `q = 1/log(e+n_steps)`
  - `H = floor(sqrt(n_steps))*terminal/n_steps`
  - absolute tolerance `1e-12`
- Production `n=10,000,000,T=50` object:
  - `finite_window_correction_numerical_profiles.csv`
- Output: `correction_profile`
- Next: rolling log-ECF rho statistic.

### 9. Compute the uncorrected rolling log-ECF statistic

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `54-67`
- Function: `rho_level4_estimate`
- Definition: `rho_level4_oracle_residual_estimator.R:119-352`
- Inputs:
  - `dR=sim$R`
  - scalar `dt=sim$dt[1]`
  - `alpha=tail$alpha_hat`
  - `sigma2=tail$sigma2_hat`
  - `beta=endpoint$beta_hat`
  - `k=floor(sqrt(n_steps))`
  - `u=1/(beta_hat*sqrt(log(e+n_steps)))`
  - chunk size `100000`
- Outputs:
  - numerator `level4$i_log_hat`
  - denominator `level4$denominator_hat`
  - uncorrected `level4$rho_hat`
  - uncorrected `level4$sigma1_hat`
  - floor/failure diagnostics
- Hidden-state behavior:
  - no latent state or truth argument
  - `return_observed_cache=FALSE`
- Next: semigroup correction.

### 10. Apply the finite-window semigroup correction

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `68-71`
- Function: `rho_fw_correct_level4`
- Definition: `finite_window_corrected_rho_estimator.R:405-456`
- Inputs:
  - `level4`
  - `correction_profile`
- Internal calls:
  - `rho_fw_invert_ratio`, lines `358-403`
  - `rho_fw_evaluation_profile`, lines `313-356`
  - `rho_fw_finest_profile`, lines `289-311`
- Calculation:
  - observed ratio is `i_log_hat / denominator_hat`
  - transfer is first-order Richardson evaluated from the finest two
    numerical resolutions
  - the strictly increasing map `rho^2 * transfer(rho)` is linearly inverted
- Outputs:
  - `corrected$rho_hat`
  - correction status/boundary
  - transfer and ratio diagnostics
- Next: recover sigma1.

### 11. Recover sigma1

- File: `finite_window_corrected_rho_estimator.R`
- Lines: `418-423`
- Exact calculation:

```text
sigma1_hat = beta_hat / rho_hat
```

- Bounds/postprocessing: none
- Failure behavior:
  - returns `NA` if corrected rho is nonfinite or nonpositive
- Next: overwrite the beta-stage internal split.

### 12. Assemble final estimates and simulation-only diagnostics

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `74-185`
- Final estimate replacement:
  - preserves internal split in diagnostic columns, lines `75-82`
  - replaces final `rho_hat` and `sigma1_hat`, lines `84-85`
- Simulation-only truth uses:
  - absolute errors and success flags, lines `86-99`
  - true-alpha sigma2 decomposition, lines `131-161`
  - first-order beta/rho contribution and cancellation diagnostics,
    lines `162-177`
- These truth calculations occur after all five estimates have been formed and
  do not feed back into them.

### 13. Write final outputs

- File: `work/run_finite_window_corrected_full_pipeline.R`
- Lines: `187-226`
- Disk writes:
  - timestamped `outputs/finite_window_corrected_full_pipeline_*.csv`
  - optional canonical `--validation-output` CSV
- Final estimate columns:
  - `alpha_hat`
  - `sigma2_hat`
  - `beta_hat`
  - `rho_hat`
  - `sigma1_hat`
- Additional outputs:
  - stage diagnostics, timings, memory, profile support, floors, statuses,
    simulation errors, and source identifiers.

## Required Non-Code Objects

| Object | Required role |
|---|---|
| `outputs/cache/conditional_cf_endpoint_kernel_production.rds` | Exact frozen 17-rho endpoint transition/additive-functional kernel. |
| One matching correction profile CSV | Frozen map from observed rho ratio to corrected rho for the exact `(n,T)` cell. |
| `frozen_run_design.csv` | Required only by the authoritative simulated-validation invocation. |
| `frozen_seed_manifest.csv` | Required by the batch controller, not by one direct runner call. |
| Checkpoint directory | Optional computational output/restart state, not a precomputed statistical input. |

