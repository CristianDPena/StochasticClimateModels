# StochasticClimateModels
## What This Repository Contains

This repository is a R research implementation and validation
bundle for estimating four parameters from residual increments of a
jump-diffusion model:

```text
alpha, sigma2, beta = sigma1 * rho, rho, sigma1 = beta / rho
```


The `Legacy Code/` folder and
`original_source/external_legacy/RhoSigma1Sigma2Alpha_EstimationVFinal.R`
are historical experiments, not part of the frozen estimator.

## Quick Start

Use R 4.5.2. no non-base R packages

Run the lightweight finite-window correction test:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' `
  'original_source/test_finite_window_corrected_rho_estimator.R'
```

The expected final line is:

```text
All finite-window correction implementation tests passed.
```

For a full frozen simulation path, work from `original_source/` and follow
the exact command in `HOW_TO_RUN.md`. A production-resolution path uses
`n = 10,000,000`, is computationally expensive, and should not be mistaken
for a smoke test.

## Required Numerical Objects

The beta likelihood depends on this frozen cache:

```text
original_source/outputs/cache/conditional_cf_endpoint_kernel_production.rds
```

Rho correction also needs a matching finite-window transfer profile. The
main production cell (`n = 10,000,000`, `T = 50`) uses:

```text
original_source/finite_window_correction_numerical_profiles.csv
```

Additional profiles in `original_source/overarching_profiles/` support only
the five other frozen `(n, T)` validation designs. Do not use a profile for a
different observation design or extrapolate the correction map.

## What Production Does and Does Not Use

Production estimation uses observed residual increments, `dt`, and estimated
upstream parameters. In particular:

- beta uses both `alpha_hat` and `sigma2_hat`;
- rho uses `alpha_hat`, `sigma2_hat`, and `beta_hat`;
- sigma1 uses only the final beta and corrected rho estimates;
- no production estimator receives the realized latent `V_t`;
- generic simulated paths in the kernel and transfer-profile tables are
  numerical approximations to model expectations, not oracle information.

The validation runners contain truth values because they simulate and score
paths. Under the frozen production setting `--tail-input=estimated`, truth is
used only after estimation for error reporting. Selectable oracle-tail modes
exist for diagnostics but are not used by the frozen controller.

## Validated Scope and Main Result

The frozen validation comprises 254 retained simulated paths across 12
parameter regimes, production, infill, long-span, stress, and untouched
panels. At the main production design (`n = 10,000,000`, `T = 50`):

- the beta/rho/sigma1 subsystem was broadly reliable: 92 of 96 development
  paths had all three absolute errors below 0.10;
- all four original parameters met that threshold on 68 of 96 paths;
- the principal full-pipeline limitations were alpha estimation and resulting
  sigma2 error propagation;
- high-clock, strong-jump paths with extreme latent excursions were the main
  downstream exception;
- extending the horizon to `T = 100` improved rho and sigma1 substantially.

Accordingly, the validation supports a restricted research-use claim rather
than an unrestricted guarantee that every four-parameter estimate is within
0.10 at `n = 10,000,000`, `T = 50`.

## Important Implementation FYI

- The frozen beta stage uses a 21-state endpoint approximation, a block
  horizon of `0.05`, three empirical-CF frequencies, 17 rho candidates, and
  beta constrained to `[0.20, 5.00]`.
- The frozen rho stage uses `k = floor(sqrt(n))`, a local ECF floor of
  `k^(-1/2)`, and correction-map inversion only within its supported rho
  domain `[0.4, 3.5]`.
- The final validation uses the memory-bounded chunked likelihood evaluator
  with 16 groups per chunk.
- The frozen source hashes, numerical caches, profiles, and design CSVs are
  the reproducibility authority. One PowerShell controller amendment affects
  logging and process auditing only; the statistical source remains frozen.
