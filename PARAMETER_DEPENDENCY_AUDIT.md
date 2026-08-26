# Parameter Dependency Audit

## Summary Matrix

`Yes` means the quantity is used in the statistical computation of that
estimate. `Evaluation only` means it is referenced only after estimation in
the simulation runner.

| Dependency | alpha_hat | sigma2_hat | beta_hat | rho_hat | sigma1_hat |
|---|---:|---:|---:|---:|---:|
| observed residual increments `R` | Yes | Yes | Yes | Yes | Indirect |
| observed `X_t` directly | No | No | No | No | No |
| known/estimated drift `U'(X_t)` | No, residual expected upstream | No | No | No | No |
| `dt` | Yes, filtering/alignment | Yes | Yes | Yes, scalar first value | Indirect |
| `alpha_hat` | self | Yes | Yes | Yes | Indirect |
| `sigma2_hat` | No | self | Yes | Yes | Indirect |
| `beta_hat` | No | No | self | Yes | Yes |
| `rho_hat` | No | No | No | self | Yes |
| true alpha | No; evaluation only | diagnostic only | No in production | No in production | No |
| true sigma2 | No | No; evaluation only | No in production | No in production | No |
| true beta | No | No | evaluation only | No | No |
| true rho | No | No | evaluation only | evaluation only | No |
| true sigma1 | No | No | evaluation only | evaluation only | evaluation only |
| realized latent `V_t` for observed path | No | No | No | No | No |
| generic unit-model simulations | No | No | Yes, cached candidate-rho kernel | Yes, transfer profile | Indirect |
| oracle-only quantities | No | No in estimate | No in production | No in production | No |

## Alpha Hat

### Statistical inputs

- Residuals and time increments:
  `work/occupation_transition_clock_prototype.R:116-122`
- Production invocation:
  `work/run_conditional_cf_endpoint_full_pipeline.R:164-170`

### Dependencies

- Uses finite, nonzero `abs(R)`.
- Uses `dt` only to keep aligned positive observations and later passes it to
  sigma2; the Hill alpha formula itself depends only on the sorted residual
  magnitudes.
- Does not use sigma2, beta, rho, sigma1, any truth, or latent state.

### Truth branches

`work/run_conditional_cf_endpoint_full_pipeline.R:165-177` contains selectable
tail modes, but the controller fixes `--tail-input=estimated` at
`work/run_overarching_validation_batch.ps1:190`. True alpha is therefore not
used in production alpha estimation.

## Sigma2 Hat

### Statistical inputs

- `alpha_hat`, sorted residual magnitudes, selected `k_hat`, and aligned `dt`:
  `work/occupation_transition_clock_prototype.R:89-114,138-139`

### Dependencies

- Uses `alpha_hat` in the stable-tail constant, power, and inverse power.
- Does not use beta, rho, sigma1, or latent state.
- The outer runner separately recomputes a true-alpha sigma2 estimate at
  `work/run_finite_window_corrected_full_pipeline.R:131-161`; this is an
  evaluation decomposition and never replaces `tail$sigma2_hat`.

## Beta Hat

### Statistical inputs

- Observed residual path, `dt`, `alpha_hat`, and `sigma2_hat`:
  `work/run_conditional_cf_endpoint_full_pipeline.R:226-244`
- The same estimated tail quantities enter the block observation:
  `work/diagnose_exact_cf_endpoint_hmm.R:991-994`
- They also enter every candidate likelihood:
  `work/diagnose_exact_cf_endpoint_hmm.R:1061-1099`

### Candidate-rho dependence

Beta is profiled separately for each candidate rho. Candidate rho determines
the generic endpoint transition/additive-functional kernel, not a truth value:

- candidate model construction:
  `work/diagnose_exact_cf_endpoint_hmm.R:961-972`
- candidate profile loop:
  `work/diagnose_exact_cf_endpoint_hmm.R:1023-1160`

The final downstream `beta_hat` is the beta coordinate interpolated at the
local quadratic maximum of that internal profile:
`work/diagnose_markov_additive_endpoint_hmm.R:1155-1193`.

### Generic versus realized latent paths

The precomputed cache contains generic `Q=1+V^2` paths simulated under each
candidate rho. The generator is
`work/diagnose_exact_cf_endpoint_hmm.R:36-161`. These paths:

- are independent of the realized residual path;
- use frozen seed 5000 and unit model scale;
- approximate parameter-indexed transition/emission expectations;
- are not the true latent `V_t` corresponding to the observed `R`.

### Truth

True beta, rho, and sigma1 appear only in simulation result columns at
`work/run_conditional_cf_endpoint_full_pipeline.R:260-275`. They do not enter
`estimate_exact_cf_endpoint_hmm`.

## Rho Hat

### Statistical inputs

The outer runner passes:

- `sim$R`
- scalar `sim$dt[1]`
- `tail$alpha_hat`
- `tail$sigma2_hat`
- `endpoint$beta_hat`

at `work/run_finite_window_corrected_full_pipeline.R:54-67`.

The stable deconvolution explicitly uses estimated alpha and sigma2 at
`rho_level4_oracle_residual_estimator.R:150-155`. Beta sets:

- frequency `u`, in the outer runner at lines `62-64`;
- denominator normalization, at
  `rho_level4_oracle_residual_estimator.R:162,235-236`.

The correction then uses a precomputed generic unit-clock transfer profile:
`work/run_finite_window_corrected_full_pipeline.R:33-52,68-71`.

### Generic transfer simulations

The transfer profile is generated from paired stationary generic unit-clock
paths in `finite_window_corrected_rho_estimator.R:99-258`. Candidate rho enters
only through `tau=rho^2 H` at line `117`. The profile does not contain the
realized latent path and is independent of alpha, sigma2, and beta after the
normalization `q=(u beta)^2=1/log(e+n)`.

### Truth and oracle code

- `rho_level4_estimate` has no truth or latent-state formal:
  `rho_level4_oracle_residual_estimator.R:119-124`.
- `rho_level4_oracle_diagnostics`, lines `354-501`, is a separately named
  unreachable function in the same file.
- True rho is used only to calculate simulation error after estimation at
  `work/run_finite_window_corrected_full_pipeline.R:86-99,162-177`.

## Sigma1 Hat

Sigma1 is recovered only after corrected rho:

```text
sigma1_hat = beta_hat / rho_hat
```

Code: `finite_window_corrected_rho_estimator.R:418-423`.

It has no direct dependence on `X`, `R`, alpha, sigma2, truth, or latent state
except through the previously estimated beta and rho. There is no bound,
clipping rule, or alternate sigma1 estimator in the authoritative final
coordinate.

## Explicit Tail-Propagation Verification

The production beta stage **does use** both `alpha_hat` and `sigma2_hat`:

- call: `work/run_conditional_cf_endpoint_full_pipeline.R:227-239`
- observation deconvolution:
  `work/diagnose_exact_cf_endpoint_hmm.R:991-994`
- candidate emission moments:
  `work/diagnose_exact_cf_endpoint_hmm.R:1068-1099`

The production rho stage **also uses** both estimates:

- call: `work/run_finite_window_corrected_full_pipeline.R:55-64`
- stable exponent:
  `rho_level4_oracle_residual_estimator.R:39-42,150-153`

## X-to-Residual Dependency

No authoritative frozen file reads `X` or computes `U'(X)`. The validated
statistical contract is an observed residual-increment vector. For an actual
`X` path, drift removal must occur upstream, with the drift function supplied
by the model or estimated separately. The only bundled source line showing
the project’s older quartic-potential formula is explicitly classified
obsolete:

- `external_legacy/RhoSigma1Sigma2Alpha_EstimationVFinal.R:26`
- `external_legacy/RhoSigma1Sigma2Alpha_EstimationVFinal.R:57-58`

