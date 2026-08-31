# Alternative And Ambiguous Implementations

## Authoritative Implementation

The only implementation used for the 254 canonical final estimates is:

`work/run_finite_window_corrected_full_pipeline.R`

Its frozen eight-source chain and hashes are listed in
`frozen_method_manifest.md:262-273`. The result rows identify the architecture
as:

```text
hill_tail__endpoint_beta__finite_window_rho__ratio_sigma1
```

and carry MD5 identifiers matching the current frozen source files.

## Ambiguous "Production" Entry Point

| Path | Purpose | Status relative to final validation |
|---|---|---|
| `work/conditional_cf_endpoint_pipeline.R` | A later modular function API for endpoint-CF estimation and numerical-accuracy controls. Its first comment calls it a production entry point. | **Not authoritative.** It is not sourced by the frozen runner and is absent from the frozen eight-source hash list. |
| `work/run_accuracy_controlled_full_pipeline_case.R` | Runs the modular accuracy-controlled endpoint pipeline. | Numerical convergence/audit runner, not the 254-path final estimator. |
| `work/run_conditional_cf_endpoint_full_pipeline.R` | Executing validation runner for tail plus endpoint likelihood. | Authoritative transitive dependency, but not final by itself because its internal rho/sigma1 split is overwritten. |
| `work/run_finite_window_corrected_full_pipeline.R` | Adds Level-4 finite-window rho and ratio sigma1. | **Authoritative final entry point.** |

The similarly named modular file is the most likely source of accidental
review confusion.

## Earlier Full-Pipeline And Pilot Runners

| Path | Purpose/difference | Classification |
|---|---|---|
| `work/run_markov_endpoint_full_pipeline.R` | Earlier Markov-additive endpoint ECF-HMM full runner; reports the endpoint likelihood’s own rho/sigma1 split. | Older, replaced for final rho/sigma1. |
| `run_corrected_full_pipeline_validation.R` | Pre-freeze held-out corrected-rho validation using beta estimates loaded from prior result files. | Validation-only; not a self-contained production estimate. |
| `run_corrected_rho_complete_pilot.R` | Pilot combining saved production beta outputs with corrected rho, including oracle/full-estimated comparisons. | Experimental/oracle pilot. |
| `run_corrected_rho_pilot_subset.R` | Small corrected-rho pilot on selected saved cases. | Experimental/oracle pilot. |
| `run_upstream_error_propagation.R` | Tail-input ablations through corrected rho. | Oracle diagnostic. |
| `run_rho_level4_oracle_residual.R` | Direct residual simulation plus detailed oracle variance diagnostics. | Oracle diagnostic; not production. |

## Alpha And Sigma2 Duplicates

The authoritative definitions are only:

- `work/occupation_transition_clock_prototype.R:46-152`

Other files contain independent or older copies:

| Path | Major difference | Status |
|---|---|---|
| `work/state_space_likelihood_y_experiment.R` | Much larger tail framework with alternate selectors and tail-likelihood options. | Retired experimental state-space pipeline. |
| `work/volatility_acf_split_experiment.R` | Embedded Hill/sigma2 copy paired with volatility-ACF rho methods. | Experimental, not frozen. |
| `work/robust_yuima_sweep.R` | Embedded Yuima simulation and earlier tail function. | Simulation experiment. |
| `work/threshold_rho_contrast_sweep.R` | Embedded tail estimation for threshold/rho contrast tests. | Experimental. |
| `work/v_recovery_iteration_harness.R` | Embedded tail estimate feeding latent-V reconstruction. | Rejected earlier direction. |
| `external_legacy/RhoSigma1Sigma2Alpha_EstimationVFinal.R` | Original Hill/sigma2 and multiple V-reconstruction/splitting variants. | Obsolete, not sourced by frozen validation. |

## Beta Implementations

| Path/family | Identification idea | Status |
|---|---|---|
| `work/diagnose_exact_cf_endpoint_hmm.R` | Conditional finite-frequency empirical-CF endpoint likelihood with additive paths. | **Authoritative beta implementation.** |
| `work/diagnose_markov_additive_endpoint_hmm.R` | Earlier linear-Qbar/Markov-additive endpoint likelihood and multiple oracle modes. | Production helper subset remains; its standalone estimator is older. |
| `work/occupation_transition_clock_prototype.R:418-513` | Scalar/vector activity occupation estimators for beta. | Rejected prototype functions; loaded but unreachable. |
| `work/occupation_transition_clock_prototype.R:1219-1477` | Qbar-state direct ECF-HMM. | Negative-result prototype; unreachable. |
| `work/state_space_likelihood_y_experiment.R` | Scaled-state particle/state-space likelihood and profile variants. | Earlier research direction, not final. |
| `work/diagnose_direct_ecf_observation_hmm.R` | Direct empirical-CF HMM using block-average activity as state. | Negative-result diagnostic. |
| `work/signed_bridge_complex_ecf_hmm.R` and `work/signed_bridge_all_path_reference.R` | Signed endpoint diffusion-bridge complex-ECF likelihood. | Later alternative research branch; not frozen final method. |
| `work/codifference_split_estimator.R` | Codifference moment split. | Experimental/failed alternative. |
| `work/multiscale_cf_increment_estimator.R` | Multiscale marginal-CF variance structure. | Experimental alternative. |

## Rho And Sigma1 Implementations

| Path/family | Major difference | Status |
|---|---|---|
| `rho_level4_oracle_residual_estimator.R:119-352` plus `finite_window_corrected_rho_estimator.R:289-456` | Rolling log-ECF numerator/denominator plus semigroup transfer inversion. | **Authoritative rho/sigma1 implementation.** |
| `work/occupation_transition_clock_prototype.R:645-926` | Point generator-GMM and vector-ECF generator-GMM. | Rejected prototypes; not called. |
| `work/occupation_transition_clock_prototype.R:1479-1782` | Block-semigroup, EIV, transition-regime, and activity-clock estimators. | Rejected prototypes; not called. |
| `work/volatility_acf_split_experiment.R` / `work/volatility_acf_split_estimator.R` | Volatility autocovariance/SMM clock split. | Earlier experimental alternative. |
| `work/ecf_variogram_estimator_experiment.R` | ECF variogram split. | Experimental. |
| `work/rv_state_space_split_estimator.R` | Realized-volatility state-space split. | Experimental/failed. |
| `observed_signed_moment_estimator.R` | Signed transition moment for rho. | Diagnostic alternative, not final. |
| `work/signed_bridge_all_path_reference.R` | Two-step signed-bridge rho conditional on beta. | Later experimental branch, not final. |
| `external_legacy/RhoSigma1Sigma2Alpha_EstimationVFinal.R` and V-recovery files | Reconstruct latent V, then estimate sigma1/rho. | Earlier discarded methodology. |

## Finite-Window Correction Implementations

The statistical implementation is unique:

- `finite_window_corrected_rho_estimator.R`

The following are generators, tests, or analyses rather than alternative
estimators:

- `run_finite_window_correction_numerics.R`
- `work/build_overarching_correction_profiles.R`
- `run_finite_window_infill_profiles.R`
- `analyze_level4_finite_window_correction.R`
- `analyze_finite_window_infill_robustness.R`
- `audit_finite_window_corrected_methodology.R`
- `test_finite_window_corrected_rho_estimator.R`

## Mixed-Library Ambiguity

`work/occupation_transition_clock_prototype.R` and
`work/diagnose_markov_additive_endpoint_hmm.R` contain many discarded
estimators alongside helpers that the final method still needs. Their
presence in the source chain does not mean those alternatives execute.
`CALL_GRAPH.md` lists the reachable subset.

No alternative implementation above was used to replace a final estimate in
the frozen 254-path results.

