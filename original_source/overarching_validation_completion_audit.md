# Overarching Validation Completion Audit

Date: 2026-07-25

## Scope Audit

| Requirement | Evidence | Status |
|---|---|---|
| Frozen estimator | `frozen_method_manifest.md`, `frozen_method_hashes.csv` | Complete |
| Frozen scenarios | 12 rows in `frozen_scenario_design.csv` | Complete |
| Frozen seeds and runs | 254 rows in both frozen manifests | Complete |
| Production benchmark | 108 paths: 96 development plus 12 untouched | Complete |
| Infill ladder | 36 paths at n=1M and 36 at n=4M; production supplies n=10M | Complete |
| Fixed-spacing long-span ladder | 18 paths each at T=25, 50, 100 | Complete |
| Heavy-tail stress | 20 paths | Complete |
| Untouched validation | 12 pre-frozen paths | Complete |
| Retain adverse paths | 254/254 result files and rows; 2 nonfinite retained | Complete |
| Alpha diagnostics | pathwise and stage tables | Complete |
| Sigma2 oracle-alpha decomposition | `sigma2_error_decomposition.csv` | Complete |
| Beta diagnostics | pathwise profile/filter diagnostics | Complete |
| Rho correction diagnostics | numerator, denominator, map, ECF fields | Complete |
| Sigma1 contribution/cancellation diagnostics | pathwise contribution fields | Complete |
| Failure taxonomy | 191 multi-label rows across 99 paths | Complete |
| Runtime/memory | internal R and external process measurements | Complete |
| Latent excursion localization | exact seed replay on 128 production-resolution/stress paths | Complete |
| Required summary tables | all named CSV artifacts present and nonempty | Complete |
| Required figures | 11 original figures plus latent-association figure | Complete |
| Final 12 questions | answered in final report | Complete |
| Fixed-n span tradeoff | optional; not performed | Not required |

## Integrity Audit

- Frozen manifest rows: `254`.
- Unique canonical result files: `254`.
- Aggregated result rows: `254`.
- Valid one-to-one joins: `254`.
- Reconciled attempts: `254`.
- Raw execution events: 254 starts, 253 complete, and one false launcher
  failure reconciled from a valid canonical result.
- Fully finite output rows: `252`.
- Retained nonfinite rows: `2`, both at `n=1M`.
- Duplicate result basenames: `0`.
- Core statistical hash matches: all.
- Transfer-profile hash matches: all.
- Frozen design hash matches: all.
- Controller hash mismatch: one documented orchestration-only amendment.
- Estimator changes after results opened: none.

## Artifact Audit

Required and present:

- `overarching_validation_protocol.md`
- `frozen_method_manifest.md`
- `frozen_scenario_design.csv`
- `frozen_seed_manifest.csv`
- `production_benchmark_results.csv`
- `infill_ladder_results.csv`
- `long_span_ladder_results.csv`
- `heavy_tail_stress_results.csv`
- `parameter_stage_diagnostics.csv`
- `sigma2_error_decomposition.csv`
- `pathwise_full_pipeline_results.csv`
- `regime_summary_results.csv`
- `n_T_scaling_summary.csv`
- `joint_success_summary.csv`
- `failure_classification.csv`
- `runtime_memory_summary.csv`
- `overarching_validation_experiment_log.csv`
- `untouched_validation_results.csv`
- `overarching_four_parameter_validation_final_report.md`

Optional and absent:

- `fixed_n_span_tradeoff_results.csv`, because Experiment D was not
  performed.

Additional diagnostic artifacts:

- `latent_excursion_diagnostics.csv`
- `latent_excursion_failure_summary.csv`
- `latent_excursion_error_correlations.csv`
- `overarching_validation_figures/latent_excursion_failure_association.png`

## Decision Audit

The final decision uses all required evidence:

- medians, means, RMSE, q75, q90, q95, and maxima;
- marginal sub-0.10 probabilities and Wilson intervals;
- joint four-parameter and downstream probabilities;
- regime sensitivity;
- n and T scaling;
- untouched validation;
- stress tests;
- numerical-failure rates;
- error propagation;
- cancellation;
- exact simulation-only latent failure localization;
- runtime and memory.

The evidence contradicts an unrestricted four-parameter production claim
and supports restricted downstream use plus further alpha/sigma2
development. No narrower passing subset was substituted for the requested
full validation program.
