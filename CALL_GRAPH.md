# Call Graph

## File-Level Graph

```mermaid
flowchart TD
  C["work/run_overarching_validation_batch.ps1<br/>orchestration only"] --> A["work/run_finite_window_corrected_full_pipeline.R<br/>authoritative entry point"]
  A --> R4["rho_level4_oracle_residual_estimator.R"]
  A --> FW["finite_window_corrected_rho_estimator.R"]
  A --> B["work/run_conditional_cf_endpoint_full_pipeline.R"]
  A --> P["matching finite-window correction profile CSV"]
  B --> E["work/diagnose_exact_cf_endpoint_hmm.R"]
  B --> K["outputs/cache/conditional_cf_endpoint_kernel_production.rds"]
  B -. "when --chunked-likelihood" .-> CH["work/chunked_exact_endpoint_likelihood.R"]
  B --> D["frozen_run_design.csv"]
  E --> M["work/diagnose_markov_additive_endpoint_hmm.R"]
  M --> O["work/occupation_transition_clock_prototype.R"]
  CH --> E
```

The source order is significant:

1. the exact-CF file sets `MARKOV_ENDPOINT_LIBRARY_ONLY=1`;
2. it sources the mixed Markov-additive file;
3. that file sources the occupation prototype;
4. the exact-CF file unsets the environment variable;
5. the endpoint runner optionally overrides the full likelihood symbol with
   the chunked exact-equivalent function.

## Reachable Statistical Functions

```mermaid
flowchart TD
  ENTRY["run_finite_window_corrected_full_pipeline.R top level"]
  INNER["run_conditional_cf_endpoint_full_pipeline.R top level"]
  SIM["simulate_residual_euler"]
  TAIL["estimate_alpha_sigma2_hill"]
  HILL["hill_alpha_path"]
  PLATEAU["choose_hill_plateau"]
  S2["estimate_sigma2_tail_robust"]
  BETA["estimate_exact_cf_endpoint_hmm"]
  OBS["block_cf_log_observations"]
  KERNEL["simulate_q_endpoint_cf_model<br/>(precomputed cache generator)"]
  QBINS["make_q_endpoint_breaks"]
  BMOM["exact_cf_conditional_moments"]
  BEMIT["conditional_exact_cf_emission_table"]
  BAGG["conditional_exact_cf_emission_aggregate_chunked_library"]
  BFWD["forward_from_emission_mixtures_chunked"]
  BOPT["global_profile_optimize"]
  BREF["refine_rho_profile_quadratic"]
  BDIAG["observed_profile_diagnostics"]
  L4["rho_level4_estimate"]
  LLOCAL["rho_level4_local_from_cumulative"]
  LH["rho_level4_h_log / rho_level4_h_delta"]
  CORR["rho_fw_correct_level4"]
  INV["rho_fw_invert_ratio"]
  EVAL["rho_fw_evaluation_profile"]
  FINE["rho_fw_finest_profile"]

  ENTRY --> INNER
  INNER --> SIM
  INNER --> TAIL
  TAIL --> HILL
  TAIL --> PLATEAU
  TAIL --> S2
  INNER --> BETA
  BETA --> OBS
  BETA --> QBINS
  BETA -. "only when no compatible cache" .-> KERNEL
  BETA --> BOPT
  BOPT --> BMOM
  BMOM --> BEMIT
  BEMIT --> BAGG
  BAGG --> BFWD
  BETA --> BREF
  BETA --> BDIAG
  ENTRY --> L4
  L4 --> LLOCAL
  L4 --> LH
  ENTRY --> CORR
  CORR --> INV
  INV --> EVAL
  EVAL --> FINE
```

## Production-Reachable Function Locations

| Function | File and lines |
|---|---|
| `simulate_residual_euler` | `work/occupation_transition_clock_prototype.R:21-44` |
| `symmetric_stable_r` | `work/occupation_transition_clock_prototype.R:8-19` |
| `estimate_alpha_sigma2_hill` | `work/occupation_transition_clock_prototype.R:116-152` |
| `hill_alpha_path` | `work/occupation_transition_clock_prototype.R:46-52` |
| `choose_hill_plateau` | `work/occupation_transition_clock_prototype.R:54-87` |
| `estimate_sigma2_tail_robust` | `work/occupation_transition_clock_prototype.R:89-114` |
| `block_cf_log_observations` | `work/occupation_transition_clock_prototype.R:264-370` |
| `stationary_v_sample` / `stationary_q_sample` | `work/occupation_transition_clock_prototype.R:154-168` |
| `stationary_q_quantile` | `work/occupation_transition_clock_prototype.R:170-184` |
| `log_sum_exp` / `regularized_chol` | `work/occupation_transition_clock_prototype.R:1197-1217` |
| `make_q_endpoint_breaks` | `work/diagnose_markov_additive_endpoint_hmm.R:220-240` |
| `state_index` | `work/diagnose_markov_additive_endpoint_hmm.R:242-245` |
| `emission_terms` | `work/diagnose_markov_additive_endpoint_hmm.R:756-788` |
| `refine_rho_profile_quadratic` | `work/diagnose_markov_additive_endpoint_hmm.R:1155-1194` |
| `observed_profile_diagnostics` | `work/diagnose_markov_additive_endpoint_hmm.R:1466-1498` |
| `simulate_q_endpoint_cf_model` | `work/diagnose_exact_cf_endpoint_hmm.R:36-161` |
| `exact_cf_conditional_moments` | `work/diagnose_exact_cf_endpoint_hmm.R:254-334` |
| `conditional_exact_cf_emission_table` | `work/diagnose_exact_cf_endpoint_hmm.R:336-351` |
| `global_profile_optimize` | `work/diagnose_exact_cf_endpoint_hmm.R:735-769` |
| `estimate_exact_cf_endpoint_hmm` | `work/diagnose_exact_cf_endpoint_hmm.R:922-1210` |
| chunked aggregation/forward functions | `work/chunked_exact_endpoint_likelihood.R:3-110` |
| `rho_level4_estimate` | `rho_level4_oracle_residual_estimator.R:119-352` |
| Level-4 local/noise helpers | `rho_level4_oracle_residual_estimator.R:39-117` |
| `rho_fw_correct_level4` | `finite_window_corrected_rho_estimator.R:405-456` |
| correction evaluation/inversion helpers | `finite_window_corrected_rho_estimator.R:289-403` |

## Unreachable Oracle Branches

The following are loaded as definitions but have no edge from the production
entry point:

- `rho_level4_oracle_diagnostics`
- `block_truth`
- known-endpoint/oracle-Qbar likelihoods
- the Markov-additive file’s main diagnostic program
- the exact-CF file’s standalone simulation program
- the occupation prototype’s generator-GMM, EIV, scalar-HMM, and regime-clock
  estimators

See `PRODUCTION_VS_ORACLE_CODE.md` for code locations.

