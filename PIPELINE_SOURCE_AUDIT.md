# Pipeline Source Audit

## Executive Finding

The authoritative 254-path methodology is the residual-level chain launched
by:

`work/run_finite_window_corrected_full_pipeline.R`

All eight frozen statistical files, the endpoint kernel, all transfer
profiles, and all frozen designs still match their frozen SHA-256 values.
The only frozen mismatch is the amended PowerShell controller. The estimator
does not use the realized true latent `V_t`.

The source bundle is complete for reviewing and rerunning a simulated
validation estimate from the existing project. It is not a standalone
observed-`X` application because:

1. the 61.3 MB kernel was intentionally not duplicated;
2. no frozen runner accepts a user-supplied `X` path;
3. correction profiles exist only for six frozen `(n,T)` cells.

## 1. Exact Authoritative Entry Point

`work/run_finite_window_corrected_full_pipeline.R`

Evidence:

- frozen manifest hash list:
  `frozen_method_manifest.md:262-273`
- controller runner selection:
  `work/run_overarching_validation_batch.ps1:14-16`
- exact controller arguments:
  `work/run_overarching_validation_batch.ps1:185-205`
- canonical result architecture:
  `work/run_finite_window_corrected_full_pipeline.R:101-104`

## 2. Scripts Implementing Each Stage

| Stage | Authoritative source |
|---|---|
| residual simulation for validation | `work/occupation_transition_clock_prototype.R:21-44` |
| alpha | `work/occupation_transition_clock_prototype.R:46-87,116-152` |
| sigma2 | `work/occupation_transition_clock_prototype.R:89-114,138-150` |
| beta observations | `work/occupation_transition_clock_prototype.R:264-370` |
| beta endpoint kernel/profile likelihood | `work/diagnose_exact_cf_endpoint_hmm.R:36-351,922-1210` |
| beta state/profile helpers | `work/diagnose_markov_additive_endpoint_hmm.R:220-245,1155-1194,1466-1498` |
| chunked beta likelihood | `work/chunked_exact_endpoint_likelihood.R:3-110` |
| uncorrected rho statistic | `rho_level4_oracle_residual_estimator.R:39-352` |
| finite-window correction | `finite_window_corrected_rho_estimator.R:289-456` |
| final orchestration and sigma1 replacement | `work/run_finite_window_corrected_full_pipeline.R:33-226` |

## 3. Bundle Completeness

For review: **yes**. Every authoritative source and small numerical/config
object is copied byte-for-byte under `original_source/`.

For execution from the same project: **yes**, because the external kernel
exists at its recorded path.

As a standalone portable archive: **no**, because the large kernel is
referenced rather than duplicated, exactly as requested.

For an arbitrary observed `X` path: **no frozen command-line interface
exists**. `HOW_TO_RUN.md` provides a transparent interactive invocation of the
unchanged residual-level functions, not a new wrapper.

## 4. Does Production Estimation Use True V?

No.

- simulator caller sets `keep_v=FALSE`:
  `work/run_conditional_cf_endpoint_full_pipeline.R:155-161`
- tail function accepts only `R,dt`:
  `work/occupation_transition_clock_prototype.R:116-117`
- beta function accepts only `R,dt,alpha_hat,sigma2_hat` plus numerical
  settings: `work/diagnose_exact_cf_endpoint_hmm.R:922-942`
- rho function accepts only `dR,dt,alpha,sigma2,beta` plus numerical settings:
  `rho_level4_oracle_residual_estimator.R:119-124`

The simulator internally evolves V to generate R, but does not return it.
Generic candidate-rho kernel paths and transfer simulations are numerical
expectations, not the realized latent path.

## 5. Stages Using Alpha Hat And Sigma2 Hat

- sigma2 uses alpha_hat:
  `work/occupation_transition_clock_prototype.R:89-114`
- beta uses both:
  `work/run_conditional_cf_endpoint_full_pipeline.R:227-239`;
  `work/diagnose_exact_cf_endpoint_hmm.R:991-994,1061-1099`
- rho uses both:
  `work/run_finite_window_corrected_full_pipeline.R:55-64`;
  `rho_level4_oracle_residual_estimator.R:150-153`
- sigma1 uses them indirectly through beta and rho.

## 6. Generic Simulated/Cached Expectations

### Endpoint cache

- path:
  `outputs/cache/conditional_cf_endpoint_kernel_production.rds`
- size:
  61,305,934 bytes
- SHA-256:
  `3fb889a85e93af3a6c2782f02b06afae05a2ccf975d6ab430f078c94b5af032f`
- object:
  17 rho-indexed Q-endpoint transition/additive-path models
- config:
  rho `0.7:2.3`, H `0.05`, 21 states, 630,000 simulation blocks, 32
  substeps, 31 path nodes, seed 5000
- read:
  `work/run_conditional_cf_endpoint_full_pipeline.R:111-137`
- mathematical generator:
  `work/diagnose_exact_cf_endpoint_hmm.R:36-161`

The current component builder and assembler are included for provenance:

- `work/build_exact_cf_kernel_component.R`
- `work/assemble_exact_cf_kernel_cache.R`

Their filesystem timestamps postdate creation of the authoritative cache, so
they are not asserted to be the literal scripts that wrote that file.
However, the current rho-0.7 component is numerically equal at zero tolerance
to the corresponding cached model. The frozen RDS, not regenerated
components, remains authoritative.

### Transfer profiles

Six small CSVs contain unit-clock finite-window transfer expectations. The
profile generator is `finite_window_corrected_rho_estimator.R:99-287`; the
frozen cell builder is `work/build_overarching_correction_profiles.R`.

## 7. Oracle Branch Reachability

Oracle tail modes exist but the controller sets `estimated`.
`rho_level4_oracle_diagnostics` and the mixed-library oracle HMM routines are
not called. See `PRODUCTION_VS_ORACLE_CODE.md`.

The result files do contain truth-based error and decomposition columns, but
these are computed after the estimates and do not alter them.

## 8. Frozen Integrity Discrepancies

### Statistical and numerical files

No discrepancies:

- all 8/8 core statistical files match;
- endpoint kernel matches;
- all 6/6 correction profiles match;
- all frozen design files match.

### PowerShell controller

| Item | Frozen | Current |
|---|---:|---:|
| bytes | 7,531 | 9,888 |
| SHA-256 | `5d51665c0e0067ffd7d71c69c494bc8bc48368eda5551ef98de9f4a26eb95b90` | `e55675359a02f641a4d0d57277446fa1d906b89daa6f02f6eef19d5f1374ceed` |

Classification: orchestration/logging only.

The documented amendment added:

- recursive process-tree memory measurement;
- rechecking of all frozen statistical/numerical/design hashes;
- `WaitForExit`;
- validation of a one-row canonical result;
- explicit handling of unavailable launcher exit codes.

Evidence: `frozen_orchestration_amendment_001.md`.

The freeze-time 7,531-byte controller is not present in the workspace, there
is no usable Git history, and no second file with its hash was found.
Therefore an exact line-by-line diff cannot be produced without fabricating
the missing side. The current amended controller is copied and the missing
frozen bytes are reported, not silently replaced.

The first launch used the pre-amendment controller and produced a valid result
that was initially misclassified by the launcher. The expansion/final
orchestration used the amended controller. The statistical command and all
estimates were unchanged.

## 9. Duplicated Implementations

Yes. The repository contains:

- a modular endpoint pipeline also labeled "production";
- the endpoint-only internal rho/sigma1 split;
- earlier Markov-additive and Qbar HMMs;
- particle/state-space, latent-V, volatility-ACF, generator-GMM, EIV,
  signed-bridge, codifference, and multiscale alternatives;
- several older full-pipeline and pilot runners;
- embedded duplicate Hill/sigma2 functions.

`ALTERNATIVE_IMPLEMENTATIONS.md` identifies the main ambiguous files and why
they are not authoritative.

## 10. Recommended Review Order

1. `work/run_finite_window_corrected_full_pipeline.R`
2. `work/run_conditional_cf_endpoint_full_pipeline.R`
3. `work/occupation_transition_clock_prototype.R`
4. `work/diagnose_exact_cf_endpoint_hmm.R`
5. `work/chunked_exact_endpoint_likelihood.R`
6. `work/diagnose_markov_additive_endpoint_hmm.R`
7. `rho_level4_oracle_residual_estimator.R`
8. `finite_window_corrected_rho_estimator.R`
9. `finite_window_correction_numerical_profiles.csv` and kernel metadata
10. `frozen_method_manifest.md`

## 11. Reported Formulas Not Implemented

No missing formula was found in the frozen residual-level five-stage
description. The Hill formula, stable-tail sigma2 scale, conditional endpoint
beta likelihood, rolling log-ECF statistic, semigroup transfer inversion, and
sigma1 ratio all have corresponding source.

Two scope qualifications matter:

1. reports call the method observed-data, but the authoritative runner does
   not implement `X -> R`; it starts from a directly simulated residual path;
2. the beta report’s jump attenuation uses a block mean of
   `dt^(1-alpha/2)`, while the model emission moments use median
   `dt_reference` as if spacing is constant. Frozen simulations are equally
   spaced, so these coincide there; irregular sampling is not the validated
   implementation.

## 12. Implemented Steps Underdescribed In Reports

The source additionally implements:

- median-to-MAD-to-SD fallback for beta frequency scaling;
- ECF modulus floor `1e-8`;
- empirical covariance diagonal floor `1e-10`;
- model moment machine-epsilon floors;
- Q-state bin construction from 300,000 stationary draws with seed 881;
- initial endpoint state `Q_0=1`;
- in-memory kernel rebuild when cache compatibility fails;
- exact-signature per-rho checkpoint reuse;
- beta-profile quadratic-refinement rejection rules;
- all-local-valid requirement in the rho statistic;
- correction-domain failure rather than extrapolation;
- intermediate endpoint CSV/RDS output before final rho replacement.

One provenance nuance is also implemented but not emphasized in the reports:
the endpoint checkpoint signature hashes the exact-CF, Markov-core,
conditional-runner, and chunked files
(`work/run_conditional_cf_endpoint_full_pipeline.R:139-146,198-223`), but it
does not directly hash `work/occupation_transition_clock_prototype.R`, the
outer corrected runner, or the two rho files. The frozen batch controller
hashes those files before launch, so the authoritative validation was
protected; a direct user reusing checkpoints outside that controller should
not mistake the checkpoint signature for a complete source-integrity check.

## 13. Hidden Guards, Floors, Truncations, And Boundaries

They exist and are fully enumerated in
`METHODOLOGY_CODE_WALKTHROUGH.md`. The most consequential are:

- Hill cap and fallback selection;
- beta bounds `[0.20,5.00]`;
- ECF and variance floors;
- rho ECF floor `k^-1/2`;
- strict positivity/all-local-valid rho requirements;
- correction profile q/H check;
- no correction-map extrapolation.

No burn-window switch, dependence reranking, base-gain guard, or alternate
fallback estimator is reachable in the authoritative chain.

## 14. Reproducibility

From the existing project: yes for a frozen simulated path, using the
documented command and existing kernel/profile.

From `methodology_pipeline_review/` alone: not fully, because the large kernel
was intentionally not copied.

For arbitrary observed data: only at the residual-function level. A supported
`X` input adapter and profiles for new `(n,T)` cells are absent.

## 15. Obstacles For An Independent Researcher

1. Obtain the exact endpoint kernel or copy it from the documented path.
2. Work at one of the six frozen observation designs or generate and validate
   a new transfer profile.
3. Supply model-consistent drift values to transform `X` into residuals.
4. Use equally spaced observations for the frozen rho function.
5. Avoid confusing the modular "production" endpoint file or endpoint
   internal split with the authoritative final runner.
6. For exact orchestration archaeology, obtain the missing pre-amendment
   7,531-byte controller; it is not needed for the estimates themselves.

## Final Classification

- Production statistical core: complete and hash-verified.
- Frozen numerical objects: present and hash-verified; large kernel referenced
  rather than copied.
- Oracle code: present in mixed files but unreachable under production
  settings.
- Observed-X interface: absent from the authoritative source.
- Statistical source discrepancy: none.
- Orchestration discrepancy: one documented controller amendment.

## Audit Verification Performed

- all 34 copied files were rehashed against their originals: 34/34 exact;
- all 15 copied R files parsed under R 4.5.2;
- all nine required Markdown review documents were present;
- all explicit file-and-line citations in the review documents resolved and
  were within source line bounds;
- all 22 frozen inventory rows were rehashed: 21 matched and the sole mismatch
  was the documented controller;
- the existing `test_finite_window_corrected_rho_estimator.R` smoke test
  passed;
- no production-resolution simulation, validation suite, or cache generation
  was run.
