# Comprehensive Frozen Four-Parameter Validation

Date: 2026-07-25

## Executive Decision

The complete estimator does **not** justify an unrestricted four-parameter
production claim at `n = 10,000,000`, `T = 50`.

It does justify a **restricted research-use claim**:

- the `beta`/`rho`/`sigma1` subsystem is broadly reliable at production
  resolution;
- the untouched panel confirms that downstream result without cancellation;
- high-clock, strong-jump paths remain vulnerable when the latent process has
  an extreme excursion;
- increasing the terminal horizon to `T = 100` materially improves rho and
  sigma1;
- alpha is the dominant remaining full-pipeline bottleneck;
- sigma2 is the second bottleneck, with most production sigma2 failures caused
  by propagation of alpha error and a smaller residual scale-estimation
  component.

The correct next research priority is the frozen Hill/tail-scale stage, not
another reconstruction or beta/rho split. The corrected downstream method
should be preserved.

## Frozen Estimator

The tested pipeline was:

1. Hill-type estimation of `alpha`;
2. stable-tail scale estimation of `sigma2`;
3. conditional empirical-CF endpoint likelihood estimation of
   `beta = sigma1 * rho`;
4. finite-window semigroup-corrected rolling log-ECF estimation of `rho`;
5. `sigma1_hat = beta_hat / rho_hat`.

No statistical setting changed after the first accuracy result was opened.
The estimator formulas, frequency rule, block rule, window rule, likelihood,
optimizer, transfer maps, parameter ordering, scenarios, and seeds remained
frozen.

Of the 22 files in the freeze-time hash inventory, 21 still match bit for bit.
The sole mismatch is the PowerShell controller, whose logging and
descendant-memory audit were amended before the main expansion. The amendment
is documented in `frozen_orchestration_amendment_001.md`. All eight core
statistical sources, the kernel cache, all transfer profiles, and all frozen
design files match their frozen hashes.

The user subsequently directed that threshold errors must not halt testing.
That execution-policy change is documented in
`validation_continuation_directive_20260724.md`; it did not change the
estimator or promotion criteria.

## Validation Program

| Experiment | Design | Scenarios | Paths |
|---|---:|---:|---:|
| Production characterization | `n=10M, T=50` | 12 | 96 |
| Infill resolution | `n=1M, T=50` | 12 | 36 |
| Infill resolution | `n=4M, T=50` | 12 | 36 |
| Fixed-spacing long span | `n=2.5M, T=25` | 6 | 18 |
| Fixed-spacing long span | `n=5M, T=50` | 6 | 18 |
| Fixed-spacing long span | `n=10M, T=100` | 6 | 18 |
| Heavy-tail stress | `n=10M, T=50` | 5 | 20 |
| Untouched validation | `n=10M, T=50` | 12 | 12 |
| **Total** |  |  | **254** |

Every manifest row has one unique result file and one reconciled completion
record. All 254 paths are retained. There was no trimming, winsorization,
truth-based deletion, or favorable-run selection.

The raw experiment log contains 254 start events, 253 `complete` events, and
one pre-amendment `failed` event for the first production path. That path has
a valid canonical result and a status record reconciled as complete; the raw
event was caused solely by an unavailable Windows launcher exit code. It was
not rerun, omitted, or counted as a statistical failure.

The optional fixed-`n` span-resolution tradeoff was not performed. The
predeclared fixed-spacing long-span experiment, stress panel, and untouched
panel were completed instead.

## Production Accuracy

### Development And Characterization Panel

There were 96 paths, eight in each of the 12 frozen regimes.

| Parameter | Median | Mean | RMSE | q90 | q95 | Maximum | P(error < 0.10) | 95% Wilson CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| alpha | 0.0388 | 0.0604 | 0.0902 | 0.1481 | 0.1874 | 0.4644 | 0.802 | [0.711, 0.869] |
| sigma2 | 0.0361 | 0.0523 | 0.0696 | 0.1289 | 0.1557 | 0.1759 | 0.833 | [0.746, 0.895] |
| beta | 0.0078 | 0.0171 | 0.0425 | 0.0340 | 0.0428 | 0.2884 | 0.979 | [0.927, 0.994] |
| rho | 0.0259 | 0.0416 | 0.0891 | 0.0689 | 0.0869 | 0.5822 | 0.958 | [0.898, 0.984] |
| sigma1 | 0.0256 | 0.0332 | 0.0494 | 0.0624 | 0.0803 | 0.2501 | 0.979 | [0.927, 0.994] |

Joint outcomes:

- all four original parameters below 0.10: `68/96 = 0.708`
  (Wilson interval `[0.611, 0.790]`);
- alpha and sigma2 both below 0.10: `70/96 = 0.729`;
- beta, rho, and sigma1 all below 0.10: `92/96 = 0.958`;
- rho and sigma1 both below 0.10: `92/96 = 0.958`;
- only alpha failed: `10/96 = 0.104`;
- only sigma2 failed: `7/96 = 0.073`;
- multiple parameter stages failed: `9/96 = 0.094`;
- numerical failures: `0/96`.

The downstream subsystem passes its broad-panel q90 criterion. The complete
pipeline fails because alpha and sigma2 both have q90 above 0.10 and success
probabilities below 0.90.

### Untouched Panel

The untouched panel was opened only after the estimator and all 12 seeds were
frozen.

| Parameter | Median | q90 | Maximum | P(error < 0.10) |
|---|---:|---:|---:|---:|
| alpha | 0.0381 | 0.0988 | 0.1043 | 0.917 |
| sigma2 | 0.0562 | 0.1518 | 0.1537 | 0.833 |
| beta | 0.0079 | 0.0291 | 0.0398 | 1.000 |
| rho | 0.0236 | 0.0720 | 0.0859 | 1.000 |
| sigma1 | 0.0140 | 0.0451 | 0.0515 | 1.000 |

Joint outcomes:

- all four below 0.10: `9/12 = 0.750`;
- alpha and sigma2 both below 0.10: `9/12 = 0.750`;
- downstream triplet below 0.10: `12/12 = 1.000`;
- rho and sigma1 below 0.10: `12/12 = 1.000`;
- sigma1 cancellation flags: `0/12`;
- numerical failures: `0/12`.

The three complete-pipeline misses were:

- `R04`: sigma2 error `0.1537`;
- `R10`: alpha error `0.1043`;
- `R12`: sigma2 error `0.1502`.

Thus the untouched result independently confirms the downstream method and
independently rejects an unrestricted full-pipeline claim.

### All Production-Resolution And Stress Paths

Pooling the 96 characterization paths, 12 untouched paths, and 20
predeclared stress paths gives 128 `n=10M, T=50` paths:

| Parameter | Median | q90 | Maximum | P(error < 0.10) |
|---|---:|---:|---:|---:|
| alpha | 0.0375 | 0.1450 | 0.4644 | 0.812 |
| sigma2 | 0.0369 | 0.1305 | 0.2166 | 0.820 |
| beta | 0.0077 | 0.0350 | 0.2884 | 0.984 |
| rho | 0.0252 | 0.0682 | 0.5822 | 0.953 |
| sigma1 | 0.0222 | 0.0614 | 0.2501 | 0.984 |

All-four success was `88/128 = 0.688` with Wilson interval
`[0.603, 0.761]`. Downstream-triplet success was `122/128 = 0.953` with
Wilson interval `[0.902, 0.978]`.

## Regime Sensitivity

The table gives q90 errors and joint success in the 96-path production
characterization panel.

| Regime | alpha q90 | sigma2 q90 | beta q90 | rho q90 | sigma1 q90 | All-four P | Downstream P |
|---|---:|---:|---:|---:|---:|---:|---:|
| R01 central | 0.141 | 0.130 | 0.013 | 0.057 | 0.048 | 0.750 | 1.000 |
| R02 low clock, weak jump | 0.412 | 0.159 | 0.012 | 0.039 | 0.065 | 0.750 | 1.000 |
| R03 high clock, strong jump | 0.193 | 0.160 | 0.282 | 0.573 | 0.246 | 0.625 | 0.625 |
| R04 sigma2-hard | 0.099 | 0.156 | 0.011 | 0.057 | 0.079 | 0.625 | 1.000 |
| R05 alpha-hard | 0.135 | 0.051 | 0.008 | 0.032 | 0.043 | 0.875 | 1.000 |
| R06 low alpha, high clock | 0.199 | 0.113 | 0.009 | 0.091 | 0.039 | 0.625 | 0.875 |
| R07 high alpha, high sigma2 | 0.043 | 0.050 | 0.039 | 0.066 | 0.085 | 1.000 | 1.000 |
| R08 high beta, high rho | 0.089 | 0.045 | 0.022 | 0.066 | 0.040 | 1.000 | 1.000 |
| R09 high beta, low rho | 0.140 | 0.052 | 0.058 | 0.055 | 0.094 | 0.500 | 1.000 |
| R10 interior high beta | 0.154 | 0.108 | 0.022 | 0.053 | 0.040 | 0.500 | 1.000 |
| R11 low beta, high rho | 0.208 | 0.150 | 0.030 | 0.070 | 0.040 | 0.750 | 1.000 |
| R12 minimum alpha, strong jump | 0.102 | 0.172 | 0.012 | 0.083 | 0.050 | 0.500 | 1.000 |

The main downstream exception is `R03`. `R06` has one marginal rho miss.
The tail-stage failures do not define a simple monotone parameter boundary:
for example, `R07` and `R08` pass every path, while alpha is unstable in
`R02`, `R09`, `R10`, and `R11`. This points to tail-selector/path
interaction, not merely one inadmissible parameter value.

Regime-level probabilities use only eight characterization paths and have
wide uncertainty. They support localization, not precise 0.90 reliability
claims for each individual regime.

## Infill Scaling At Fixed T

The table compares the same 12 scenarios using independent frozen seed
blocks. It is an empirical resolution ladder, not a pathwise coupled
convergence experiment.

| Parameter | n=1M q90 / pass | n=4M q90 / pass | n=10M q90 / pass |
|---|---:|---:|---:|
| alpha | 0.297 / 0.611 | 0.239 / 0.639 | 0.148 / 0.802 |
| sigma2 | 0.224 / 0.639 | 0.152 / 0.778 | 0.129 / 0.833 |
| beta | 0.197 / 0.806 | 0.061 / 0.944 | 0.034 / 0.979 |
| rho | 0.134 / 0.639 | 0.093 / 0.917 | 0.069 / 0.958 |
| sigma1 | 0.162 / 0.750 | 0.088 / 0.917 | 0.062 / 0.979 |

Joint all-four success rises from `0.389` at `n=1M`, to a higher but still
tail-limited level at `n=4M`, and to `0.708` at production resolution.

Two `n=1M` paths produced nonfinite rho and sigma1 because the local spot
variance became invalid. These are the only two numerical failures in the
254-path program. No `n>=4M` path had a numerical failure.

The q90 local-ECF floor fraction declines from `0.00117` at `n=1M`, to
`0.000487` at `n=4M`, to `0.000317` at `n=10M`. This supports higher
frequency for numerical resolution, although the probability of any floor
activation is not monotone.

Increasing n has the clearest effect on beta, rho, and sigma1. Alpha improves
more slowly, and sigma2 inherits both that slow improvement and its own
scale error.

## Long-Span Scaling At Approximately Fixed Delta

| Parameter | T=25 q90 / pass | T=50 q90 / pass | T=100 q90 / pass |
|---|---:|---:|---:|
| alpha | 0.148 / 0.722 | 0.148 / 0.778 | 0.125 / 0.778 |
| sigma2 | 0.144 / 0.833 | 0.102 / 0.889 | 0.128 / 0.833 |
| beta | 0.040 / 0.944 | 0.031 / 1.000 | 0.026 / 1.000 |
| rho | 0.103 / 0.889 | 0.096 / 0.889 | 0.067 / 1.000 |
| sigma1 | 0.079 / 0.944 | 0.073 / 1.000 | 0.059 / 1.000 |

At `T=100`, all 18 beta, rho, and sigma1 errors are below 0.10. The
high-clock/strong-jump `R03` rho errors fall from two misses at `T=50` to a
maximum of `0.0892` at `T=100`.

All-four success is `0.611`, `0.611`, and `0.667` at `T=25`, `50`, and `100`.
Long span therefore helps rho and sigma1 more clearly than it helps alpha or
sigma2. The tail stages remain dominated by heavy-tail sampling variability.

Because each horizon cell contains only 18 independent paths, small
nonmonotonic changes, such as sigma2 q90 between `T=50` and `T=100`, should
not be interpreted as structural reversal.

## Heavy-Tail Stress Panel

Across 20 predeclared stress paths:

| Parameter | Median | q90 | Maximum | P(error < 0.10) |
|---|---:|---:|---:|---:|
| alpha | 0.0347 | 0.2342 | 0.3943 | 0.800 |
| sigma2 | 0.0325 | 0.1541 | 0.2166 | 0.750 |
| beta | 0.0058 | 0.0592 | 0.0655 | 1.000 |
| rho | 0.0160 | 0.0998 | 0.1429 | 0.900 |
| sigma1 | 0.0222 | 0.0572 | 0.0679 | 1.000 |

All-four success is `11/20 = 0.550`; tail-pair success is `13/20 = 0.650`;
downstream-triplet success is `18/20 = 0.900`.

The two downstream stress misses are rho errors in `R03`. The stress panel
again localizes the broad failure to alpha/sigma2 and the rare downstream
failure to high-clock/strong-jump paths.

## Stage Diagnostics

### Alpha And Sigma2

At production characterization resolution:

- alpha fails in `19/96` paths;
- sigma2 fails in `16/96` paths;
- alpha and sigma2 errors have correlation `0.521`;
- `13/16` sigma2 failures would pass if the true alpha were supplied;
- sigma2 q90 is `0.1289` with estimated alpha and `0.0889` with true alpha;
- median absolute alpha-to-sigma2 propagation is `0.0356`;
- q90 absolute alpha-to-sigma2 propagation is `0.1183`;
- median absolute residual sigma2-stage error is `0.0271`;
- q90 residual sigma2-stage error is `0.0889`.

This is not a single isolated sigma2 failure. It is a recurring but
moderate-scale problem across several regimes. Most of it is upstream alpha
propagation, but residual sigma2 sensitivity remains visible in `R04` and
`R12`, including the untouched and stress panels.

There were no small-sample Hill flags in the full study, only one
tail-candidate boundary flag, and no production boundary flag. The failures
therefore are statistical tail-selection/scale errors rather than an obvious
optimizer or coding defect.

### Beta

Beta is the most stable stage at production resolution:

- q90 error `0.0340`;
- success probability `0.979`;
- no optimizer convergence failures;
- no production profile-boundary failure that generated a beta miss.

Its two large production errors occur in `R03` paths with extreme latent
excursions. This is rare-path finite-horizon variability, not ordinary beta
performance.

### Corrected Rho

Rho has q90 `0.0689` and success probability `0.958` in the 96-path
production panel. However, the correction does not erase all high-clock
path variability.

Across all finite validation paths, Spearman correlation between rho error
and the effective window `rho^2 H` is `0.512`. Rho q90 rises from `0.0425`
in the lowest effective-window quartile to `0.1825` in the highest. This is
consistent with residual difficulty when substantial latent clock time
elapses inside a local ECF window.

The `T=100` result shows that longer span supplies the missing dynamic
information: all 18 long-span rho estimates pass, including every `R03`
path.

### Sigma1 And Cancellation

In the production characterization panel:

- sigma1 q90 is `0.0624`;
- success probability is `0.979`;
- only `2/96` sigma1 successes are flagged as beta/rho cancellation;
- the untouched panel has `0/12` cancellations;
- median absolute first-order beta contribution is `0.0055`;
- median absolute first-order rho contribution is `0.0180`;
- median absolute nonlinear remainder is only `0.00043`.

Thus production sigma1 accuracy is mostly genuine downstream coherence, not
accidental cancellation.

The theoretical result remains conditional:

```text
beta_hat -> beta and rho_hat -> rho > 0
imply sigma1_hat = beta_hat / rho_hat -> sigma1.
```

This is a continuous-mapping argument, not an independently proved
sigma1-first theorem. A dummy value `r0` would produce
`beta_hat/r0 -> sigma1*rho/r0`, which is inconsistent unless `r0` itself
converges to rho. Feeding that value back into a second rho estimate is
circular. The current rho stage is already a one-parameter inversion
conditional on beta, so the proposed dummy-rho split does not simplify the
identification problem.

## Exact Latent-Excursion Audit

The estimator never used latent `V`. After all observable diagnostics were
recorded, a simulation-only diagnostic regenerated the exact Euler latent
path from each frozen seed for all 128 `n=10M, T=50`
production/stress/untouched paths. An independent small-path check showed
bit-for-bit equality with the original simulator.

Predeclared excursion thresholds gave:

| Definition | Frequency | All-four failure if extreme | All-four failure otherwise | Downstream failure if extreme | Downstream failure otherwise |
|---|---:|---:|---:|---:|---:|
| max `|V| > 10` | 23/128 = 0.180 | 0.565 | 0.257 | 0.217 | 0.010 |
| max `|V| > 20` | 3/128 = 0.023 | 0.667 | 0.304 | 0.333 | 0.040 |
| above-stationary occupancy at `|V| > 10` | 15/128 = 0.117 | 0.600 | 0.274 | 0.200 | 0.027 |
| top-decile mean latent activity | 13/128 = 0.102 | 0.462 | 0.296 | 0.077 | 0.043 |

Maximum `|V|` has Spearman correlations `0.337`, `0.172`, `0.321`, `0.159`,
and `0.016` with alpha, sigma2, beta, rho, and sigma1 errors respectively.

Five of the six production-resolution/stress downstream failures have
`max |V| > 10`; the sixth has `max |V| = 9.20`. All four downstream failures
in the 96-path production characterization panel have `max |V| > 10`.

The interpretation is specific:

- extreme excursions substantially increase beta/rho failure risk;
- they explain the rare large downstream outliers, especially in `R03`;
- they do not explain every alpha or sigma2 miss;
- they do not make sigma1 generally unstable because beta and rho errors
  often remain coherent, but the two most extreme `R03` paths defeat both.

This is an association diagnostic, not a causal theorem or an observable
real-data gate.

## Numerical Reliability And Cost

- valid manifest/result joins: `254/254`;
- finite complete parameter vectors: `252/254`;
- nonfinite vectors: two `n=1M` resolution paths;
- raw log events: 253 complete plus one documented false launcher-failure,
  all 254 reconciled from canonical outputs;
- production, stress, long-span, and untouched numerical failures: zero;
- optimizer nonconvergence: zero;
- non-successful rho correction status: two `n=1M` spot-variance failures;
- total measured serial CPU time: `99.93` hours;
- median production path runtime: `1403` seconds;
- maximum production path runtime: `1661` seconds;
- median production memory: `993 MB`;
- maximum observed worker memory: `1.60 GB`.

The conditional-CF endpoint beta stage dominates runtime. Median stage times
at production resolution are approximately:

- simulation: `62.6` seconds;
- tail estimation: `2.9` seconds;
- endpoint beta stage: `1332.0` seconds;
- corrected rho stage: `6.3` seconds.

Runtime is almost flat from `n=1M` to `n=10M` at fixed `T=50` because the
endpoint likelihood uses a nearly fixed number of blocks. It grows roughly
with the number of blocks/terminal span in the fixed-spacing ladder:
`702`, `1300`, and `2598` median seconds for `T=25`, `50`, and `100`.

## Multi-Label Failure Taxonomy

The 99 paths with at least one flagged condition generate 191 diagnostic
labels:

| Category | Labels |
|---|---:|
| alpha tail-selection failure | 64 |
| alpha-to-sigma2 propagation | 45 |
| sigma2 scale-estimation failure | 17 |
| beta/rho cancellation | 16 |
| upstream-to-rho propagation | 14 |
| beta scale failure | 12 |
| corrected-rho finite-window failure | 10 |
| local ECF resolution failure | 10 |
| numerical nonconvergence | 2 |
| unexplained statistical variation | 1 |

These counts are deliberately multi-label and should not be added as if they
were mutually exclusive path counts. No optimizer, parameter-boundary, or
model-code failure was detected.

The single unexplained label is ordinary finite-sample variation after the
standard observable diagnostics; it is not evidence of an unidentified
systematic implementation error.

## Answers To The Required Questions

### 1. Reliability at n=10M, T=50

In the 96-path characterization panel, alpha/sigma2/beta/rho/sigma1 pass
probabilities are `0.802/0.833/0.979/0.958/0.979`. Beta, rho, and sigma1 have
q90 below 0.10. Alpha and sigma2 do not.

### 2. Joint probability all four errors are below 0.10

It is `0.708` in characterization, `0.750` in untouched validation, and
`0.688` across all 128 production-resolution/stress paths.

### 3. Dominant bottleneck

Alpha is now the dominant stage. Its production q90 is `0.148` and its pass
probability is `0.802`. Sigma2 is second; much of its error is propagated
from alpha.

### 4. Is sigma2 failure systematic or rare?

It is neither a single rare path nor a universal collapse. It recurs in
`16/96` characterization paths and in the sigma2-hard/low-alpha regimes.
`13/16` production failures would pass with true alpha, so upstream
propagation is the main mechanism. Residual sigma2-stage sensitivity remains
in `R04` and `R12`.

### 5. Improvement with n

Downstream improvement is strong and monotone in the aggregate ladder.
Alpha and sigma2 improve much more slowly. At `n=1M`, numerical failures are
possible; none occurs at `n>=4M`.

### 6. Improvement with T

At approximately fixed observation spacing, beta, rho, and sigma1 improve
with T. Alpha and sigma2 remain variable. All 18 downstream estimates pass at
`T=100`.

### 7. Does rho benefit from infill, long span, or both?

Both. Infill resolves local ECF/activity estimation: rho q90 falls from
`0.134` at `n=1M` to `0.069` at `n=10M`, fixed `T=50`. Long span improves
clock information: rho q90 falls to `0.067` at `T=100`, and all 18 paths pass.
For high-clock paths, long span is the more decisive final improvement.

### 8. Does sigma1 remain accurate without cancellation?

Yes at production resolution. Only `2/96` characterization successes are
flagged as cancellation, and none of the 12 untouched successes is.
Consistency is nevertheless conditional on beta and rho consistency.

### 9. Hardest parameter combinations

`R03` high-clock/strong-jump is the hardest downstream regime. `R04` and
`R12` are hardest for sigma2. Alpha is unstable in several otherwise
different regimes (`R02`, `R09`, `R10`, `R11`), so its difficulty is not
captured by one simple parameter boundary.

### 10. Frequency and role of extreme latent excursions

At production resolution, `17.97%` of audited paths reach `|V|>10` and
`2.34%` reach `|V|>20`. Downstream failure is `21.7%` among the former versus
`1.0%` otherwise. Extreme excursions are the principal explanation for rare
large downstream errors, but not for most tail-stage misses.

### 11. Supported operating region

The supported downstream region is:

- `n >= 10M`, `T >= 50`;
- model-correct observations;
- ordinary latent paths;
- rho inside the frozen transfer domain;
- no local-ECF spot-variance failure.

For high-clock/strong-jump work, `T=100` is better supported than `T=50`.
There is not yet a parameter-only operating region that guarantees full
alpha/sigma2 accuracy.

### 12. Final use decision

- **Unrestricted research use:** no.
- **Restricted use:** yes, especially for the downstream subsystem.
- **Further development:** yes, first alpha tail selection and then residual
  sigma2 scaling.
- **Longer data requirement:** yes for robust high-clock rho; `T=100` is
  strongly supported.
- **Reject the current downstream method:** no.
- **Reject the unrestricted full-pipeline production claim:** yes.

## Publication-Level Interpretation

The evidence supports a coherent methodological story:

- stable-tail statistics identify alpha and sigma2 but remain the
  finite-sample bottleneck;
- the endpoint empirical-CF likelihood identifies beta accurately;
- the finite-window semigroup correction recovers rho from dynamic
  information;
- sigma1 follows by a transparent ratio and is accurate when beta and rho are
  accurate;
- rare heavy-tail latent excursions create finite-horizon outliers that no
  stationary correction can eliminate pathwise.

What is not yet proved is equally important:

- there is no complete joint asymptotic theorem with estimated alpha,
  sigma2, and beta;
- there is no formal heavy-tail failure-probability theorem;
- there is no independent sigma1-first consistency theorem;
- simulation correctness does not test real-data model misspecification.

The publishable claim should therefore be a statistically coherent
downstream estimator with explicit finite-sample operating restrictions and
an honestly identified upstream tail-estimation limitation.

## Artifact Index

Primary data and summaries:

- `pathwise_full_pipeline_results.csv`
- `production_benchmark_results.csv`
- `untouched_validation_results.csv`
- `infill_ladder_results.csv`
- `long_span_ladder_results.csv`
- `heavy_tail_stress_results.csv`
- `parameter_stage_diagnostics.csv`
- `sigma2_error_decomposition.csv`
- `regime_summary_results.csv`
- `n_T_scaling_summary.csv`
- `joint_success_summary.csv`
- `failure_classification.csv`
- `runtime_memory_summary.csv`
- `latent_excursion_diagnostics.csv`
- `latent_excursion_failure_summary.csv`
- `latent_excursion_error_correlations.csv`

Design and integrity:

- `overarching_validation_protocol.md`
- `frozen_method_manifest.md`
- `frozen_method_hashes.csv`
- `frozen_scenario_design.csv`
- `frozen_seed_manifest.csv`
- `frozen_run_design.csv`
- `frozen_orchestration_amendment_001.md`
- `validation_continuation_directive_20260724.md`
- `overarching_validation_experiment_log.csv`

Figures are in `overarching_validation_figures/`, including the added
simulation-only `latent_excursion_failure_association.png`.
