# Overarching Four-Parameter Validation Protocol

## Purpose

This is a frozen evaluation of the existing sequential estimator:

\[
\widehat\alpha
\rightarrow \widehat\sigma_2
\rightarrow \widehat\beta
\rightarrow \widehat\rho_{\rm finite-window}
\rightarrow \widehat\sigma_1=\widehat\beta/\widehat\rho.
\]

No accuracy result generated under this protocol may be used to alter an
estimator formula, threshold, frequency, window, optimization setting,
transfer construction, or parameter ordering.

The only implementation changes made for this validation are:

- deterministic scenario and seed generation;
- observation-design-specific transfer-profile generation and validation;
- batch orchestration with a three-worker ceiling;
- quoted Windows process arguments;
- checkpoint and completed-run resume handling;
- immutable failed-attempt records and explicit retry authorization;
- pathwise diagnostic, runtime, and process-memory logging;
- result aggregation, uncertainty calculations, classification, and figures.

These changes do not alter any statistic or estimate.

## Frozen Scenario Panel

The panel contains twelve regimes, all inside parameter bounds already used
by the project:

```text
alpha  in [1.08, 1.58]
sigma2 in [0.60, 1.40]
sigma1 in [0.70, 1.60]
rho    in [0.90, 2.20]
beta   in [0.8925, 1.92]
```

It includes the required central, low-clock, high-clock, B07 sigma2-hard,
H09 alpha-hard, and extreme-ECF anchors, plus the required marginal and
interaction regimes. Exact values are in `frozen_scenario_design.csv`.

## Frozen Observation Experiments

### A. Production

```text
T = 50
n = 10,000,000
```

Eight development paths and one untouched path are frozen per regime.
Development replications 1 through 3 form the feasibility stage.

### B. Fixed-T Infill

```text
T = 50
n in {1,000,000, 4,000,000, 10,000,000}
```

Three paths per regime are frozen at the lower two resolutions. The
10-million arm reuses production development replications 1 through 3.

### C. Approximately Fixed-Delta Long Span

```text
T = 25,  n = 2,500,000
T = 50,  n = 5,000,000
T = 100, n = 10,000,000
```

Three paths per cell are frozen for six regimes: central, low clock, high
clock, sigma2-hard, alpha-hard, and minimum-alpha strong-jump.

### D. Fixed-n Span Tradeoff

This optional experiment is not scheduled in the initial manifest. It may
be added only after A through C if resource use remains feasible, and must
receive a new frozen seed supplement before execution.

### E. Heavy-Tail Stress

Four additional production-design paths are frozen in five high-risk
regimes. No path is selected in response to new outcomes.

## Replication And Untouched Structure

All seeds are unique and predeclared in `frozen_seed_manifest.csv`.

- production development: eight paths per regime;
- production untouched: one new path per regime;
- infill: three paths per lower-resolution cell;
- long span: three paths per cell;
- stress: four extra paths per high-risk regime.

Untouched production paths have priority 6 and cannot be opened before the
development and stress panels are classified. Any statistical estimator
change invalidates them.

## Numerical Transfer Profiles

The finite-window transfer depends on

\[
q=1/\log(e+n),\qquad
H=\lfloor\sqrt n\rfloor T/n.
\]

Each distinct `(n,T)` cell therefore receives its own transfer profile,
using the unchanged construction:

- rho grid `0.4, 0.5, ..., 3.5`;
- primary 50,000 paired paths, seed `941001`;
- independent 30,000 paired paths, seed `941002`;
- nested 64, 128, and 256 time resolutions;
- first-order Richardson evaluation from 128 and 256;
- chunk size 1,000.

A run aborts if its profile `q` or `H` does not match its observation design.

## Execution Order And Stopping

1. production feasibility;
2. production primary expansion;
3. fixed-T infill;
4. fixed-spacing long span;
5. heavy-tail stress;
6. untouched production validation.

Stop before later priorities if a decisive structural failure, invalid
numerical behavior, or infeasible resource requirement appears. No failure
may trigger estimator repair inside this program.

## Computational Envelope

The host has eight physical cores, sixteen logical processors, and 15.92 GB
of physical memory. A production worker has previously required up to about
1.36 GB, while only about 7 GB was free at freeze time.

The worker limit is therefore three, leaving an explicit operating-system
and application margin. Nine-worker execution is prohibited on this host.

Checkpoints are stored on drive `D:` because the workspace drive had only
about 10 GB free at freeze time. Final results, logs, manifests, and summaries
remain in the workspace.

Completed canonical results are never rerun. Failed attempts are retained
with attempt-numbered logs and status records and are not retried unless the
batch controller is invoked explicitly with `-RetryFailed`.

## Evaluation Rules

Every path, failure, warning, and nonfinite result is retained. No trimming,
winsorization, truth-based selection, fallback estimate, or per-case setting
is permitted.

The required summaries use median, mean, RMSE, q75, q90, q95, maximum,
success probability, Wilson interval, numerical-failure probability,
runtime, and memory.

Downstream and complete four-parameter success are reported separately.
Sigma1 cancellation is explicitly flagged.

## Failure Taxonomy

Errors above `0.10` are assigned, after standard diagnostics, to:

- alpha tail-selection failure;
- sigma2 scale-estimation failure;
- alpha-to-sigma2 propagation;
- beta scale failure;
- corrected-rho finite-window failure;
- local ECF resolution failure;
- extreme latent-excursion failure;
- upstream-to-rho propagation;
- beta/rho cancellation;
- optimizer failure;
- numerical nonconvergence;
- parameter-boundary failure;
- unexplained statistical variation.

Latent `V` may be regenerated only after observable diagnostics are recorded
and only for failure localization.
