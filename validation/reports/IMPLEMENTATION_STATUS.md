# Validation package status

The repository contains the complete reduced study design and execution/reporting
framework. The frozen estimator files are unchanged. There are 216 principal fits
and 40 separately labeled diagnostic evaluations, using 168 independent trajectory
IDs. Seven named batches run through a two-worker, resumable wall-clock scheduler.
This is a multi-session study, not a claim that all experiments finish in 12 hours.

## Completed software verification

- R source, production-kernel and profile integrity checks.
- Exact design counts, unique job IDs and shared per-trajectory simulation seeds.
- Native coupled Euler simulation agrees with an R reference on identical Brownian
  increments to 1e-12; coarser observations are exact increment sums.
- Chunked and original forward likelihoods agree to 1e-8 on the same summaries.
- Both tail benchmarks tested against a simulated stable law, including scale units.
- Failure denominators verified separately from finite-error moments.
- Scheduler rejects more than two workers or a budget above eleven hours; a
  near-deadline test starts no expensive fit and cleans its lock.
- End-to-end small-sample ordinary and oracle smoke jobs completed. The ordinary
  job retained a tail-domain failure; the oracle job reached all downstream stages
  with finite output. These are computational checks, not scientific evidence.
- Summary CSVs, paired contrasts, vector PDF figures and a LaTeX recovery table
  generated from saved records; final calibration figure layouts visually checked.

Test environment: R 4.6.1, evir 1.7-4, Apple arm64 C compiler, Python 3.

## Scientific findings

No full-study robustness or comparative-superiority finding is available yet.
See the run-specific `summary/validation_report.md` for actual completion counts,
failures and pending jobs. A two-path, full-resolution production calibration has completed. Both ordinary
fits returned finite estimates; fixed-estimator times were 278.9 and 286.7 seconds,
with the complete two-worker session taking 296.2 seconds. See
[the calibration report](calibration/README.md). These two baseline paths do not
justify robustness or superiority claims across the parameter grid. The remaining
214 principal and 40 diagnostic evaluations are pending.

The strongest planned manuscript evidence is ordinary per-regime recovery and
paired information-growth comparisons. Tail benchmark comparisons and oracle
contrasts support specific stage-level conclusions only. Numerical sensitivity
and materially adverse outcomes should be reported even when they weaken claims.
