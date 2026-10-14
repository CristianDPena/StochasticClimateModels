# Fixed-method validation on a 16 GB MacBook Air

**Completed study:** [results and interpretation](reports/completed-study/README.md)
cover all 256 planned evaluations. Regenerate the compact report with
`python3 validation/export_results.py`; this also audits the retained raw run.


This package evaluates the existing estimator; it does not redesign or tune it.
The authoritative implementation is `original_source/`, pinned to upstream commit
`e6eb6da6d4db1325c9ca086dd5df2e09abceb915`. Every estimator source and production
kernel is checked against `config/frozen_sources.csv` before execution. Checksums
are file-integrity checks, not a security guarantee.

## Revised study (replaces the proposed 750-fit study)

| Regime | alpha | sigma2 | sigma1 | rho | beta |
|---|---:|---:|---:|---:|---:|
| B | 1.30 | 1.00 | 1.00 | 1.50 | 1.50 |
| Aminus | 1.10 | 1.00 | 1.00 | 1.50 | 1.50 |
| Aplus | 1.60 | 1.00 | 1.00 | 1.50 | 1.50 |
| Jminus | 1.30 | 0.50 | 1.00 | 1.50 | 1.50 |
| Cminus | 1.30 | 1.00 | 1.875 | 0.80 | 1.50 |
| Sminus | 1.30 | 1.00 | 0.50 | 1.50 | 0.75 |
| H | 1.55 | 1.40 | 0.70 | 2.20 | 1.54 |

All seven regimes have 24 independent replications at T=50, n=10 million:
168 ordinary fits, split across four recovery batches (six seeds per regime in
each batch). B, Jminus, and H are fixed diagnostic choices, not selected by error.

For B and H, the first 12 trajectories also supply n=2.5 million at T=50 and
n=20 million at T=100 with spacing 0.000005: 48 additional fits in two information
batches. All comparisons use the same paths, aggregation, or prefixes. Baseline
summaries in the information figure use the same first 12 trajectory IDs, not all
24 reference replications. There are only two levels per information axis; do not
fit or claim convergence rates.

The diagnostic batch contains 40 evaluations:
- True alpha/sigma2, then true alpha/sigma2/beta, on four fixed datasets each in
  B, Jminus, H (24 evaluations). The latter skips the likelihood.
- Double within-block kernel substeps (32 to 64), holding 21 states, 31 path
  representatives, clock grid, and 30,000 transitions/state fixed, on two datasets
  each in B and H (4). This is a targeted refinement check, not exhaustive
  convergence of every likelihood quadrature dimension. Kernel paths are not
  Brownian-coupled across step counts; shifts include Monte Carlo variation.
- Latent simulator substeps 4 to 8 on the same four datasets (4).
- Refined transfer (512 steps, 100,000 pairs) and independent standard-resolution
  transfer (seed 941002) on those four datasets (8); reuse the original beta fit.

The old three-level information ladders, broad state/node/grid expansions, and
redundant regimes are omitted. The fixed estimator's own statistical grids,
frequencies, thresholds, bounds, and tuning rules are **not** reduced. There is no
true-alpha-only ablation in this bounded study. Additional investigations require
a separately declared extension, not opportunistic retuning.

## Expected resources and limitations

Two full-resolution baseline fits were run concurrently on this MacBook Air.
The fixed estimator took 278.9 and 286.7 seconds; the complete two-job session took
296.2 seconds. A mid-run resident-memory snapshot was about 1.4 GB per worker.
These are short-run calibration measurements, not peak-memory guarantees or
measurements of all regimes and refinement settings.

At that throughput a recovery batch projects to about 1.7 hours, and an information
batch to about 1.4 hours if cost scales with horizon. Budget 2–4 hours per such
session for regime variation, generation, and sustained thermal load. Diagnostics
also include new cache construction and are less predictable. A rough full-study
projection is 10–16 hours over the seven sessions, not a guaranteed complete
12-hour study. The scheduler's bounded sessions remain the reliable constraint.

The scheduler runs at most two R processes, limits numerical-library threads,
and stops after a default 10.5-hour calendar-time budget (including cache preparation
and system sleep). On macOS, keep the lid open and consider `caffeinate -i` while
running the study; a suspended computer cannot execute the deadline handler until
it wakes. Elapsed fit times can include suspension and should not be interpreted
as CPU cost. It does
not launch a new fit in the final 45 minutes and interrupts unfinished work at the
deadline. Completed clock-profile rows are checkpointed. Repeat the SAME batch
command to resume; completed job IDs are never rerun. The seven named batches
are scientifically non-overlapping but a slow batch may require another session.
Summary generation and final figure rendering each have a separate 15-minute ceiling. A time-bounded session does
not guarantee that every planned fit in the batch finishes.

Raw storage is approximately 16 GB of residual doubles, plus small additional
refinement trajectories, kernels, checkpoints and outputs (budget 16–25 GB).
Full latent and stable-component arrays are not retained. Seeds, simulator source,
raw observed increments, per-horizon latent extrema and checksums are retained.
Two simultaneous likelihood fits are intended to fit within 16 GB; measure actual
peak resident memory during the first batch. Use one worker if system pressure is
high. Do not run other batches concurrently.

## Setup and commands

Run commands from the repository root. Requires R (tested version recorded in the
verification report), a C compiler usable by `R CMD SHLIB`, Python 3, and `evir`
1.7-4. No global R library is modified. On macOS install Command Line Tools if
needed. The core estimator remains Base R.

```sh
Rscript validation/setup.R
Rscript validation/tests/test_validation.R
Rscript validation/run.R --command=smoke --out=validation/results/smoke-v1
Rscript validation/run.R --command=summarize --mode=smoke --out=validation/results/smoke-v1

# First production timing check: two real, full-resolution fits, retained in study.
python3 validation/schedule.py --batch=recovery1 --max-jobs=2

# Then run these one session at a time, repeating a command if it stops pending.
python3 validation/schedule.py --batch=recovery1
python3 validation/schedule.py --batch=recovery2
python3 validation/schedule.py --batch=recovery3
python3 validation/schedule.py --batch=recovery4
python3 validation/schedule.py --batch=information1
python3 validation/schedule.py --batch=information2
python3 validation/schedule.py --batch=diagnostics
Rscript validation/run.R --command=summarize
Rscript validation/figures.R
```

Use `--out=/path/on/local/disk` on every command to avoid intensive cloud-sync I/O.
All batches use one shared output directory with disjoint job IDs. `--hours=8`
changes the session budget, not the scientific design. Do not change the generated
CSV to fit a time limit. Framework changes require a new output directory; source
and design checks reject mixed revisions. An interrupted scheduler may leave a
`scheduler.lock`; inspect its PID and confirm no child R process remains before
removing that lock. Worker-owned locks are otherwise cleaned automatically.

## Outputs and reproducibility

- `design.csv`, `regimes.csv`, `provenance.rds`, `session.txt`, `git_commit.txt`:
  exact design, seeds, software versions and source hashes.
- `raw/*.rds` + `.md5`: longest required residual path per independent trajectory;
  coarse observations and shorter horizons are derived, not separately stored.
- `results/*.rds`: estimates, full stage diagnostics, likelihood profile, timings,
  warnings, benchmark fits, source-data checksum and completion time.
- `cache/`: exact reusable numerical kernels and design-specific transfer tables.
- `checkpoints/`: completed likelihood clock rows for interrupted fits.
- `logs/`: per-job logs, execution errors and per-session scheduler timing records.
- `summary/`: long-form estimates, metrics, failures/joint recovery, paired
  contrasts with Monte Carlo intervals, pending IDs, LaTeX table and report.
- `figures/`: vector PDF recovery, information-growth and benchmark figures.

Publication figures are produced from completed data only. Empty or smoke outputs
are not manuscript evidence. Raw large arrays stay out of Git; preserve the full
output directory in a research archive. Reproduction requires both this repository
revision and the completed run directory, not just a screenshot or summary table.

## Interpretation

Error moments use finite parameter estimates and explicitly report their counts.
Joint tolerance recovery includes failed complete estimates in the denominator.
Thresholds are absolute 0.10 for alpha and relative 10% for positive parameters.
Bias/RMSE Monte Carlo SEs are empirical; paired intervals resample independent
trajectory IDs within a regime/design. No inference uses overlapping windows as
independent replications. With 24 paths/regime, tail quantiles and rare failures
are imprecise. Zero failures does not establish 99% reliability.

Keep oracle and numerical-refinement results separate. Preserve all poor cases.
An invalid ECF variance or unsupported transfer inversion is an estimator
breakdown, not automatically a coding/numerical error. Diagnose numerical causes
only from controlled refinements; leave unresolved cases unresolved. Execution
errors are logged as pending, not fabricated statistical estimates. No experiment
here validates drift removal, measurement error or irregular observations.

Main manuscript: ordinary recovery, paired information plots, difficult-regime
results and a brief numerical-stability assessment. Supplement: complete tables,
partial tail benchmarks, oracle contrasts and numerical diagnostics. Internal:
optimization and execution logs. No claim of superiority over a four-parameter
competitor is possible from the partial benchmarks below.

See `BENCHMARKS.md` for eligibility, formulas, defaults, and limitations.
