# Initial production calibration — two paths only

**The full validation study is not complete.** These records establish successful
execution and provide an initial resource estimate. They must not be used as the
paper's numerical validation or as evidence of comparative superiority.

Both paths use regime B: alpha=1.30, sigma2=1.00, sigma1=1.00, rho=1.50,
T=50 and n=10,000,000. They are the first two prespecified production trajectories,
not paths selected for favorable accuracy, and remain part of the study.

| Replication | alpha estimate | sigma2 estimate | beta estimate | rho estimate | sigma1 estimate | Fit seconds |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 1.343596 | 0.995557 | 1.490425 | 1.526611 | 0.976297 | 286.655 |
| 2 | 1.340848 | 0.969795 | 1.492238 | 1.426415 | 1.046146 | 278.894 |

Both fixed-method fits returned success. The two-worker session elapsed time was
296.2 seconds including simulation, fit, benchmark, and summary work. A mid-run
snapshot measured about 1.4 GB resident memory per worker; this is not peak memory.
The later standalone final plot rendering is not included in that scheduler time.
R 4.6.1 and evir 1.7-4 were used on the user's Apple arm64 MacBook Air.

The fixed-rank Hill comparator produced alpha=1.6231 and 2.1118; the second scale
estimate is undefined under stable-tail conversion. GPD ML at the fixed 99%
threshold produced alpha=2.5328 and 3.1969, so both stable scale conversions are
undefined. Those outputs are retained. No threshold or method was changed after
seeing these outcomes. At this small sample of paths the comparisons provide no
basis for a general ranking; threshold contamination is a possible explanation
that the complete prespecified study must assess.

The full numerical design has 216 principal jobs and 40 diagnostics. Only 2/216
principal jobs and 0/40 diagnostic jobs are complete. `pending.csv` lists every
remaining job. Oracle smoke tests elsewhere are computational checks, not completed
scientific oracle diagnostics.

- `estimates.csv`: all ordinary and benchmark estimates for these same datasets.
- `metrics.csv`: finite-error summaries, counts, and descriptive Monte Carlo measures.
- `failures_and_joint_recovery.csv`: full-estimator status and joint tolerance outcomes.
- `paired_*.csv`: paired records and summaries; with two paths, intervals are unstable.
- `recovery.pdf`, `benchmarks.pdf`: vector figure output, with finite/total counts.
- `information_growth.pdf`: explicitly states that these results are pending.
- `recovery_table.tex`: a compilable table fragment requiring booktabs.
- scheduler JSON: measured session timing (per-job elapsed includes waiting for
  the other worker when collecting finished processes; use Fit seconds for exact
  estimator runtime).

The large raw residuals and complete per-job RDS records are retained locally under
`validation/results/laptop-v1/`, excluded from Git. Their checksums are stored with
the records. Archive that directory after completion; the CSV summaries alone are
not a substitute for the reproducible study record.

Resume with `python3 validation/schedule.py --batch=recovery1`. It skips these two
completed jobs. The README at `validation/README.md` gives the remaining batches.
