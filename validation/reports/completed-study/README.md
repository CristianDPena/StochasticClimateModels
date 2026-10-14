# Completed numerical validation study

All 256 planned evaluations were run: 168 main recovery fits, 48 paired
sample-size/horizon fits, and 40 prespecified diagnostics. All raw outcomes,
including failures and poor estimates, are retained. The estimator was not tuned
or changed in response to these experiments. Of 216 principal fits, 215 returned
all four estimates and one reported a downstream variance failure. All 40
diagnostic evaluations returned estimates. These counts describe computation;
accuracy is assessed separately below.

## Design and interpretation

The seven regimes, seeds, horizons, sample sizes, and numerical variants are in
[design.csv](design.csv) and [regimes.csv](regimes.csv). The main study uses
T=50 and n=10 million, with 24 independent trajectories per regime. The paired
information study uses the first 12 trajectory IDs in B and H, comparing
n=2.5 million to n=10 million at fixed T=50, and T=50 to T=100 at fixed spacing
5e-6. Noise and latent simulation streams are fixed and shared across the
prescribed paired comparisons. The simulator starts at V(0)=0. The frozen
estimator configuration and its numerical approximations are documented in
[the protocol](../../README.md); no theoretical consistency, rate, or
identifiability assertion follows from this experiment.

Error moments use finite estimates, with finite counts explicitly retained.
Undefined downstream estimates are failures in tolerance and joint-recovery
calculations. Joint recovery means absolute alpha error below 0.10 and relative
error below 10% for each of sigma2, sigma1, and rho. The CSVs include bias,
standard deviation, RMSE, MAE, median and 90th-percentile absolute errors,
Monte Carlo standard errors, and Wilson intervals for tolerance fractions.
With 24 independent trajectories per main regime, rare failures and tail error
quantiles remain imprecise. Paired uncertainty uses 999 trajectory-level
bootstrap draws within each regime/design; paired finite counts are reported.

## Main recovery results

The table gives absolute RMSE for alpha and relative RMSE for the other
parameters. All four estimates are finite in every main-study fit.

| Regime | alpha RMSE | sigma2 rel. RMSE | sigma1 rel. RMSE | rho rel. RMSE | Joint recovery |
| --- | --- | --- | --- | --- | --- |
| B | 0.0657 | 7.23% | 2.24% | 2.22% | 18/24 |
| Aminus | 0.1093 | 10.29% | 2.46% | 2.04% | 14/24 |
| Aplus | 0.0659 | 3.58% | 5.82% | 3.73% | 20/24 |
| Jminus | 0.1181 | 8.94% | 2.21% | 2.13% | 12/24 |
| Cminus | 0.0442 | 5.07% | 2.66% | 2.71% | 23/24 |
| Sminus | 0.0442 | 5.67% | 3.45% | 3.07% | 22/24 |
| H | 0.0913 | 2.70% | 7.99% | 5.32% | 17/24 |

Baseline amplitude and clock estimates have relative RMSE near 2.2%; tail-index
and jump-scale recovery are more variable. Alpha RMSE is largest in Jminus and
Aminus, while the mixed regime H has the largest relative sigma1 and rho RMSE.
These are observed finite-sample difficulties, not proof of structural weak
identification. Cminus and Sminus perform well under the prescribed settings.
The [recovery figure](figures/recovery.pdf) displays every finite estimate,
including outliers; complete parameter-wise results are in
[metrics.csv](summary/metrics.csv).

## Sample size and observation horizon

The following are absolute RMSEs on the same 12 trajectory IDs per regime.
The last column gives the finite count for sigma1 and rho; alpha and sigma2
have 12 finite estimates in every row.

| Regime | Design | alpha | sigma2 | sigma1 | rho | Downstream finite |
| --- | --- | --- | --- | --- | --- | --- |
| B | T=50, n=2.5M | 0.0937 | 0.0611 | 0.0278 | 0.0462 | 12/12 |
| B | T=50, n=10M | 0.0695 | 0.0711 | 0.0225 | 0.0342 | 12/12 |
| B | T=100, n=20M | 0.1166 | 0.1300 | 0.0183 | 0.0250 | 12/12 |
| H | T=50, n=2.5M | 0.1205 | 0.0441 | 0.1067 | 0.2204 | 12/12 |
| H | T=50, n=10M | 0.0946 | 0.0352 | 0.0539 | 0.1087 | 12/12 |
| H | T=100, n=20M | 0.1292 | 0.6280 | 0.0421 | 0.0967 | 11/12 |

Increasing observation density improves alpha, sigma1, and rho RMSE in both
regimes; sigma2 improves in H but not B. Increasing horizon improves the
reported finite-output RMSE for sigma1 and rho, but worsens alpha and sigma2
RMSE in both regimes. The H horizon comparison is conditional on 11 finite
downstream outputs, whereas its reference row has 12, so it cannot support an
unqualified claim of improved reliability. The [information figure](figures/information_growth.pdf)
shows marginal RMSE +/- 1.96 Monte Carlo SE, not confidence intervals for the
paired differences. Those differences and their bootstrap intervals are in
[paired_summary.csv](summary/paired_summary.csv). Two information levels and
12 paths per regime do not establish a convergence rate or consistency.

The failed case is H trajectory 10 at T=100 (job v0238). It returned
alpha=1.969058, sigma2=3.573359 and beta=5.000000, followed by
`rho:uncorrected_spot_variance_failure`; sigma1 and rho are undefined.
Its upstream estimates remain in their error moments, so the large sigma2
error is not omitted. Beta reaches the configured upper bound. This record
shows where the pipeline failed; it does not alone distinguish statistical
error, approximation error, or a combination. No post hoc retuning or
selectively chosen rescue fit was performed. In particular, the prespecified
oracle/refinement subset does not include this trajectory, and its causal
failure mechanism remains unresolved.

## Targeted numerical checks

Each row below compares the four prescribed B/H datasets with their ordinary
fit. Entries are maximum absolute changes over these four datasets.

| Check | beta change | rho change | sigma1 change |
| --- | --- | --- | --- |
| Kernel substeps 32 to 64 | 0.002421 | 0.012381 | 0.002526 |
| Simulator substeps 4 to 8 | 0.001033 | 0.003460 | 0.001532 |
| Refined transfer | 0.000000 | 0.010016 | 0.004937 |
| Independent transfer seed | 0.000000 | 0.008108 | 0.003925 |

The shifts are small on this subset, supporting local numerical stability of
the implemented configuration. Kernel refinement changes within-block substeps
only; states, path representatives, clock grid, and simulated transitions per
state are held fixed. Kernel paths are not coupled across resolutions, so the
shift includes Monte Carlo variation. Simulator refinement uses coupled Brownian
increments. Transfer refinements reuse beta, making its zero change intentional.
These checks do not establish convergence of every likelihood/quadrature
approximation, validate all regimes, or explain the failed horizon fit.
A broader state/node/transition-budget study would be needed for stronger
numerical-convergence claims.

## Upstream-error localization

The following downstream absolute RMSEs use the same four prescribed
trajectories in each diagnostic regime. Oracle rows are diagnostics, not
practically available competitors.

| Regime | Upstream quantities | sigma1 RMSE | rho RMSE |
| --- | --- | --- | --- |
| B | Estimated upstream | 0.0318 | 0.0469 |
| B | True alpha, sigma2 | 0.0292 | 0.0444 |
| B | True alpha, sigma2, beta | 0.0258 | 0.0377 |
| Jminus | Estimated upstream | 0.0235 | 0.0311 |
| Jminus | True alpha, sigma2 | 0.0158 | 0.0262 |
| Jminus | True alpha, sigma2, beta | 0.0141 | 0.0218 |
| H | Estimated upstream | 0.0440 | 0.0791 |
| H | True alpha, sigma2 | 0.0100 | 0.0328 |
| H | True alpha, sigma2, beta | 0.0119 | 0.0380 |

Replacing estimated alpha and sigma2 by truth reduces downstream error on these
small subsets, most visibly in H. Supplying true beta provides a further
reduction in B and Jminus but does not monotonically improve H. The results
therefore support an upstream-error contribution without suggesting that all
downstream error is upstream error or that an oracle must improve every finite
sample. Residual variation remains when all upstream quantities are supplied.
Four paths per regime permit localization, not a precise error decomposition.

## Tail-estimator comparisons

Hill uses k=floor(sqrt(n)); generalized-Pareto maximum likelihood uses the
99th percentile threshold and evir 1.7-4 defaults. Both receive the same observed
residual increments as the fixed method. Thresholds were declared in advance
and were not changed after observing outcomes. They estimate the tail index
and, through the stable-tail conversion, sigma2; neither estimates sigma1 or
rho. Their tail approximation is applied to dependent composite residuals,
not directly observed iid stable jumps. See [benchmark specifications](../../BENCHMARKS.md).

| Regime | Fixed alpha RMSE | Hill alpha RMSE | GPD alpha RMSE | Hill finite sigma2 | GPD finite sigma2 |
| --- | --- | --- | --- | --- | --- |
| B | 0.0657 | 1.0765 | 2.8672 | 11/24 | 1/24 |
| Aminus | 0.1093 | 1.5422 | 4.5693 | 5/24 | 0/24 |
| Aplus | 0.0659 | 0.4708 | 1.2229 | 16/24 | 1/24 |
| Jminus | 0.1181 | 1.8497 | 2.9111 | 0/24 | 0/24 |
| Cminus | 0.0442 | 0.7592 | 1.9510 | 18/24 | 1/24 |
| Sminus | 0.0442 | 0.6139 | 1.5282 | 20/24 | 17/24 |
| H | 0.0913 | 0.3449 | 0.5197 | 21/24 | 12/24 |

The fixed method has lower empirical alpha RMSE than these two declared
implementations in every main regime. Many competitor tail-index estimates
fall outside the admissible stable interval, making the stable-scale conversion
undefined. These indices are retained in alpha error moments; undefined scales
remain in tolerance denominators. Conditional scale RMSEs must be read with
their finite counts, rather than ranked as if all methods succeeded equally.
This is evidence against these particular generic tail procedures at the
specified settings for these observations, not a universal comparison against
all threshold choices or established four-parameter estimators. The benchmark
figure and all paired results belong in the supplement.

## Recommended scholarly presentation

Use the seven-regime recovery table and distribution figure as the central
validation. Include the paired information-growth results and the failed
horizon case in the main text: a blanket claim that accuracy improves with
information is contradicted by this study. Summarize the numerical shifts
briefly in the main text and provide their complete paired table in the
supplement. Keep oracle results, full benchmark tables, Monte Carlo uncertainty,
and all finite/failure counts in the supplement. Execution logs remain internal.

The defensible conclusion is useful four-parameter recovery at the principal
high-frequency design, with stronger downstream precision than tail-parameter
precision in typical regimes, accompanied by nonuniform information-growth
behavior and a documented downstream failure. Further work should investigate
that failure with a separately declared diagnostic study, increase replication
for rare failures, assess additional quadrature dimensions, and evaluate drift
removal, measurement noise, irregular sampling, and eventual observational-data
fitting. This validation does not warrant theoretical or universal robustness
claims.

## Reproducibility and resources

The [integrity audit](audit/README.md) verifies 168 raw trajectory files and
sidecars, all 256 design-matched result records, their raw-data links, frozen
source/framework hashes, and summary statuses. The raw run directory is
`validation/results/laptop-v1` (approximately 29 GiB including checkpoints and
caches); it is retained locally and excluded from Git. This compact published
folder contains every summary row, figures, design, software versions, source
provenance, and the raw/result checksum inventory. It is not a substitute for
archiving the full raw directory. No external raw-data archive was uploaded.

Two workers were used throughout fitting. After correcting the scheduler's
handling of macOS suspension, recorded bounded sessions took approximately
1.38 to 6.55 hours. These elapsed times include suspension and are not CPU
benchmarks. The initial scheduler used a clock that could exclude suspension;
its timing must not be included in a precise total-runtime claim. The initial
session was interrupted once, and two unfinished jobs resumed from saved
likelihood rows. All completed estimates were retained. See
[execution notes](execution/execution_notes.md). Actual storage exceeded the
initial 16-25 GB estimate because retained checkpoints and caches also occupy
space. The study required multiple sessions, not a single sub-12-hour run.

To regenerate the published package from the retained raw run:

```sh
Rscript validation/run.R --command=summarize
Rscript validation/figures.R
python3 validation/export_results.py
```

The exporter invokes the complete integrity audit before copying results.
