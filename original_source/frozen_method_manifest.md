# Frozen Four-Parameter Validation Manifest

## Freeze Declaration

Freeze date: 2026-07-24

This manifest fixes the statistical estimator, simulation panel, observation
designs, numerical transfer profiles, seeds, software environment, and
computational limits before any new accuracy result in the overarching
validation is opened.

The estimator order is:

\[
\widehat\alpha
\longrightarrow \widehat\sigma_2
\longrightarrow \widehat\beta
\longrightarrow \widehat\rho_{\rm FW}
\longrightarrow \widehat\sigma_1
=\widehat\beta/\widehat\rho_{\rm FW}.
\]

After this freeze, changes are limited to orchestration, logging,
checkpointing, aggregation, uncertainty calculations, classification, and
figures. Statistical settings and numerical transfer tables are immutable.

## Model And Simulator

The frozen simulation model is:

\[
dV_t=-\rho^2V_tdt+\rho\sqrt{1+V_t^2}\,dB_t,
\qquad
dR_t=\sigma_1dV_t+\sigma_2dL_t,
\]

where \(L_t\) is symmetric alpha-stable with unit characteristic exponent.
Simulation uses the existing Base R Euler implementation, `V_0 = 0`, a
Gaussian Brownian increment, and the existing Chambers-Mallows-Stuck stable
generator. The simulator seed is reset once per path from the frozen seed
manifest.

## Alpha Stage

Let \(A_{(1)}\ge\cdots\ge A_{(m)}\) be sorted nonzero absolute residual
increments. The Hill path is:

\[
\widehat\gamma_k
=\frac1k\sum_{i=1}^k\log A_{(i)}-\log A_{(k+1)},
\qquad
\widehat\alpha_k=1/\widehat\gamma_k.
\]

Frozen settings:

- running-median width: 5;
- alpha cap used in plateau selection: 1.98;
- 60 logarithmically spaced candidate values;
- candidate lower bound: 20;
- candidate upper bound: `min(m/100, 20000, m-10)`;
- stability score: standard deviation over indices `floor(k/2):k`;
- selected `k`: lowest score among admissible candidates;
- no truth-based or case-specific selection.

## Sigma2 Stage

For \(c_\alpha=\Gamma(\alpha)\sin(\pi\alpha/2)/\pi\), order-statistic
threshold \(u_k=A_{(k)}\), and \(p_k=(k-0.5)/m\), the frozen scale estimate at
each local tail index is:

\[
\widehat\sigma_{2,k}
=\left\{
\frac{p_ku_k^{\widehat\alpha}}
{2c_{\widehat\alpha}\,\operatorname{median}(\Delta)}
\right\}^{1/\widehat\alpha}.
\]

The reported estimate is the median over `k_hat - 20` through `k_hat + 20`,
truncated only by the available order statistics. The oracle-alpha version is
computed only as a diagnostic decomposition and never replaces the production
estimate.

## Beta Stage

The frozen beta stage is the validated conditional empirical-CF endpoint
likelihood. It integrates the block additive activity through a precomputed
endpoint-state Markov kernel and profiles beta at each rho grid point.

Observation settings:

- block horizon: 0.05;
- block size: `max(10, round(0.05 / median(dt)))`;
- standardized increments: `R / sqrt(dt)`;
- robust frequency scale: `median(abs(R / sqrt(dt)))`;
- frequency multipliers: 0.25, 0.50, 0.75;
- actual frequencies: multipliers divided by the robust scale;
- stable jump attenuation:
  `sigma2_hat^alpha_hat * abs(u)^alpha_hat *
  mean(dt^(1-alpha_hat/2))`;
- observation covariance: `conditional_diag`;
- likelihood evaluator: exact-equivalent chunked evaluator;
- groups per chunk: 16.

Endpoint and transition settings:

- rho profile grid: `0.7, 0.8, ..., 2.3`;
- endpoint states: 21;
- transition simulation blocks: 630,000;
- transition substeps per block: 32;
- maximum path nodes per endpoint pair: 31;
- transition kernel seed: 5000;
- beta bounds: `[0.20, 5.00]`;
- beta optimization: 25-point global grid in log beta followed by
  one-dimensional `optimize` inside the neighboring grid interval;
- optimization tolerance: `1e-4`;
- rho refinement: the existing local quadratic profile refinement;
- returned downstream coordinate: the refined `beta_hat`.

The production transition kernel has 17 rho models, is 58.466 MB, and has
SHA-256:

```text
3fb889a85e93af3a6c2782f02b06afae05a2ccf975d6ab430f078c94b5af032f
```

## Finite-Window Rho Stage

With \(n\) residual increments, the frozen rolling log-ECF rules are:

\[
k=\lfloor\sqrt n\rfloor,\qquad
H=k\Delta,\qquad
u=\{\widehat\beta\sqrt{\log(e+n)}\}^{-1}.
\]

The stable deconvolution exponent is:

\[
\lambda
=\widehat\sigma_2^{\widehat\alpha}|u|^{\widehat\alpha}
\Delta^{1-\widehat\alpha/2}.
\]

The local empirical-CF floor is exactly \(k^{-1/2}\). The frozen numerator
uses adjacent backward and forward log-ECF variance estimates, the existing
first-order local ECF noise correction, and chunk size 100,000. The
denominator is:

\[
\widehat D
=4\Delta\sum_j
\left(1-\frac{\widehat\beta^2}{\widehat C_{j,+}}\right).
\]

The finite-window estimate solves the one-dimensional equation:

\[
\frac{\widehat I_{\log}}{\widehat D}
=r^2a(q,r^2H),
\qquad
q=\{\log(e+n)\}^{-1},
\]

by inversion of the frozen monotone numerical transfer on
`r in [0.4, 3.5]`. No guard, fallback, reranking rule, or truth-based
selection is permitted.

Finally:

\[
\widehat\sigma_1=\widehat\beta/\widehat\rho_{\rm FW}.
\]

## Transfer Profiles

Every distinct observation cell has a matching transfer profile. Each profile
uses:

- rho grid `0.4, 0.5, ..., 3.5`;
- 50,000 paired primary paths, seed 941001;
- 30,000 paired independent paths, seed 941002;
- nested 64, 128, and 256 step resolutions;
- first-order Richardson evaluation from 128 and 256;
- chunk size 1,000.

| n | T | q | H | profile status |
|---:|---:|---:|---:|---|
| 1,000,000 | 50 | 0.07238240 | 0.05000 | pass |
| 2,500,000 | 25 | 0.06788036 | 0.01581 | pass |
| 4,000,000 | 50 | 0.06578166 | 0.02500 | pass |
| 5,000,000 | 50 | 0.06483004 | 0.02236 | pass |
| 10,000,000 | 50 | 0.06204207 | 0.01581 | pass |
| 10,000,000 | 100 | 0.06204207 | 0.03162 | pass |

Across the six profiles:

- every corrected map is strictly increasing;
- maximum fixed-point inversion error is zero on the interpolation map;
- maximum independent-replication discrepancy is 0.26642 standard errors;
- maximum 128-to-256 transfer difference is 0.005714;
- maximum smallest-tau distance from one is 0.02601.

The authoritative checks are in
`overarching_transfer_profile_validation.csv`.

## Scenario And Observation Bounds

The twelve frozen regimes use:

```text
alpha  in [1.08, 1.58]
sigma2 in [0.60, 1.40]
sigma1 in [0.70, 1.60]
rho    in [0.90, 2.20]
beta   in [0.8925, 1.92]
```

The 254 frozen run specifications contain:

- 96 production development paths;
- 12 untouched production paths;
- 72 lower-resolution infill paths;
- 54 fixed-spacing long-span paths;
- 20 targeted heavy-tail stress paths.

All 254 run keys, case identifiers, and simulation seeds are unique. Exact
values and execution priorities are fixed by:

- `frozen_scenario_design.csv`;
- `frozen_seed_manifest.csv`;
- `frozen_run_design.csv`.

## Software And Host

```text
R:                 4.5.2 (2025-10-31 ucrt)
R platform:        x86_64-w64-mingw32
yuima installed:   1.15.34
statistical code:  Base R only; yuima is not called by the frozen estimator
OS:                Microsoft Windows 10 Home 10.0.19045
CPU:               Intel Core i9-9900K, 8 physical / 16 logical
physical memory:   15.92 GB
worker ceiling:    3
checkpoint drive:  D:
```

The three-worker ceiling was selected from the smaller of the host core limit
and memory-safe capacity. A prior production worker peak was approximately
1.36 GB, while about 7 GB was free at freeze.

## Source Integrity

The complete SHA-256 inventory is in `frozen_method_hashes.csv`. Its own
SHA-256 is:

```text
001e2b5bee09eff6e3fdfb1888b8e0003cf3c57f89f650be7a6dffbbd9b54bab
```

Core statistical source hashes:

| file | SHA-256 |
|---|---|
| `work/occupation_transition_clock_prototype.R` | `218bcd323ba2e9b188b21cd182d316b99f1ed6a1d55ba0bd2e2c4b449f6fed6e` |
| `work/diagnose_markov_additive_endpoint_hmm.R` | `35009392761613f2c61fe569001382492de9add14a09289216ce84c19b5c5334` |
| `work/diagnose_exact_cf_endpoint_hmm.R` | `6a80f94705f34a114359261c5b91f5b53c215f27e0a1b13db16eb431fc07ef08` |
| `work/chunked_exact_endpoint_likelihood.R` | `55ea5431dff46f4f0fc000c8c0496eb504e0103d02cf86e635e763008380d6b3` |
| `work/run_conditional_cf_endpoint_full_pipeline.R` | `f09d5229e07bfec3e4a509241baeb3fe38df85548a9af7604ecb5d406ccda07d` |
| `rho_level4_oracle_residual_estimator.R` | `c8690657da478294bfcade7fd5cc0df2d849ea6ebe54b40a38deb58d685fd87c` |
| `finite_window_corrected_rho_estimator.R` | `d269782cf88542953d7526eb0e136fddb7a9a0a29886470f6ab84937787e23b9` |
| `work/run_finite_window_corrected_full_pipeline.R` | `57042881c272bfa8bdfdd9a93179145aade44aa3cd3a0b38d135f22d6969b28e` |

## Execution And Retention

Runs execute in priority order:

1. production feasibility;
2. production expansion;
3. fixed-T infill;
4. fixed-spacing long span;
5. heavy-tail stress;
6. untouched production.

Completed canonical results are never repeated. Failed attempts receive
attempt-numbered start, status, stdout, and stderr records and are not retried
without explicit `-RetryFailed` authorization. Every statistical failure,
nonfinite output, warning, and extreme path is retained.

The untouched priority-6 panel cannot be opened until the development and
stress panels are classified. Any statistical code or profile change after
this freeze invalidates that panel.
