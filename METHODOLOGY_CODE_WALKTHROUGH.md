# Methodology Code Walkthrough

This document describes what the frozen source computes, including its
boundary rules and numerical approximations. It does not replace the source.

## A. Alpha Stage

### Input

The estimator receives residual increments `R` and aligned time increments
`dt`, not `X`:

- call: `work/run_conditional_cf_endpoint_full_pipeline.R:164-170`
- function: `work/occupation_transition_clock_prototype.R:116-152`

It retains indices where:

```text
R is finite, R != 0, dt is finite, dt > 0
```

and sorts `A=abs(R)` in decreasing order
(`work/occupation_transition_clock_prototype.R:118-122`).

### Hill path

For descending `A[1],...,A[m]`, code lines `46-52` calculate, for
`k=1,...,m-1`:

```text
gamma_hat[k] =
  (sum_{i=1}^k log(A[i]) - k*log(A[k+1])) / k

alpha_raw[k] = 1 / gamma_hat[k]
```

The raw alpha path is smoothed with a running median of default width 5,
forced to an odd width, using `runmed(...,endrule="median")` at lines
`123-131`.

### Candidate grid and stability selection

`choose_hill_plateau`, lines `54-87`, uses:

- `k_min=20`
- `k_max=min(floor(m/100),20000,m-10)`
- 60 logarithmically spaced candidates, rounded to unique integers
- segment `floor(k/2):k`
- stability score equal to the population standard deviation computed from
  segment first and second moments
- candidates with `alpha_path[k] > 1.98` are skipped
- first strictly lowest score wins

Boundary behavior:

- if `k_max<20`, select
  `min(max(2,floor(m/10)),m)` and flag `small_sample=TRUE`;
- if all candidates are inadmissible, select `k=20`;
- the outer estimator requires at least 100 finite nonzero residuals;
- it stops if final alpha is nonfinite or outside `(0,2)`.

### Returned values

- `alpha_hat`
- `k_hat`
- `hill_score`
- candidate range
- small-sample flag
- sigma2 diagnostics

Code: `work/occupation_transition_clock_prototype.R:143-150`.

## B. Sigma2 Stage

### Formula

`estimate_sigma2_tail_robust`,
`work/occupation_transition_clock_prototype.R:89-114`, computes:

```text
c_alpha = Gamma(alpha_hat) sin(pi alpha_hat/2) / pi
p_k     = (k - 0.5) / m
u_k     = A[k]

sigma2_k =
  {p_k u_k^alpha_hat /
   (2 c_alpha median(dt))}^{1/alpha_hat}
```

The reported estimate is the median of finite `sigma2_k` over:

```text
k_lo = max(2, k_hat-20)
k_hi = min(m-1, k_hat+20)
```

### Dependence and corrections

- Uses `alpha_hat` directly.
- Uses the median aligned positive `dt`.
- Uses no latent-state or dependence correction.
- Uses no truth in production.
- The caller stops if the result is nonfinite or nonpositive.

### Diagnostics

Returns `k_hat`, `k_lo`, `k_hi`, robust `dt`, `c_alpha`, and window size.

The outer runner’s true-alpha recomputation at
`work/run_finite_window_corrected_full_pipeline.R:131-161` is evaluation only.

## C. Beta Stage

### Block observations

Function:
`block_cf_log_observations`,
`work/occupation_transition_clock_prototype.R:264-370`.

Frozen settings:

- `block_horizon=0.05`
- `block_size=max(10,round(0.05/median(dt)))`
- at least five complete blocks
- `Z=R/sqrt(dt)`
- scale `s=median(abs(Z))`
- fallbacks if `s` is invalid: MAD, then standard deviation
- frequencies `u_m=(0.25,0.50,0.75)/s`

For block `j` and frequency `u_m`:

```text
phi_hat[j,m] = mean_k exp(i u_m Z_k)
y_raw[j,m]   = log(max(|phi_hat[j,m]|, 1e-8))
jump_att[j,m] =
  sigma2_hat^alpha_hat |u_m|^alpha_hat
  mean_k dt_k^(1-alpha_hat/2)
y_dejumped[j,m] = y_raw[j,m] + jump_att[j,m]
```

The function also delta-method estimates an empirical covariance and a
log-modulus bias. In the frozen `conditional_diag` likelihood:

- `y_dejumped`, not `y_dejumped_bias_corrected`, is passed to the likelihood
  at `work/diagnose_exact_cf_endpoint_hmm.R:1074-1078`;
- a model-based finite-sample log-modulus bias and diagonal variance are
  instead included in `exact_cf_conditional_moments`, lines `254-334`;
- the empirical covariance is used in the projected beta initialization,
  lines `995-1005`, but not as the production emission covariance.

### State representation

The latent endpoint state is a 21-bin discretization of
`Q=1+V^2`:

- `make_q_endpoint_breaks`:
  `work/diagnose_markov_additive_endpoint_hmm.R:220-240`
- bins are stationary quantiles from 300,000 draws with default seed 881
- lower and upper cut bounds are forced to 1 and infinity

For each candidate rho, the cached model contains:

- endpoint transition probabilities;
- active start/end state pairs;
- up to 31 representative within-block Q paths per active pair;
- an initial transition from `Q_0=1`;
- 32 within-block substeps.

Generator:
`work/diagnose_exact_cf_endpoint_hmm.R:36-161`.

### Conditional observation model

For a candidate within-block path `q_s`, beta, frequency `u`, and jump
attenuation `lambda_u`, lines `254-334` use:

```text
phi_s(u) = exp(-0.5 beta^2 u^2 q_s - lambda_u)
mean_phi = average_s phi_s(u)

model mean =
  log(mean_phi) + lambda_u + log-modulus-bias
```

The diagonal model variance is calculated from `phi(u)` and `phi(2u)` and
divided by block size. Values are floored at machine epsilon.

The likelihood mixes representative path-node Gaussian emissions within each
endpoint pair, then performs an HMM forward recursion. The production
memory-bounded implementation is:

- `work/chunked_exact_endpoint_likelihood.R:3-110`

It is intended to be arithmetically equivalent to the full-matrix evaluator.

### Beta profile

`estimate_exact_cf_endpoint_hmm`,
`work/diagnose_exact_cf_endpoint_hmm.R:922-1210`:

1. evaluates candidate rho values `0.7,0.8,...,2.3`;
2. profiles beta separately at each rho;
3. searches 25 equally spaced points in log beta over `[0.20,5.00]`;
4. if the best grid point is interior, calls `optimize` in the neighboring
   interval with tolerance `1e-4`;
5. stores the maximum likelihood and beta for each rho;
6. performs a three-point quadratic fit in log rho around the best profile
   point;
7. if the quadratic is concave and its vertex lies inside those neighbors,
   interpolates beta to that vertex; otherwise retains the best grid point.

The returned `beta_hat` is the refined beta coordinate. The endpoint
likelihood also returns an internal rho/sigma1 split, but the authoritative
outer runner saves it only as a diagnostic and replaces it.

### Bounds and failure behavior

- beta is bounded to `[0.20,5.00]`;
- nonfinite objectives are mapped to a very large finite penalty;
- a boundary rho maximum produces a warning, not an alternate estimator;
- missing/incompatible kernel cache triggers in-memory kernel construction;
- exact-signature checkpoints may be reused.

## D. Corrected Rho Stage

### Frozen normalization

Outer call:
`work/run_finite_window_corrected_full_pipeline.R:54-67`.

```text
k = floor(sqrt(n))
H = k dt
u = 1 / {beta_hat sqrt(log(e+n))}
q = (u beta_hat)^2 = 1/log(e+n)
lambda =
  sigma2_hat^alpha_hat |u|^alpha_hat
  dt^(1-alpha_hat/2)
```

Stable exponent code:
`rho_level4_oracle_residual_estimator.R:39-42,150-153`.

The rho stage requires one positive scalar `dt`; the runner supplies
`sim$dt[1]`.

### Local empirical CF and deconvolution

At every valid center:

- backward window: the preceding `k` cosine observations;
- forward window: the following `k` cosine observations;
- empirical cosine CF floor: `k^(-1/2)`;
- stable deconvolution multiplier: `exp(lambda)`.

For either side:

```text
CF_effective = max(CF_raw, k^(-1/2))
CF_tilde     = exp(lambda) CF_effective
C_hat        = -2 log(CF_tilde) / u^2
```

Code:
`rho_level4_oracle_residual_estimator.R:84-117,150-186`.

The code records, but does not clip away:

- nonpositive raw ECFs;
- floor activations;
- deconvolved ECFs above one;
- nonpositive C estimates;
- C estimates below beta squared.

### Numerator

For valid positive backward/forward C:

```text
dlogC = log(C_forward) - log(C_backward)

h_delta(u,C,lambda,alpha) =
  2/u^4 [
    expm1(u^2 C + 2 lambda) +
    expm1(-u^2 C + 2 lambda - 2^alpha lambda)
  ]

h_log = h_delta / C^2

I_log_hat =
  sum_j [
    3/(2k) dlogC_j^2 -
    3/k^2 h_log(u,C_forward_j,lambda,alpha)
  ]
```

Code:
`rho_level4_oracle_residual_estimator.R:44-57,211-234`.

The second term is the analytical local ECF sampling-noise subtraction.

### Denominator

```text
D_hat =
  4 dt sum_j (1 - beta_hat^2/C_forward_j)
```

Code:
`rho_level4_oracle_residual_estimator.R:235-246`.

The implementation demands that every local C and every required correction
term be valid. Any invalid local value makes numerator/denominator final
values `NA` and assigns a failure status.

### Uncorrected coordinate

```text
rho_uncorrected = sqrt(I_log_hat / D_hat)
```

only if status is success and both numerator and denominator are positive
(`rho_level4_oracle_residual_estimator.R:251-285`).

### Finite-window transfer object

The generic unit-clock transform uses stationary
`Z~t_3/sqrt(3)` and `X=asinh(Z)`. It evolves:

```text
dX_s = -1.5 tanh(X_s) ds + dW_s
```

and evaluates:

```text
A_q,tau = tau^(-1) integral exp{-q cosh(X_s)^2/2} ds
Q_q,tau = -2 log(A_q,tau)/q
L_q,tau = log(Q_q,tau)
tau     = rho^2 H
```

Implementation:
`finite_window_corrected_rho_estimator.R:24-97,99-258`.

The transfer estimate is:

```text
numerator_rate_unit   = (3/tau) E[0.5(L_1-L_2)^2]
denominator_rate_unit = 4 E[0.5{1-1/Q_1 + 1-1/Q_2}]
transfer              = numerator_rate_unit/denominator_rate_unit
corrected_ratio       = rho^2 transfer
```

The paired paths share the same stationary initial state and have independent
future Brownian innovations.

Frozen profiles use:

- rho `0.4:3.5` by `0.1`;
- primary 50,000 path pairs, seed 941001;
- independent 30,000 path pairs, seed 941002;
- 64/128/256 nested time steps;
- chunk size 1,000.

### Richardson evaluation and inversion

`rho_fw_evaluation_profile`,
`finite_window_corrected_rho_estimator.R:313-356`, takes the finest and next
coarsest resolution and sets:

```text
transfer_evaluation =
  transfer_fine + (transfer_fine - transfer_coarse)
```

It requires a final 2:1 resolution ratio, positive finite transfer, and a
strictly increasing corrected map.

`rho_fw_invert_ratio`, lines `358-403`, computes:

```text
observed_ratio = I_log_hat / D_hat
rho_hat = linear interpolation inverse of
          rho -> rho^2 transfer_evaluation(rho)
```

It does not extrapolate or clip to the profile boundary. Values outside the
map return `NA` with `correction_domain_failure`.

## E. Sigma1 Stage

The final estimate is exactly:

```text
sigma1_hat = beta_hat / corrected_rho_hat
```

Code:
`finite_window_corrected_rho_estimator.R:418-423`.

There are no sigma1 bounds, shrinkage, fallback, or independent
postprocessing steps. The beta/rho cancellation fields in
`work/run_finite_window_corrected_full_pipeline.R:162-177` are
simulation-evaluation diagnostics and do not alter the estimate.

## Implemented Floors, Bounds, and Failure Rules

| Stage | Rule | Code |
|---|---|---|
| Alpha | discard zero/nonfinite R and nonpositive/nonfinite dt | occupation prototype `118-122` |
| Alpha | at least 100 observations | `122` |
| Alpha | cap 1.98 during plateau selection | `54-84` |
| Alpha | small-sample and no-admissible-candidate fallbacks | `58-62,84` |
| Sigma2 | local k window clipped to `[2,m-1]` | `95-102` |
| Beta | block size at least 10 and at least 5 blocks | `264-276` |
| Beta | residual-scale fallback median to MAD to SD | `280-290` |
| Beta | ECF modulus floor `1e-8` | `315` |
| Beta | empirical covariance diagonal floor `1e-10` | `347-350` |
| Beta | beta bounds `[0.20,5.00]` | exact endpoint `926,1082` |
| Beta | machine-epsilon model variance/log denominators | exact endpoint `297-303` |
| Beta | quadratic rho refinement only if concave/interior | Markov file `1171-1193` |
| Rho | local ECF floor `k^-1/2` | rho Level-4 `84-117,177` |
| Rho | all local C/h values must be valid | rho Level-4 `240-265` |
| Rho | numerator and denominator must be positive | rho Level-4 `261-270` |
| Correction | exact q/H match tolerance `1e-12` | outer runner `42-52` |
| Correction | monotone map and 2:1 resolution required | correction `313-355` |
| Correction | no out-of-domain extrapolation | correction `358-403` |
| Sigma1 | finite positive rho required; otherwise NA | correction `418-423` |

