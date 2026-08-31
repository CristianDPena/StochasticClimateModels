# Finite-Window Semigroup Correction For The Rolling Log-ECF Statistic

## 1. Scope And Claim Status

This document derives a finite-window target from the statistic implemented
in `rho_level4_oracle_residual_estimator.R`. It does not replace the rolling
statistic with a simpler proxy.

The derivation separates:

- exact identities;
- an infill limit with nonzero physical window width;
- a stationary finite-window estimating equation;
- a numerical semigroup expectation;
- approximations and proof obligations.

The proposed correction is model based. It contains no coefficients fitted
to rho estimation errors.

## 2. Continuous Functional Targeted By The Local ECF

Let `Delta -> 0`, `k -> infinity`, and `H = k Delta > 0` be held fixed.
For a window `[a,a+H]`, stable deconvolution and the local martingale law of
large numbers give the finite-window target

\[
\mathcal C_{u,H}[C;a]
=-\frac{2}{u^2}\log\left[
\frac1H\int_a^{a+H}
\exp\left\{-\frac{u^2C_s}{2}\right\}ds
\right].
\tag{2.1}
\]

The stable term is absent from (2.1) because its conditional characteristic
factor is removed before the logarithm. The finite-variation drift of `Y`
contributes a phase of order `sqrt(Delta)` to each normalized increment and
vanishes in this limit.

Equation (2.1) is exact for the limiting local statistic at fixed `u` and
fixed `H`. It is a log-Laplace or entropic average of the activity path, not
pointwise activity.

If additionally `u -> 0`,

\[
\mathcal C_{u,H}[C;a]
\longrightarrow
\frac1H\int_a^{a+H}C_sds.
\tag{2.2}
\]

Thus even at small frequency the finite-window object is block-average
activity, not `C_a`.

Define

\[
G_{u,H}(a)=\log\mathcal C_{u,H}[C;a].
\]

Ignoring the separately estimated ECF sampling-noise term, the continuous
finite-window limit of the implemented squared component is

\[
J_{u,H,T}
=\frac{3}{2H}
\int_H^{T-H}
\{G_{u,H}(t)-G_{u,H}(t-H)\}^2dt.
\tag{2.3}
\]

The discrete analytical `h_log` subtraction remains necessary at the
implemented sampling frequency. It removes local empirical-CF noise before
the statistic is interpreted through (2.3).

The denominator has the corresponding finite-window limit

\[
D_{u,H,T}
=4\int_H^{T-H}
\left\{1-\frac{\beta^2}{\mathcal C_{u,H}[C;t]}\right\}dt.
\tag{2.4}
\]

The correction branch leaves the implemented denominator unchanged.

## 3. Unit-Clock Reduction

Let the unit-clock diffusion satisfy

\[
dZ_s=-Z_sds+\sqrt{1+Z_s^2}\,dW_s,
\tag{3.1}
\]

so that

\[
V_t\overset{d}=Z_{\rho^2t}.
\]

Let

\[
q=(u\beta)^2,\qquad \tau=\rho^2H.
\]

For a unit-clock path define

\[
A_{q,\tau}
=\frac1\tau\int_0^\tau
\exp\left\{-\frac q2(1+Z_s^2)\right\}ds,
\tag{3.2}
\]

\[
Q_{q,\tau}
=-\frac2q\log A_{q,\tau},
\qquad
L_{q,\tau}=\log Q_{q,\tau}.
\tag{3.3}
\]

Then

\[
\mathcal C_{u,H}[C;0]
\overset{d}=\beta^2Q_{q,\tau},
\qquad
G_{u,H}(0)
\overset{d}=\log(\beta^2)+L_{q,\tau}.
\tag{3.4}
\]

The additive `log(beta^2)` cancels from adjacent differences. Because the
frozen frequency satisfies `u beta = 1/sqrt(log(e+n))`, the correction
depends on `n` through `q`, on the clock through `tau = rho^2 H`, and not on
the regime-specific value of beta.

## 4. Reversibility And Adjacent Windows

The invariant density of (3.1) is

\[
\pi(z)=\frac{2}{\pi}(1+z^2)^{-2}.
\tag{4.1}
\]

The stationary probability current is zero, so the diffusion is reversible.
Conditional on the common endpoint `Z_0=z`, the past and future path
segments are independent. By reversibility, their additive functionals have
the same conditional law.

Let `L^-` and `L^+` denote (3.3) on the adjacent past and future unit-clock
windows of width `tau`. Under stationarity,

\[
E\{(L^+-L^-)^2\mid Z_0\}
=2\operatorname{Var}(L_{q,\tau}\mid Z_0).
\tag{4.2}
\]

Therefore

\[
E_\pi\{(L^+-L^-)^2\}
=2E_\pi\left[
\operatorname{Var}(L_{q,\tau}\mid Z_0)
\right].
\tag{4.3}
\]

This identity incorporates overlap geometry correctly: the two windows do
not overlap, but they share the latent endpoint and are conditionally
independent rather than marginally independent.

## 5. Finite-Window Transfer Function

From (2.3), (3.4), and (4.3), the stationary expected numerator rate is

\[
\frac{E[J_{u,H,T}]}{T-2H}
=\frac{3\rho^2}{\tau}
E_\pi\left[
\operatorname{Var}(L_{q,\tau}\mid Z_0)
\right].
\tag{5.1}
\]

The stationary expected denominator rate is

\[
\frac{E[D_{u,H,T}]}{T-2H}
=4E_\pi\left(1-\frac1{Q_{q,\tau}}\right).
\tag{5.2}
\]

Define the dimensionless finite-window transfer

\[
a(q,\tau)
=
\frac{
\dfrac3\tau
E_\pi[\operatorname{Var}(L_{q,\tau}\mid Z_0)]
}{
4E_\pi(1-Q_{q,\tau}^{-1})
}.
\tag{5.3}
\]

By construction,

\[
E\left[
J_{u,H,T}
-\rho^2a(q,\rho^2H)D_{u,H,T}
\right]=0
\tag{5.4}
\]

under stationary initialization, apart from discrete ECF estimation
remainders already addressed by the analytical `h_log` subtraction.

Equation (5.4), not an expectation of a ratio, defines the corrected
estimating equation.

The observed-data prototype will solve

\[
\widehat I_{\log}
=r^2a(q,r^2H)\widehat D
\tag{5.5}
\]

over a predeclared positive rho domain. The old estimator is the special
case `a = 1`.

## 6. Vanishing-Window Limit

Let

\[
\ell_s=\log(1+Z_s^2).
\]

Its unit-clock martingale coefficient is

\[
\gamma(Z_s)=\frac{2Z_s}{\sqrt{1+Z_s^2}},
\qquad
\gamma^2(Z_s)=4\left(1-\frac1{1+Z_s^2}\right).
\tag{6.1}
\]

For small `tau`, the log of either local additive functional has the leading
conditional expansion

\[
L_{q,\tau}
=\ell_0+
\int_0^\tau\left(1-\frac{s}{\tau}\right)
\gamma(Z_0)dW_s+O_p(\tau).
\tag{6.2}
\]

Hence

\[
E_\pi[\operatorname{Var}(L_{q,\tau}\mid Z_0)]
=\frac{\tau}{3}E_\pi[\gamma^2(Z_0)]+o(\tau).
\tag{6.3}
\]

Also,

\[
4E_\pi(1-Q_{q,\tau}^{-1})
\longrightarrow
4E_\pi\left(1-\frac1{1+Z_0^2}\right)
=E_\pi[\gamma^2(Z_0)].
\tag{6.4}
\]

Substitution into (5.3) gives

\[
\boxed{\lim_{\tau\downarrow0}a(q,\tau)=1.}
\tag{6.5}
\]

Thus the corrected equation reduces to the current square-root ratio as the
physical window vanishes.

## 7. Why High Clock Is Attenuated

The relevant smoothing argument is

\[
\tau=\rho^2H.
\]

At fixed physical `H`, a larger rho makes the latent process traverse more
unit-clock time within every local ECF window. The additive functional in
(3.2) then averages over more latent evolution. Adjacent transformed window
averages fluctuate less than the pointwise `log(C)` process whose quadratic
variation appears in the old inversion.

For the frozen `n=10000000` design:

\[
\tau_{\rm low}\simeq0.0128,\quad
\tau_{\rm case2}\simeq0.0332,\quad
\tau_{\rm high}\simeq0.0765.
\]

The correction must therefore be negligible first in the low-clock regime
and largest in the high-clock regime. This ordering is a prediction of the
time-change model, not a calibration to the observed errors.

The numerical phase must verify whether `a(q,tau) < 1` over the relevant
domain. The derivation does not assume global monotonicity without checking
it.

## 8. Numerical Evaluation Strategy

Use the Lamperti state

\[
X_s=\operatorname{asinh}(Z_s),
\]

which solves

\[
dX_s=-\frac32\tanh(X_s)ds+dW_s.
\tag{8.1}
\]

Its diffusion coefficient is constant and its invariant density is
proportional to `sech(x)^3`. Exact stationary draws are obtained from

\[
Z_0=T_3/\sqrt3,\qquad X_0=\operatorname{asinh}(Z_0).
\]

For every candidate `tau`:

1. draw stationary `X_0`;
2. simulate two conditionally independent paths from the same `X_0`;
3. compute their two `L_{q,tau}` values and `Q_{q,tau}` values;
4. estimate
   \[
   E_\pi[\operatorname{Var}(L\mid Z_0)]
   =\frac12E[(L^{(1)}-L^{(2)})^2];
   \]
5. estimate the denominator expectation by averaging both replicate
   values of `1 - 1/Q`;
6. evaluate (5.3).

Common random numbers must be used across the predeclared rho grid.
Successive time-step resolutions must use a nested Brownian construction.
The implementation must report:

- numerator and denominator expectation estimates;
- Monte Carlo standard errors;
- discretization resolution;
- successive-resolution differences;
- transfer `a(q,tau)`;
- corrected map `rho^2 a(q,rho^2 H)`;
- map monotonicity;
- runtime and memory.

The numerator integrand has finite moments under (4.1) because logarithms
replace the polynomial activity tail, but it is not bounded. Monte Carlo
standard errors and independent replication are therefore mandatory.

Deterministic PDE or quadrature evaluation may later replace Monte Carlo,
but it is not required for the first controlled prototype.

## 9. Assumptions

The finite-window estimating equation requires:

1. the model dynamics for `V` are correctly specified;
2. `V` is stationary, or the nonstationary boundary contribution is
   asymptotically negligible relative to the observation horizon;
3. the diffusion is reversible under (4.1);
4. `alpha`, `sigma2`, and beta are fixed at their Level-4 oracle values;
5. local ECF sampling noise is removed to first order by the existing
   stable-adjusted `h_log` correction;
6. local empirical CFs remain statistically resolved;
7. the finite-window Riemann approximation is accurate at the chosen
   observation mesh;
8. the corrected map in (5.5) is one-to-one on the predeclared rho domain.

Stationarity is a substantive model assumption. The current direct
simulation starts at `V_0=0`, not at stationarity. At `T=50`, this creates a
transient contribution not represented by (5.4). It must be measured as a
separate approximation error; it must not be hidden by a burn-window switch.

## 10. What Is Exact And What Is Approximate

Exact:

- time change `V_t = Z_{rho^2 t}` in law;
- invariant density and reversibility;
- finite-window local ECF functional (2.1) in the fixed-window infill limit;
- adjacent-window conditional variance identity (4.2);
- stationary transfer definition (5.3);
- vanishing-window limit (6.5).

Approximate in the first prototype:

- replacing the discrete rolling statistic by (2.3);
- the first-order ECF sampling-noise subtraction;
- stationary averaging for paths initialized at zero;
- numerical simulation of the semigroup functional;
- interpolation and inversion of a finite rho grid.

The correction is an unbiased stationary estimating equation, not a
pathwise identity. Consequently, it can remove systematic finite-window
attenuation but cannot remove finite-horizon path variability.

## 11. Acceptance Criteria Before Level-4 Validation

The correction is implementable only if:

- independent numerical replications agree within reported uncertainty;
- successive time-step resolutions agree;
- the transfer is close to one as `tau` approaches zero;
- the corrected map is one-to-one over the frozen rho grid;
- interpolation error is negligible relative to Monte Carlo error;
- roots are interior and stable across numerical resolutions;
- no true path or true rho enters the correction.

If these conditions fail, the correction must stop at a negative numerical
decision rather than advancing to the Level-4 simulation gate.
