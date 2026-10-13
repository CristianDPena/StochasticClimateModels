# Benchmark protocol

All competitors see exactly the same observed residual vector and sampling
interval as the fixed estimator. Neither sees V, jump increments, or simulation
truth. Benchmark settings are fixed before production. They run on every
ordinary main/infill/horizon dataset, and their full outputs and warnings are
retained. They are not fed back into the final estimator.

## Classical Hill with a fixed intermediate rank

For descending absolute residuals A, k=floor(sqrt(n)):
alpha_hat = [mean(log(A[1:k]/A[k+1]))]^-1.
The fitted Pareto tail constant is K=(k/n)*A[k+1]^alpha_hat.
The stable scale is [K/(2*c_alpha*dt)]^(1/alpha_hat), where
c_alpha=Gamma(alpha)*sin(pi*alpha/2)/pi.

This is the classical Hill statistic with a transparent deterministic rank rule,
not a claim that sqrt(n) is optimal or universally standard. It compares adaptive
threshold selection against a fixed rule. It is a tail-only comparator, not a
four-parameter estimator. Reference: Hill (1975), doi:10.1214/aos/1176343247.

## Generalized Pareto peaks-over-threshold maximum likelihood

Use CRAN `evir` 1.7-4, `gpd(..., method="ml")`, at the empirical 99th percentile
(default R type-7 quantile). Retain the package's optimizer/default initial values.
For numerical units, divide absolute residuals and the threshold by their median
before fitting, then multiply the fitted GPD scale back. This deterministic scale
change is likelihood-equivariant, not truth-dependent tuning. Record convergence;
no alternative optimizer or threshold is selected after failures.

For positive fitted shape xi and scale s, alpha_hat=1/xi. The asymptotic GPD tail
constant is K=p_u*(s/xi)^(1/xi), with p_u the observed exceedance fraction. Apply
the same stable-tail conversion above. This scale conversion is model-specific
postprocessing, not a claim that evir itself estimates a stable jump scale.
A fitted alpha outside (0,2) is retained and marked outside the stable domain;
its scale estimate is undefined. No clipping, winsorization, or hidden refitting.

Source: https://CRAN.R-project.org/package=evir and its gpd reference manual.
The 99% threshold is prespecified, not theoretically optimal. Since residuals are
dependent and include a continuous component, the POT approximation is not exact.
Do not use package iid Hessian standard errors as valid model confidence intervals;
uncertainty is assessed across independent simulated trajectories.

## What is deliberately not compared

Pure iid stable-law quantile/ECF/ML fits assume that the entire residual has a
stable law, which this model does not. OU estimators assume constant diffusion
and/or observed latent states. Neither is an appropriate four-parameter comparator
here. Oracle inputs in separate diagnostics are labeled stronger-information
experiments and are never treated as fair observed-data competitors.

Both retained benchmarks estimate alpha and sigma2 only. Comparisons can support
claims about tail recovery at these thresholds. They cannot establish superiority
of the complete alpha/sigma2/sigma1/rho pipeline over all established methods.
