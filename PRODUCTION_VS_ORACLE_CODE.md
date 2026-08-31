# Production Versus Oracle Code

## Production Mode

The authoritative controller fixes:

```text
--tail-input=estimated
--chunked-likelihood
--groups-per-chunk=16
```

at `work/run_overarching_validation_batch.ps1:186-198`.

Under these settings:

1. alpha and sigma2 are estimated from `R`;
2. beta uses those estimates in the conditional endpoint likelihood;
3. rho uses those estimates and beta in the rolling log-ECF statistic;
4. sigma1 is beta divided by corrected rho;
5. the realized latent state is never supplied to an estimator.

## Selectable Tail Modes

Defined at `work/run_conditional_cf_endpoint_full_pipeline.R:12-19`:

| Mode/flag | alpha passed downstream | sigma2 passed downstream | Production? |
|---|---|---|---|
| `--tail-input=estimated` | Hill estimate | stable-tail estimate | **Yes** |
| `--tail-input=true-alpha` | true alpha overwrites Hill alpha | the already computed estimated-alpha sigma2 remains | No |
| `--tail-input=true-sigma2` | Hill estimate | true sigma2 overwrites estimate | No |
| `--tail-input=oracle` | true alpha | true sigma2 | No |
| legacy `--oracle-tail` | selects `oracle` by default | selects `oracle` by default | No |

The overwrite logic is at
`work/run_conditional_cf_endpoint_full_pipeline.R:164-177`. In particular,
`true-alpha` does not recompute sigma2 using true alpha; it only replaces the
alpha field. This branch is diagnostic and was not selected by the frozen
controller.

## Truth Used to Simulate and Score

The validation runner necessarily uses the design truth to generate a path:

- configuration: `work/run_conditional_cf_endpoint_full_pipeline.R:30-64`
- simulation: lines `155-161`

It later uses truth to calculate error columns:

- endpoint-runner result columns: lines `249-281`
- final corrected errors: `work/run_finite_window_corrected_full_pipeline.R:84-99`
- tail error decomposition and sigma1 cancellation diagnostics:
  lines `131-177`

These uses are downstream evaluation only. They do not feed back into a
production estimate when `tail_input_mode == "estimated"`.

## Realized Latent V

`simulate_residual_euler` constructs `v` internally to generate the simulated
residual:

- `work/occupation_transition_clock_prototype.R:28-38`

The authoritative caller sets `keep_v=FALSE`:

- `work/run_conditional_cf_endpoint_full_pipeline.R:155-161`

The simulator therefore returns only `R` and `dt`. No estimator receives the
realized `v`, and the outer runner has no `V` object.

This is different from saying the simulator does not have a latent state. It
does; the estimator does not access that realized state.

## Generic Simulated Paths Used by Production

### Endpoint kernel

The endpoint cache contains candidate-rho simulations of generic unit-model
`Q=1+V^2` transitions and within-block paths:

- generator: `work/diagnose_exact_cf_endpoint_hmm.R:36-161`
- cache read: `work/run_conditional_cf_endpoint_full_pipeline.R:111-137`
- likelihood use:
  `work/diagnose_exact_cf_endpoint_hmm.R:956-990,1023-1160`

These are numerical quadrature/Monte Carlo objects for a likelihood. They are
not the latent path corresponding to the observed residuals.

### Finite-window transfer

The correction CSV contains expectations from paired generic stationary
unit-clock simulations:

- generator: `finite_window_corrected_rho_estimator.R:99-287`
- profile generation settings:
  `work/build_overarching_correction_profiles.R:17-59`
- read and selection:
  `work/run_finite_window_corrected_full_pipeline.R:33-52`

Again, this is a generic expectation table, not pathwise oracle information.

## Oracle Definitions Present but Unreachable

| File/location | Oracle content | Reachable from production entry point? |
|---|---|---:|
| `rho_level4_oracle_residual_estimator.R:354-501` | `rho_level4_oracle_diagnostics` accepts oracle local variance, true rho, and true sigma1 | No |
| `work/diagnose_markov_additive_endpoint_hmm.R:247-264` | `block_truth(V,...)` | No |
| `work/diagnose_markov_additive_endpoint_hmm.R:836-911` | known-endpoint/oracle-Qbar likelihoods | No |
| `work/diagnose_markov_additive_endpoint_hmm.R:1679-2191` | standalone diagnostic program with oracle and observed modes | No; suppressed by `MARKOV_ENDPOINT_LIBRARY_ONLY=1` |
| `work/diagnose_exact_cf_endpoint_hmm.R:1212-1278` | standalone simulation test with supplied truth | No; `sys.nframe()!=0` when sourced |
| `work/occupation_transition_clock_prototype.R:979-1010` | oracle activity-noise mode | Function definitions are loaded, but not called |

## Mixed Files

Several frozen files contain both production functions and old diagnostic
functions. File inclusion is therefore not evidence of oracle use. The
reachable production call chain is listed in `CALL_GRAPH.md`, and it excludes
all of the oracle routines above.

The most misleading filename is
`rho_level4_oracle_residual_estimator.R`. Its production function
`rho_level4_estimate` uses only `dR`, scalar `dt`, alpha, sigma2, beta, and
numerical settings. The word `oracle` refers to the residual-level validation
stage and to a separately named diagnostic function, not to the inputs of the
production estimator.

## Production Reachability Conclusion

- True alpha/sigma2 branches: present, but disabled by the frozen controller.
- True rho/sigma1 error calculations: present after estimation, evaluation
  only.
- Realized latent `V_t`: generated internally for simulation, discarded, and
  not used by estimation.
- Generic simulated paths and transfer expectations: used as frozen numerical
  likelihood/correction objects.
- Oracle diagnostics: defined in mixed source files, not called by the
  authoritative production entry point.

