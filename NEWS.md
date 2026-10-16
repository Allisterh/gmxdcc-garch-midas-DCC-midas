# gmxdcc 0.1.2

* New optional analytic gradients. **Default behavior is unchanged**: with
  `analytic_grad = FALSE` (the default) every estimate, log-likelihood,
  standard error and output field other than the two new `spec` entries
  is identical to 0.1.1.
  - New exported functions `gm_score()` / `gm_expected_info()`: analytic
    per-observation score and expected information matrix of the Gaussian
    GARCH-MIDAS quasi-likelihood (any number of covariates including
    realized volatility, any number of regimes, Beta or Almon weights;
    `distribution = "std"` is not supported and raises an error).
  - New exported functions `dcc_midas_score()` / `dcc_expected_info()`:
    the same for the DCC-MIDAS correlation step, for two assets with
    exogenous covariates only (`rc = FALSE`), with or without breaks. The
    recursion starts where `dcc_midas_loglik()` starts it (`burn_in`
    option respected); rows before the start are zero.
  - New argument `analytic_grad = FALSE` (last position) of
    `fit_garch_midas()` and `fit_dcc()`. `TRUE` passes the analytic score
    to `maxLik` as `grad`. In `fit_dcc()` it acts on the first step
    (`univariate = "garch_midas"`) and on the second step (bivariate
    DCC-MIDAS, `rc = FALSE`); other correlation models, more than two
    assets or `rc = TRUE` fall back to the numerical gradient in the second
    step with a one-time message. With `TRUE`, `gradientObs` and the
    numerical Hessian (a difference of the analytic gradient) also change,
    so robust standard errors differ slightly from the default fit.
  - The choice is recorded in `fit$spec$analytic_grad`; `fit_dcc()` also
    stores `fit$spec$analytic_grad_used` (which steps actually used it).
  - Internal: `run_multistart()` gains `grad_fn = NULL` (last position);
    with `NULL` the `maxLik` call is the same as before.
* No change to the C++ kernels, the likelihoods, the start-value draw, or
  the constraints.

# gmxdcc 0.1.1

* New `burn_in` argument of `fit_dcc()` controlling where the correlation
  recursion starts:
  - `"dccmidas"` (default) keeps the convention of the dccmidas package and
    the legacy script — **results are identical to 0.1.0**. For DCC-MIDAS
    the first `K_c` days are placeholders; for cDCC/aDCC/DECO the first day
    (or the first `K_c` days if `K_c` is supplied).
  - `"exact"` starts as soon as every input of the long-run correlation is
    genuinely available: after `N_c + K_c` days when the realized-correlation
    term is present (its first `N_c` daily correlations are zero-filled, so
    under `"dccmidas"` the early long-run correlations are partly built from
    zeros), after one day otherwise (monthly covariates are already complete
    after calendar alignment; `K_c` counts months there, not days).
* New output field `burn_in` (integer): the number of leading in-sample days
  whose `R_t` / `H_t` / `R_t_bar` slices are placeholders — identity
  correlation, `H_t` equal to the covariance of the *standardized* residuals
  (not on the return scale), zero likelihood contribution. Drop them before
  computing conditional betas or forecast losses:
  `fit$H_t[, , -seq_len(fit$burn_in)]`. Arrays keep their full length, as
  before. The chosen option is recorded in `fit$spec$burn_in`, and
  `summary()` prints a burn-in line.
* Implemented on the R side only; the C++ kernels are unchanged.

# gmxdcc 0.1.0 (updates since the initial GitHub release)

* Fixed the routing of optimizer start-value injection in `fit_dcc()`:
  `control$start_matrix` now targets the second (correlation) step only.
  Previously the same `control` list was forwarded to the univariate stage,
  so a correlation-sized start matrix leaked into `fit_garch_midas()` and
  raised a column-count error whenever `univariate = "garch_midas"`.
* New: `univariate_args$start_matrix` injects first-step optimizer starts —
  a single matrix applied to every asset, or a list with one matrix per
  asset (columns in the canonical `gm_par_index()` order). Together with
  `control$start_matrix` this enables cheap warm restarts across rolling /
  recursive re-estimation windows: seed each window's optimizer with the
  previous window's estimates instead of a fresh random multistart.
* Injected starts that sit on (or, after the 6-decimal rounding of
  `coef()`, marginally outside) an inequality constraint boundary are now
  nudged minimally into the strict interior before optimization.
  `constrOptim`'s log-barrier requires a strictly interior start, so warm
  restarts from a near-boundary optimum previously failed with
  "initial value in 'vmmin' is not finite".

# gmxdcc 0.1.0

First packaged release, refactored from the research script
`integrated noskew GARCH-DCC-MIDAS fit.R` (frozen under `dev/legacy/`).

## New interface

* `fit_garch_midas()` replaces `garchmidas_noskew_fit()` and its wrapper:
  the long-run component now accepts an arbitrary named list of
  low-frequency covariates (`X`), each with its own lag order (`K`), plus
  an optional internally generated realized-volatility covariate (`rv`).
* `fit_dcc()` replaces `dccmidas_noskew_longfocus_fit()` and its wrapper,
  with `correlation = "cDCC" / "aDCC" / "DECO" / "DCC-MIDAS"`, an `rc`
  toggle, a `corr_X` covariate list and `K_c`.
* Structural breaks are arguments (`break_dates`, `corr_break_dates`)
  instead of `_sb` model-name suffixes; intercepts and slopes become
  piecewise-constant across the implied regimes.
* `seed` makes the random multi-start reproducible (the legacy script drew
  starting values without a seed, so fits were not reproducible).
* Optimizer bounds/start ranges are collected in documented tables
  (`gm_bounds_default`, `dcc_bounds_default`) and can be overridden via
  `control$bounds`.

## Behavior changes vs. the legacy script

* **Out-of-sample warm start.** Held-out conditional variances and
  correlations are obtained by filtering the full sample with frozen
  in-sample parameters. The legacy script restarted the recursions cold at
  the first out-of-sample day (short-run variance component reset to 1,
  DCC quasi-correlation reset to the identity) and, for cDCC/aDCC/DECO,
  computed the unconditional target from *out-of-sample* residuals
  (a look-ahead). The C++ kernels gained optional `S_in`/`N_in` arguments
  so all unconditional targets now come from in-sample residuals; the
  `N_c_res_for_oos` patch mechanism is gone.
* In-sample results are numerically identical to the legacy script; the
  regression-test suite verifies likelihoods and components at fixed
  parameters to 1e-12 and end-to-end fits from identical starting values.

## Bug fixes folded into the port

* `integrated_garchmidas_noskew_output_function()` referenced the global
  `univ_model` instead of its own `model` argument (worked only by global
  leakage); the new API has no model strings, removing the bug class.
* `stop(cat(...))` calls (which raised empty error messages) replaced by
  proper `stop()` messages via `assert_xts()` and validators.
* `align_sample()` (né `cut_func`) no longer errors when `period`/`macro`
  are `NULL`.
* MIDAS weights are computed internally (`midas_weights()`), numerically
  identical to `rumidas::beta_function()`/`exp_almon()` (verified by unit
  test); `rumidas` is no longer a hard dependency.
