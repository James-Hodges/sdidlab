# Fit one synthetic-DiD specification and return its full result set

Wraps
[`synthdid::synthdid_estimate()`](https://rdrr.io/pkg/synthdid/man/synthdid_estimate.html)
(or `sc_estimate()` / `did_estimate()`) with the inference and
diagnostics a single-treated-unit design needs on every specification.

## Usage

``` r
fit_sdid(
  Y,
  N0,
  T0,
  X = NULL,
  replications = 200L,
  conf_level = 0.95,
  trim_k = 5,
  se_method = c("placebo", "none"),
  estimator = c("sdid", "sc", "did"),
  keep = c("none", "est", "all")
)
```

## Arguments

- Y:

  Outcome matrix, units x periods, controls first and the treated unit
  last – exactly what
  [`synthdid::panel.matrices()`](https://rdrr.io/pkg/synthdid/man/panel.matrices.html)
  returns.

- N0:

  Number of control units (the first `N0` rows of `Y`).

- T0:

  Number of pre-treatment periods (the first `T0` columns).

- X:

  Optional unit x period x covariate array, handed to `synthdid`, which
  estimates the covariate coefficients jointly with the weights. To
  residualise on untreated cells first instead, use
  [`project_out_covariates()`](https://james-hodges.github.io/sdidlab/reference/project_out_covariates.md)
  and pass the result as `Y`.

- replications:

  Placebo replications. 200 is `synthdid`'s default.

- conf_level:

  Confidence level for the interval.

- trim_k:

  Passed to
  [`placebo_se_trimmed()`](https://james-hodges.github.io/sdidlab/reference/placebo_se_trimmed.md).
  Placebos fitting worse than `trim_k` times the treated unit's
  pre-period trend RMSE are dropped from the reference distribution.
  `Inf` reproduces `synthdid`'s untrimmed placebo SE.

- se_method:

  `"placebo"` (default) or `"none"`. Placebo inference refits the model
  once per replication and dominates the runtime, so sweeps that only
  need point estimates should pass `"none"`, which returns `NA` for
  `se`, `t`, `ci_lo` and `ci_hi`.

- estimator:

  `"sdid"` (default), `"sc"` for synthetic control, or `"did"` for plain
  difference-in-differences, all as `synthdid` defines them. The placebo
  distribution is refitted with the same estimator.

- keep:

  How much to return beyond the summary. `"none"` (default) returns the
  summary fields only. `"est"` adds the `synthdid_estimate` object,
  which every plot needs. `"all"` adds the individual placebo draws.
  Neither consumes extra randomness, so the estimate and standard error
  are unchanged.

## Value

A list with `att`, `se`, `t`, `ci_lo`, `ci_hi`, `trend_rmse`,
`level_gap`, `N0`, `T0`, `n_periods`, `omega` (unit weights, named by
the control rownames of `Y`), `lambda` (time weights), `estimator`, and
the placebo trimming record `trim_k`, `placebos_kept`, `placebos_total`.
With `keep = "est"`, also `est`; with `keep = "all"`, also
`placebo_draws`, `treated_rmse`, `placebo_rmse_median`, `unit_ids` and
`periods`.

## Details

Three things it standardises:

- **Placebo inference.** With a single treated unit the jackknife
  variance is undefined – it would require dropping the treated unit –
  so the placebo variance is the default. With `trim_k = Inf` it is
  exactly `synthdid`'s own `vcov(est, method = "placebo")`; see
  [`placebo_se_trimmed()`](https://james-hodges.github.io/sdidlab/reference/placebo_se_trimmed.md)
  for what trimming changes.

- **The right pre-period fit statistic.** SDID carries unit fixed
  effects and therefore tolerates a *constant* level gap between treated
  and synthetic control; absorbing it is what the DiD half does. Raw
  RMSE against the synthetic path penalises exactly what the estimator
  is designed to ignore. `trend_rmse` is the RMSE of the treated-minus-
  synthetic gap **after removing its pre-period mean**, which is the
  parallel-trends assumption; `level_gap` reports the absorbed offset
  separately. Both are computed on the outcome the estimator actually
  fits: when covariates are passed as `X`, that is `Y - X beta`, not
  `Y`.

- **A confidence interval.** When every estimate is null the bound is
  the result, and a bound cannot be read off a p-value.

Placebo standard errors are simulated: the same specification moves by a
few percent between runs depending on the RNG state. Set a seed before
calling if results must be reproducible.

## References

Arkhangelsky, D., Athey, S., Hirshberg, D. A., Imbens, G. W. and Wager,
S. (2021). Synthetic Difference-in-Differences. *American Economic
Review* 111(12), 4088-4118.

## Examples

``` r
# \donttest{
data("california_prop99", package = "synthdid")
setup <- synthdid::panel.matrices(california_prop99)
set.seed(1)
r <- fit_sdid(setup$Y, setup$N0, setup$T0, replications = 50)
c(att = r$att, lo = r$ci_lo, hi = r$ci_hi)
#>        att         lo         hi 
#> -15.603828 -37.387689   6.180034 
# }
```
