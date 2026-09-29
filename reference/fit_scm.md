# Fit a bias-corrected synthetic control

The Abadie, Diamond and Hainmueller synthetic control, matched on an
arbitrary predictor matrix rather than on pre-period outcomes alone,
with the Abadie-L'Hour / Ben-Michael-Feller-Rothstein bias correction
applied on top. This is the estimator Orlova, Rury and Wiltshire (2023)
use as their preferred specification, via the Stata package `allsynth`
(Wiltshire 2022);
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
covers the synthetic-DiD and plain-DiD comparisons they report alongside
it.

## Usage

``` r
fit_scm(
  Y,
  N0,
  T0,
  X = NULL,
  outcome_periods = NULL,
  bias_correct = TRUE,
  bc_ridge = 1e-04,
  bc_donors = c("auto", "positive", "all"),
  v_method = c("equal", "mspe", "synth"),
  v_starts = 8L,
  v_maxit = 1000L,
  penalty = 0
)
```

## Arguments

- Y:

  Outcome matrix, units x periods, controls first and the treated unit
  last – the layout
  [`synthdid::panel.matrices()`](https://rdrr.io/pkg/synthdid/man/panel.matrices.html)
  returns.

- N0:

  Number of control units (the first `N0` rows of `Y`).

- T0:

  Number of pre-treatment periods (the first `T0` columns).

- X:

  Optional units x predictors matrix of time-invariant predictors, rows
  in the same order as `Y`. Columns that are constant across units are
  dropped with a warning: they cannot discriminate between donors.

- outcome_periods:

  Integer indices of pre-treatment columns of `Y` to include as
  predictors, on top of `X`. `NULL` (the default) uses every
  pre-treatment period. Pass a subset to reproduce a specification that
  matches on selected pre-treatment years.

- bias_correct:

  Apply the bias correction. `TRUE` by default.

- bc_ridge:

  Ridge penalty for the bias-correction regressions, on standardised
  predictors. Small and positive; only binds when the positive-weight
  donors are few.

- bc_donors:

  Which donors the bias-correction regressions are fitted on.
  `"positive"` uses only donors carrying positive weight, as `allsynth`
  does; `"all"` uses the whole pool, as augmented SCM does; `"auto"`
  (the default) uses the positive-weight donors when there are
  comfortably more of them than predictors and falls back to the whole
  pool otherwise, reporting when it does. The fallback matters: a sparse
  synthetic control can easily name fewer donors than there are
  predictors, and the regression is then rank-deficient exactly when the
  correction is being asked to do the most work.

- v_method:

  How the diagonal of `V` is chosen. `"equal"` (default) weights
  standardised predictors equally. `"mspe"` searches for the `V`
  minimising pre-period MSPE, from `v_starts` deterministic starting
  points. `"synth"` does the `"mspe"` search AND asks the `Synth`
  package for its own solution, then keeps whichever fits the pre-period
  better – see the note below on why that is not redundant. `"synth"`
  needs the suggested `Synth` package and falls back to `"mspe"` with a
  warning if it is missing.

- v_starts:

  Restarts for the `"mspe"` search, the first from equal weights and the
  rest scattered. One start is not enough – see the note in the code.
  Ignored when `v_method = "equal"`.

- v_maxit:

  Nelder-Mead iterations per restart.

- penalty:

  Abadie-L'Hour (2021) penalty `lambda`, non-negative. Zero (the
  default) is the classic synthetic control. A positive value adds a
  charge for each donor's own distance from the treated unit, which
  makes the weights **unique** and favours donors that resemble the
  treated unit individually rather than only on average. Its main use
  here is breaking ties when the pre-period is fitted exactly: any small
  positive value then picks, among the exact fits, the one built from
  the nearest donors. Larger values trade pre-period fit for locality,
  and as `penalty` grows the synthetic control becomes the single
  nearest neighbour. The scale is that of the standardised predictors.

## Value

A list with `att` (mean corrected gap over the post-treatment periods),
`gap` (the corrected gap in every period, named by the columns of `Y`),
`gap_raw` (before correction), `pre_rmspe`, `post_rmspe`, `rmspe_ratio`
(`post_rmspe / pre_rmspe`), the same three computed on the uncorrected
gaps with a `_raw` suffix – plus `rmspe_ratio_hybrid` (corrected post
over uncorrected pre) – see the code for which to hand
[`rmspe_pvalue()`](https://james-hodges.github.io/sdidlab/reference/rmspe_pvalue.md)
– `omega` (donor weights, named by the control rownames), `v` (the
predictor weights used), `n_predictors`, `bias_correct`, `N0` and `T0`,
plus two flags. `exact_match` is `TRUE` when the predictors are matched
exactly: the weights are then generally not unique and the estimate
depends on which optimal weights are returned, so it warns (class
`sdidlab_exact_fit`) unless a `penalty` is set. `exact_fit` is `TRUE`
when the whole pre-period is reproduced exactly, which also leaves the
RMSPE ratios undefined (`NA`).

## Details

The weight solve is checked against an exact quadratic-programming
solution in the test suite: at convergence it reaches the same
objective, which is at least as low as `Synth`'s interior-point solve.

It is deliberately separate from
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
rather than another `estimator =` branch there.
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
delegates to `synthdid`, whose `sc_estimate()` matches on pre-period
outcomes only and has no bias correction; bolting predictor matching
onto it would change the meaning of the existing `"sc"` option and put
every stored run's numbers in play.

## What the two steps do

**Matching.** Donor weights `w` sit on the simplex and minimise \\(X_1 -
X_0 w)' V (X_1 - X_0 w)\\, where the columns of `X` are the predictors
and `V` is diagonal. Predictors are standardised by their cross-unit
standard deviation first, so `V = I` means "weight every predictor
equally after standardisation" – the special case of ADH that does not
require the nested optimisation. `v_method = "mspe"` instead searches
over diagonal `V` for the one minimising pre-period MSPE of the outcome,
which is what `Synth` does by default.

**Neither search dominates, which is why `"synth"` runs both.** `Synth`
optimises `V` with `optimx` from a single equal-weight start and solves
for `w` with an interior-point QP; this package uses several
deterministic starts and an exact active-set solve (see
[`sc_weights_simplex()`](https://james-hodges.github.io/sdidlab/reference/sc_weights_simplex.md)).
The `V` surface is rough, so the two searches can land in *different*
local minima, and either can find the better one. The objective is the
same and is directly comparable, so taking the better of the two is free
accuracy and the only defensible default when the surface is this rough.

When several convex combinations of donors fit equally well the
minimiser is not unique, and different solvers can return very different
comparison groups with the same loss. `penalty` resolves such ties;
without it, report `omega` rather than assuming it.

**Bias correction.** Even a good synthetic control leaves a residual
discrepancy in the predictors, and that discrepancy biases the estimate.
The correction regresses the outcome on the predictors *among the donor
units*, one regression per period, and subtracts the model's prediction
of the discrepancy. A ridge penalty keeps it defined when the positive-
weight donors are fewer than the predictors.

Note the diagnostic consequence, which is a feature and not a bug: **in
any period whose outcome is itself a predictor, the corrected gap is
mechanically zero**, because the period's regression fits that column
exactly. Seeing a row of exact zeros in the gap path is confirmation
that the correction ran, not evidence of a perfect match. It is what
produces the flat segments at matched years in bias-corrected SCM gap
figures.

## References

Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic Control
Methods for Comparative Case Studies. *Journal of the American
Statistical Association* 105(490), 493-505.

Abadie, A. and L'Hour, J. (2021). A Penalized Synthetic Control
Estimator for Disaggregated Data. *Journal of the American Statistical
Association* 116(536), 1817-1834.

Ben-Michael, E., Feller, A. and Rothstein, J. (2021). The Augmented
Synthetic Control Method. *Journal of the American Statistical
Association* 116(536), 1789-1803.

## See also

[`rmspe_pvalue()`](https://james-hodges.github.io/sdidlab/reference/rmspe_pvalue.md)
for the permutation inference this feeds,
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
for the synthetic-DiD and DiD comparisons.

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(11 * 8), nrow = 11,
            dimnames = list(c(paste0("d", 1:10), "treated"), 2010:2017))
Y[11, 7:8] <- Y[11, 7:8] + 3
fit_scm(Y, N0 = 10, T0 = 6)$att
#> Bias correction fitted on all 10 donors, not the 3 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> [1] 2.394834
```
