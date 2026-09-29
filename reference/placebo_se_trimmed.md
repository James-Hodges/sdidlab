# Placebo standard error with poorly-fitting placebos trimmed

`synthdid`'s placebo variance permutes the control units, treats the
last one as if it were treated, re-estimates, and takes the spread of
those estimates. In a large heterogeneous donor pool most controls look
nothing like the treated unit, so many placebo draws are fitted badly
and produce large spurious "effects". Those draws fatten the tails of
the reference distribution and inflate the standard error.

## Usage

``` r
placebo_se_trimmed(est, replications = 200L, trim_k = 5, return_draws = FALSE)
```

## Arguments

- est:

  A `synthdid_estimate` object.

- replications:

  Number of placebo permutations.

- trim_k:

  Keep placebos whose pre-period trend RMSE is at most `trim_k` times
  the treated unit's. `Inf` disables trimming.

- return_draws:

  Also return the per-draw estimates. Consumes no extra randomness, so a
  run with and without it gives the same `se`.

## Value

A list: `se`, `n_kept`, `n_total`, `treated_rmse`,
`placebo_rmse_median`, and – when `return_draws` is `TRUE` – `draws`, a
data frame of `draw`, `att`, `rmse`, `keep` and `placebo_unit` (the
control that stood in for the treated unit).

## Details

The symptom is measurable with in-time placebos, where the true effect
is zero by construction: in large heterogeneous pools `|t|` comes out
far more concentrated near zero than the `|N(0, 1)|` a calibrated
interval implies, while small homogeneous pools show no such
compression.

The fix follows the spirit of Abadie, Diamond and Hainmueller (2010):
restrict the placebo reference distribution to units the estimator can
actually fit, by discarding placebos whose pre-treatment fit is worse
than `trim_k` times the treated unit's. A placebo that cannot be fitted
says nothing about how variable a well-fitted unit is – it is evidence
of a bad match, not of uncertainty.

**This trims the inference distribution, not the donor pool.** The
counterfactual is still built from every donor, and the criterion is
pre-treatment fit, evaluated without reference to any unit's estimated
effect, so it cannot be used to select donors that produce a preferred
estimate.

Fit is measured the way
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
measures it: RMSE of the treated- minus-synthetic gap after removing its
pre-period mean, on the covariate-adjusted outcome (see
[`adjusted_outcome()`](https://james-hodges.github.io/sdidlab/reference/adjusted_outcome.md)).

## Relationship to synthdid

Each placebo is refitted exactly as `synthdid:::placebo_se()` refits it:
`synthdid_estimate()` on the permuted panel, with the original fit's
options and its unit weights renormalised over the placebo controls as
the starting point. Because an `sc_estimate()` or `did_estimate()`
stores its fixed weights and options on the estimate, this reproduces
those estimators too. The draws consume the RNG exactly as `synthdid`
does, so under the same seed `trim_k = Inf` returns `synthdid`'s
`sqrt(vcov(est, method = "placebo"))` to machine precision.

## References

Abadie, A., Diamond, A. and Hainmueller, J. (2010). Synthetic Control
Methods for Comparative Case Studies. *Journal of the American
Statistical Association* 105(490), 493-505.

## Examples

``` r
# \donttest{
data("california_prop99", package = "synthdid")
setup <- synthdid::panel.matrices(california_prop99)
est <- synthdid::synthdid_estimate(setup$Y, setup$N0, setup$T0)
set.seed(1)
placebo_se_trimmed(est, replications = 50, trim_k = 5)$se
#> [1] 11.11442
# }
```
