# The outcome a synthdid estimate was actually fitted on

`synthdid` estimates covariate coefficients `beta` jointly with the
weights and computes the effect on `Y - X beta`. Any diagnostic of the
fit – the pre-period gap, the trend RMSE, a gap-path plot – has to be
computed on that adjusted outcome, or it describes a fit the estimator
never made. Without covariates this is `Y` itself.

## Usage

``` r
adjusted_outcome(est)
```

## Arguments

- est:

  A `synthdid_estimate` object.

## Value

A units x periods matrix.

## Examples

``` r
# \donttest{
data("california_prop99", package = "synthdid")
setup <- synthdid::panel.matrices(california_prop99)
est <- synthdid::synthdid_estimate(setup$Y, setup$N0, setup$T0)
identical(adjusted_outcome(est), setup$Y)
#> [1] TRUE
# }
```
