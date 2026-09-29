# Residualise the outcome on covariates using untreated cells only

The covariate adjustment of Kranz (2022): regress the outcome on the
covariates with unit and time fixed effects, using only the
**untreated** unit-periods – every control-unit cell plus the treated
unit's pre-period – then subtract the covariates' fitted contribution
from the whole outcome matrix. Fitting on untreated cells keeps the
treatment effect out of the covariate coefficients. The residualised
matrix is then passed to
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md)
as `Y` with no `X`.

## Usage

``` r
project_out_covariates(Y, X, N0, T0)
```

## Arguments

- Y:

  Outcome matrix, units x periods, controls first and the treated unit
  last.

- X:

  Units x periods x covariates array aligned with `Y` (a units x periods
  matrix is accepted for a single covariate).

- N0:

  Number of control units.

- T0:

  Number of pre-treatment periods.

## Value

`Y` with the covariates' estimated contribution removed, with the
coefficients attached as the `"beta"` attribute. A covariate that is
collinear with the fixed effects gets an `NA` coefficient and is left
out of the adjustment.

## Details

This is the alternative to handing `X` to `synthdid`, which estimates
the coefficients jointly with the synthetic-control weights. The two can
give materially different answers when the covariates are strongly
related to the outcome; comparing them is a useful robustness check.

The coefficients are estimated once, on the real treatment assignment.
Placebo inference on the residualised matrix therefore holds them fixed
rather than re-estimating them for each placebo.

## References

Kranz, S. (2022). Synthetic Difference-in-Differences with Time-Varying
Covariates. Working paper.

## Examples

``` r
# \donttest{
data("california_prop99", package = "synthdid")
setup <- synthdid::panel.matrices(california_prop99)
set.seed(1)
X <- setup$Y * 0.01 + rnorm(length(setup$Y))
Yres <- project_out_covariates(setup$Y, X, setup$N0, setup$T0)
attr(Yres, "beta")
#>        x1 
#> 0.7019451 
# }
```
