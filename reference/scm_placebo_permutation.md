# Permute treatment in space across every donor

Refits
[`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md)
once per unit, each time pretending that unit is the treated one and
that the rest of the pool are its donors. The result is the reference
distribution
[`rmspe_pvalue()`](https://james-hodges.github.io/sdidlab/reference/rmspe_pvalue.md)
ranks against, and the set of grey placebo paths behind an Abadie-style
gap figure.

## Usage

``` r
scm_placebo_permutation(Y, N0, T0, X = NULL, ...)
```

## Arguments

- Y:

  Outcome matrix, units x periods, controls first and the treated unit
  last.

- N0:

  Number of control units.

- T0:

  Number of pre-treatment periods.

- X:

  Optional units x predictors matrix, rows aligned with `Y`.

- ...:

  Passed to
  [`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md)
  – `outcome_periods`, `bias_correct`, `v_method`, `bc_donors` and
  friends. The same arguments are used for the treated fit and every
  placebo, which is the point.

## Value

A list with `results`, a tibble of one row per unit (`unit`, `treated`,
`att`, `pre_rmspe`, `post_rmspe`, `rmspe_ratio` and their `_raw`
counterparts computed before bias correction), `gaps`, a units x periods
matrix of the corrected gap paths, and `gaps_raw`, the same before
correction.

## Details

The treated unit is always excluded from the donor pool of a placebo
fit. If it were left in, a placebo could reproduce itself out of the
genuinely treated unit's post-treatment trajectory, which would import
the effect being tested into the null distribution.

Everything is deterministic: the simplex solve is Frank-Wolfe from a
fixed start and the bias correction is a regression, so no seed is
involved and repeated calls agree exactly. That is the difference
between this and the placebo variance in
[`fit_sdid()`](https://james-hodges.github.io/sdidlab/reference/fit_sdid.md),
which is simulated.

## See also

[`rmspe_pvalue()`](https://james-hodges.github.io/sdidlab/reference/rmspe_pvalue.md)
for the ranking that turns this into a p-value, and
[`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md)
for the fit it repeats.

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(6 * 8), nrow = 6,
            dimnames = list(c(paste0("d", 1:5), "treated"), 2010:2017))
Y[6, 7:8] <- Y[6, 7:8] + 3
scm_placebo_permutation(Y, N0 = 5, T0 = 6)$results
#> Bias correction fitted on all 4 donors, not the 2 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> Bias correction fitted on all 4 donors, not the 3 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> Bias correction fitted on all 4 donors, not the 3 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> Bias correction fitted on all 4 donors, not the 3 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> Bias correction fitted on all 4 donors, not the 3 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> Bias correction fitted on all 5 donors, not the 3 with positive weight.
#> ℹ 6 predictors need more donors than that to be identified.
#> # A tibble: 6 × 12
#>   unit    treated    att pre_rmspe post_rmspe rmspe_ratio pre_rmspe_raw
#>   <chr>   <lgl>    <dbl>     <dbl>      <dbl>       <dbl>         <dbl>
#> 1 d1      FALSE   -0.162     0.342      0.950        2.78         0.595
#> 2 d2      FALSE    0.190     0.437      0.774        1.77         0.783
#> 3 d3      FALSE    0.185     0.371      1.16         3.12         0.728
#> 4 d4      FALSE   -0.365     0.942      1.02         1.08         1.05 
#> 5 d5      FALSE   -0.297     0.566      1.52         2.69         0.916
#> 6 treated TRUE     3.23      1.06       3.23         3.04         1.15 
#> # ℹ 5 more variables: post_rmspe_raw <dbl>, rmspe_ratio_raw <dbl>,
#> #   rmspe_ratio_hybrid <dbl>, exact_fit <lgl>, exact_match <lgl>
```
