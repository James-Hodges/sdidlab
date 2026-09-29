# Ranked-RMSPE permutation p-value for a synthetic control

The inference Abadie, Diamond and Hainmueller propose and Orlova, Rury
and Wiltshire (2023) use: fit the same synthetic control with each donor
standing in as the treated unit, rank the treated unit's ratio of
post-treatment to pre-treatment RMSPE against that distribution, and
report the rank as a fraction of the units.

## Usage

``` r
rmspe_pvalue(
  ratio,
  att,
  treated,
  alternative = c("two.sided", "greater", "less")
)
```

## Arguments

- ratio:

  Numeric vector of RMSPE ratios, one per unit – the `rmspe_ratio`
  element of
  [`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md),
  for the treated unit and every in-space placebo. Named, ideally.

- att:

  Numeric vector of estimated effects in the same order, used only to
  sign the statistic for a one-sided test.

- treated:

  Index or name identifying the treated unit within `ratio`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

## Value

A list with `p_value`, `rank` (the treated unit's rank, 1 being the most
extreme), `n_units`, `statistic` (the treated unit's signed or unsigned
ratio), `floor` (`1 / n_units`), and `alternative`.

## Details

Dividing by the pre-period RMSPE is what makes the statistic comparable
across units: a placebo that never had a good synthetic control in the
first place has a large post-period gap for reasons that have nothing to
do with treatment, and the ratio nets that out. Poorly fitting placebos
are therefore **kept**, not trimmed – trimming them is what would make
the test anti-conservative.

## What one-sided means here

The RMSPE ratio is unsigned, so a two-sided test can rank it directly. A
one-sided test needs the sign back: the statistic becomes
`sign(att) * ratio`, so a placebo that moved hard in the *opposite*
direction to the hypothesised effect ranks below the treated unit rather
than beside it. Where the direction of any effect is genuinely known in
advance, this is a real power gain with a small number of donors; where
it is not, `"two.sided"` is the honest default and is what this function
does unless told otherwise.

## The floor, which matters with few donors

The smallest attainable p-value is `1/n_units`: even a treated unit that
out-ranks every placebo gets `1/n`, never zero. With 40 units that floor
is 0.025, and a p-value of 0.05 means "second largest of forty", not
"significant at five percent" in the usual sampling sense. Read these as
ranks, and do not report them to three decimals as if they were
continuous.

## See also

[`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md),
and
[`scm_placebo_permutation()`](https://james-hodges.github.io/sdidlab/reference/scm_placebo_permutation.md)
which produces the vectors this takes.

## Examples

``` r
rmspe_pvalue(c(wm = 6.2, a = 1.1, b = 2.0, c = 0.7),
             att = c(28.9, -2, 4, 1), treated = "wm")
#> $p_value
#> [1] 0.25
#> 
#> $rank
#> [1] 1
#> 
#> $n_units
#> [1] 4
#> 
#> $statistic
#> [1] 6.2
#> 
#> $floor
#> [1] 0.25
#> 
#> $alternative
#> [1] "two.sided"
#> 
```
