# What each donor contributes to the estimate and to the standard error

Sets each donor's share of the counterfactual (its unit weight) against
its share of the placebo variance behind the standard error – the draws
in which it stood in for the treated unit. A donor far above the
diagonal supplies more uncertainty than estimate.

## Usage

``` r
donor_influence(run)
```

## Arguments

- run:

  An `sdid_run` fitted with `keep = "all"`.

## Value

A tibble with one row per donor: `unit`, `label`, `weight`,
`weight_share`, `n_placebo`, `n_trimmed`, `att_mean`, `abs_att_max`,
`rmse_ratio` and `se_share`, largest `se_share` first.

## Details

Dropping a donor because it inflates the standard error is selection on
the result; this table is for understanding an interval, not for
choosing a pool.
