# Summarise fitted runs as a table

One row per run: what was estimated, the estimate, its interval, and the
fit and inference records. `NULL` entries (specifications that failed in
[`run_specs()`](https://james-hodges.github.io/sdidlab/reference/run_specs.md))
are skipped.

## Usage

``` r
sdid_results_table(runs, digits = 3)
```

## Arguments

- runs:

  A list of `sdid_run`s.

- digits:

  Rounding for the estimate columns.

## Value

A data frame.

## Examples

``` r
# \donttest{
runs <- run_specs(list(sdid_spec(label = "SDID", inference = list(replications = 20))),
                  prop99_panel(), progress = FALSE)
sdid_results_table(runs)
#>   label estimator        outcome    treated treat_time    window N0 T0   ATT
#> 1  SDID      sdid PacksPerCapita California       1989 1970-2000 38 19 -15.6
#>     SE     t  ci_lo ci_hi trend_rmse placebos    spec_hash
#> 1 10.6 -1.48 -36.33 5.122       1.73    19/20 13393a31bcc0
# }
```
