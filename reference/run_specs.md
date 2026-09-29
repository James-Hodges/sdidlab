# Fit a list of specifications

Fits each specification with
[`fit_spec()`](https://james-hodges.github.io/sdidlab/reference/fit_spec.md),
in order. A specification that cannot be fitted (too few donors, a
too-short window) is reported and skipped, and the rest carry on.

## Usage

``` r
run_specs(
  specs,
  panel,
  store = NULL,
  keep = c("all", "est", "none"),
  seed_each = TRUE,
  progress = TRUE
)
```

## Arguments

- specs:

  A list of
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md)s,
  for example from
  [`read_specs()`](https://james-hodges.github.io/sdidlab/reference/spec_yaml.md).

- panel:

  An
  [`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md).

- store:

  Directory of a run store to save into, or `NULL`.

- keep:

  How much of the fit to retain: `"all"` (default) keeps the fitted
  object and every placebo draw, which every plot needs; `"est"` keeps
  only the fitted object; `"none"` keeps the summary alone.

- seed_each:

  Seed each fit with
  [`seed_for()`](https://james-hodges.github.io/sdidlab/reference/seed_for.md),
  so every specification's standard error is reproducible on its own.
  `FALSE` inherits one stream across the whole sweep instead, so a
  result depends on its position.

- progress:

  Print one line per fitted specification.

## Value

A list of `sdid_run`s, with `NULL` for specifications that failed.

## Examples

``` r
# \donttest{
specs <- list(
  sdid_spec(label = "SDID", inference = list(replications = 30)),
  sdid_spec(label = "SC", estimator = "sc", inference = list(replications = 30))
)
runs <- run_specs(specs, prop99_panel())
#> ✔ 1/2 SDID: ATT -15.6, SE 9.55
#> ✔ 2/2 SC: ATT -19.62, SE 8.25
sdid_results_table(runs)
#>   label estimator        outcome    treated treat_time    window N0 T0    ATT
#> 1  SDID      sdid PacksPerCapita California       1989 1970-2000 38 19 -15.60
#> 2    SC        sc PacksPerCapita California       1989 1970-2000 38 19 -19.62
#>     SE     t  ci_lo  ci_hi trend_rmse placebos    spec_hash
#> 1 9.55 -1.63 -34.32  3.115       1.73    30/30 8dc336c466ef
#> 2 8.25 -2.38 -35.78 -3.459       1.66    29/30 3ed8e77dff0b
# }
```
