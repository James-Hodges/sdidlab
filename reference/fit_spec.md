# Fit one specification end to end

Specification in, run out: builds the estimator inputs with
[`build_setup()`](https://james-hodges.github.io/sdidlab/reference/build_setup.md),
fits, runs inference, and – given a store – saves the whole run so it
can be reloaded and re-plotted without refitting.

## Usage

``` r
fit_spec(
  spec,
  panel,
  store = NULL,
  keep = c("all", "est", "none"),
  seed = NULL
)
```

## Arguments

- spec:

  An
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md).

- panel:

  An
  [`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md).

- store:

  Directory of a run store to save into, or `NULL`.

- keep:

  How much of the fit to retain: `"all"` (default) keeps the fitted
  object and every placebo draw, which every plot needs; `"est"` keeps
  only the fitted object; `"none"` keeps the summary alone.

- seed:

  `NULL` to inherit the caller's RNG stream, or an integer seed applied
  to this fit only.

## Value

An `sdid_run`: `spec`, `spec_hash`, `run_id`, `fit`, `pool`, `balance`,
`periods`, `timings` and `meta`.

## Reproducibility of the standard error

Placebo standard errors are simulated, so they depend on the state of
the RNG stream when the fit begins. `seed = NULL` inherits the caller's
stream: reproducible for a whole sweep run after one
[`set.seed()`](https://rdrr.io/r/base/Random.html), but a spec refitted
on its own gets a different SE. Passing `seed = seed_for(spec)` makes a
spec reproducible in isolation. The caller's own stream is left
untouched either way.

## Examples

``` r
# \donttest{
spec <- sdid_spec(label = "Prop 99", inference = list(replications = 50))
run <- fit_spec(spec, prop99_panel(), seed = 1)
run
#> 
#> ── sdid_run Prop 99 
#> • ATT -15.6 95% CI [-37.4, 6.18]
#> • SE 11.1 (t = -1.40)
#> • fit sdid, N0 = 38, T0 = 19, pre-period RMSE 1.73
#> • placebos 47/50 kept at trim_k 5
#> • treated California from 1989
#> • spec bf7b1f5e74c4 run 20260929T053936906476-bf7b1f5e
# }
```
