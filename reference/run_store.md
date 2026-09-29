# Save, load and index fitted runs

A run store is a directory of fitted runs, one folder per run, grouped
by specification hash. Each folder holds the specification as YAML, the
fitted object and placebo draws, the pool and balance records, and a
one-row summary that `refresh_runs_index()` collects into an index. A
run reloaded with `load_run()` can be re-plotted without refitting.

## Usage

``` r
save_run(run, store)

load_run(path)

refresh_runs_index(store, write = TRUE)

prune_runs(store, keep_latest_per_spec = 3L, dry_run = FALSE)
```

## Arguments

- run:

  An `sdid_run` from
  [`fit_spec()`](https://james-hodges.github.io/sdidlab/reference/fit_spec.md).

- store:

  Directory of the run store.

- path:

  A single run's directory (the `path` column of the index).

- write:

  Cache the index as `index.csv` in the store.

- keep_latest_per_spec:

  How many runs to keep for each specification.

- dry_run:

  Report what would be deleted without deleting it.

## Value

`save_run()` returns the run with its `path` set; `load_run()` an
`sdid_run`; `refresh_runs_index()` a tibble with one row per run;
`prune_runs()` the rows it deleted (or would delete).

## Examples

``` r
# \donttest{
store <- tempfile("runs")
run <- fit_spec(sdid_spec(inference = list(replications = 20)), prop99_panel(),
                store = store, seed = 1)
idx <- refresh_runs_index(store)
idx[, c("spec_hash", "att", "se")]
#> # A tibble: 1 × 3
#>   spec_hash      att    se
#>   <chr>        <dbl> <dbl>
#> 1 13393a31bcc0 -15.6  12.3
again <- load_run(idx$path[1])
all.equal(again$fit$att, run$fit$att)
#> [1] TRUE
# }
```
