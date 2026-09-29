# Build a specification

A specification is a complete, declarative description of one model: the
outcome, the treated unit and timing, the donor pool, covariates, the
estimator and the inference. Specifications hash to a stable identifier
([`spec_hash()`](https://james-hodges.github.io/sdidlab/reference/spec_hash.md)),
round-trip through YAML
([`write_spec_yaml()`](https://james-hodges.github.io/sdidlab/reference/spec_yaml.md)),
and are fitted by
[`fit_spec()`](https://james-hodges.github.io/sdidlab/reference/fit_spec.md).

## Usage

``` r
sdid_spec(...)
```

## Arguments

- ...:

  Named sections to override: `data`, `design`, `pool`, `covariates`,
  `estimator`, `scm`, `inference`, `options`, and the metadata fields
  `id`, `label`, `note`, `tags`.

## Value

An `sdid_spec`.

## Details

Pass any subset of sections; each merges onto
[`spec_defaults()`](https://james-hodges.github.io/sdidlab/reference/spec_defaults.md),
one level deep (two for `design$window`), so
`design = list(treat_time = 1989)` leaves the other design fields at
their defaults. Unknown fields are an error rather than silently
ignored.

## Examples

``` r
spec <- sdid_spec(
  label = "Prop 99, SDID",
  design = list(treat_time = 1989, window = list(start = 1970)),
  inference = list(replications = 100)
)
spec
#> 
#> ── sdid_spec Prop 99, SDID 
#> • hash ec87073104d3
#> • outcome (panel default)
#> • treated from the panel, treated from 1989
#> • window 1970-(panel end)
#> • pool all
#> • estimator sdid
#> • inference placebo, 100 reps, trim_k 5
spec_hash(spec)
#> [1] "ec87073104d3"
```
