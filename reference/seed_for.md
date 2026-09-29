# A specification's own seed

Folds a specification's hash into a base seed, so refitting one
specification alone reproduces the standard error a sweep reported.

## Usage

``` r
seed_for(spec, base = 20110901L)
```

## Arguments

- spec:

  An
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md).

- base:

  Base seed to fold the hash into.

## Value

A single integer.

## Examples

``` r
seed_for(sdid_spec())
#> [1] 141890334
```
