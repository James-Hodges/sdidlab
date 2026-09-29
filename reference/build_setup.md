# Build the estimator inputs for a specification

Resolves a specification against a panel: the outcome series, the
treated unit and its timing, the donor pool, the periods (window,
dropped periods, blackout), a balanced panel, and the covariate array.
The result is what
[`fit_spec()`](https://james-hodges.github.io/sdidlab/reference/fit_spec.md)
hands to the estimator, and is useful on its own for checking what a
specification actually selects.

## Usage

``` r
build_setup(spec, panel)
```

## Arguments

- spec:

  An
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md).

- panel:

  An
  [`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md).

## Value

An `sdid_setup`: `Y`, `N0`, `T0` (as
[`synthdid::panel.matrices()`](https://rdrr.io/pkg/synthdid/man/panel.matrices.html)
returns them), `X` (covariate array or `NULL`), `pool`, `balance`,
`periods`, `blackout`, `treated_unit`, `treated_label`, `treat_time` and
`spec`.

## Details

Every exclusion is recorded rather than silent: `pool` lists each
candidate donor with whether it was included and why not, and `balance`
lists donors dropped for missing the outcome in a kept period.

## Examples

``` r
setup <- build_setup(sdid_spec(), prop99_panel())
dim(setup$Y)
#> [1] 39 31
c(N0 = setup$N0, T0 = setup$T0)
#> N0 T0 
#> 38 19 
```
