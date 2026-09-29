# California's Proposition 99 as an sdid_panel

The example panel from Abadie, Diamond and Hainmueller (2010) and
Arkhangelsky et al. (2021): per-capita cigarette sales in 39 US states,
1970-2000, with California treated from 1989. The data ship with
`synthdid`; this wraps them for use with
[`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md)
and
[`fit_spec()`](https://james-hodges.github.io/sdidlab/reference/fit_spec.md).

## Usage

``` r
prop99_panel()
```

## Value

An
[`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md).

## Examples

``` r
prop99_panel()
#> 
#> ── sdid_panel 
#> • 39 units (State) x 31 periods (Year, 1970-2000)
#> • outcome PacksPerCapita
#> • treatment treated
#> • columns PacksPerCapita and treated
```
