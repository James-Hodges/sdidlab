# Resolve a specification's donor pool

Every candidate donor with whether it is included and, if not, why. The
candidates are the panel's other units (`pool$base = "all"`), the units
listed in `pool$units`, or the result of a registered pool; then
`drop_units` are excluded and `add_units` forced in.

## Usage

``` r
resolve_pool(spec, panel, treated)
```

## Arguments

- spec:

  An
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md).

- panel:

  An
  [`sdid_panel()`](https://james-hodges.github.io/sdidlab/reference/sdid_panel.md).

- treated:

  The treated unit, which is never a candidate.

## Value

A tibble with `unit`, `label`, `included` and `reason`.

## Examples

``` r
pool <- resolve_pool(sdid_spec(pool = list(drop_units = "Nevada")),
                     prop99_panel(), treated = "California")
pool[!pool$included, ]
#> # A tibble: 1 × 4
#>   unit   label  included reason         
#>   <chr>  <chr>  <lgl>    <chr>          
#> 1 Nevada Nevada FALSE    dropped by hand
```
