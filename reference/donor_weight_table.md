# Donor weights, largest first

Donor weights, largest first

## Usage

``` r
donor_weight_table(run)
```

## Arguments

- run:

  An `sdid_run`.

## Value

A tibble of `unit`, `label`, `weight` and `cumulative` weight.

## Examples

``` r
# \donttest{
run <- fit_spec(sdid_spec(inference = list(se_method = "none")), prop99_panel())
head(donor_weight_table(run))
#> # A tibble: 6 × 4
#>   unit          label         weight cumulative
#>   <chr>         <chr>          <dbl>      <dbl>
#> 1 Nevada        Nevada        0.124       0.124
#> 2 New Hampshire New Hampshire 0.105       0.230
#> 3 Connecticut   Connecticut   0.0783      0.308
#> 4 Delaware      Delaware      0.0704      0.378
#> 5 Colorado      Colorado      0.0575      0.436
#> 6 Illinois      Illinois      0.0534      0.489
# }
```
