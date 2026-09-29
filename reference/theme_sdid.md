# The sdidlab ggplot2 theme

The sdidlab ggplot2 theme

## Usage

``` r
theme_sdid(base_size = 12)
```

## Arguments

- base_size:

  Base font size.

## Value

A ggplot2 theme.

## Examples

``` r
ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_sdid()
```
