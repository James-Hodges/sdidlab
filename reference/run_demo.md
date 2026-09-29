# Launch the interactive demo

A Shiny app for building a specification, fitting it, inspecting every
diagnostic plot, comparing specifications, and downloading them as YAML
plus an R script that reproduces the fits. It runs on the built-in
Proposition 99 panel or an uploaded CSV.

## Usage

``` r
run_demo(...)
```

## Arguments

- ...:

  Passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), e.g.
  `port` or `launch.browser`.

## Value

Called for its side effect of running the app.

## Examples

``` r
if (interactive()) run_demo()
```
