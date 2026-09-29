# Declare a long panel for sdidlab

Wraps a long data frame – one row per unit and period – with the names
of its unit, time and outcome columns, so that specifications can refer
to columns by name. Any other columns (further outcomes, covariates, a
display label) stay available: a specification's `data$outcome` and
`covariates$vars` name them.

## Usage

``` r
sdid_panel(data, unit, time, outcome = NULL, treatment = NULL, label = NULL)
```

## Arguments

- data:

  A data frame in long format.

- unit, time:

  Names of the unit and time columns.

- outcome:

  Name of the default outcome column, used when a specification does not
  name one. Optional.

- treatment:

  Name of a 0/1 treatment column. Optional. When given, a specification
  that leaves `design$treated_unit` or `design$treat_time` unset takes
  them from it: the one unit that is ever treated, and its first treated
  period.

- label:

  Name of a column holding a display name for each unit, used in plots
  and tables. Optional.

## Value

An `sdid_panel` object: a list with `data` (the data frame, with the
unit column coerced to character) and the column names.

## Details

Units are stored as character, so numeric identifiers and names behave
the same way everywhere, and times must be numeric (years, or any
ordered numeric period index).

## Examples

``` r
data("california_prop99", package = "synthdid")
panel <- sdid_panel(california_prop99, unit = "State", time = "Year",
                    outcome = "PacksPerCapita", treatment = "treated")
panel
#> 
#> ── sdid_panel 
#> • 39 units (State) x 31 periods (Year, 1970-2000)
#> • outcome PacksPerCapita
#> • treatment treated
#> • columns PacksPerCapita and treated
```
