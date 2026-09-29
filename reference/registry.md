# Register derived outcomes and named donor pools

Extension points for data that a plain column cannot express. A
specification names an outcome in `data$outcome` and a pool in
`pool$base`; when the name is not a column (or `"all"`), sdidlab looks
it up here.

## Usage

``` r
register_outcome(name, fn, label = name, units = NA_character_)

register_pool(name, fn)

registered_outcomes()

registered_pools()
```

## Arguments

- name:

  The name specifications will use.

- fn:

  The function, as described above.

- label, units:

  For outcomes: how axes and tables should describe it.

## Value

The name, invisibly.

## Details

- An **outcome** function takes `(panel, spec)` and returns a data frame
  with columns `unit`, `time` and `y` – for example a share computed
  from two columns, or a series built from another table.

- A **pool** function takes `(panel, spec)` and returns either a
  character vector of donor units, or a data frame with columns `unit`,
  `included` (logical) and optionally `reason`, which is how a pool
  records *why* a candidate was left out.

Both receive the whole specification, so a spec's free-form `options`
section can parameterise them. Registrations last for the R session; put
them in a package's `.onLoad()` or at the top of a script.

## Examples

``` r
register_outcome("packs_log", function(panel, spec) {
  d <- panel$data
  data.frame(unit = d[[panel$unit]], time = d[[panel$time]],
             y = log(d$PacksPerCapita))
}, label = "Log packs per capita")

register_pool("west", function(panel, spec) {
  c("Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico",
    "Oregon", "Utah", "Washington", "Wyoming")
})

registered_outcomes()
#> [1] "packs_log"
registered_pools()
#> [1] "west"
```
