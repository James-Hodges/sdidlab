# Read and write specifications as YAML

`write_spec_yaml()` writes one specification; `write_specs_yaml()`
writes a *sweep* – several specifications in one file, with the settings
they share hoisted into a `defaults:` block so the file shows exactly
what differs between them. `read_specs()` reads either form (or every
YAML file in a directory) and `read_spec_yaml()` reads a single-spec
file.

## Usage

``` r
write_spec_yaml(spec, path)

read_spec_yaml(path)

read_specs(path)

write_specs_yaml(specs, path, title = "sdidlab specifications")
```

## Arguments

- spec:

  An
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md).

- path:

  File to write, or `NULL` to return the text instead. For
  `read_specs()`, a file or a directory.

- specs:

  A list of
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md)s.

- title:

  Comment line at the top of a sweep file.

## Value

The writers return the path invisibly (or the YAML text when
`path = NULL`); `read_specs()` returns a list of specifications and
`read_spec_yaml()` a single one.

## Details

Writing then reading is an identity: the specifications that come back
have the same
[`spec_hash()`](https://james-hodges.github.io/sdidlab/reference/spec_hash.md)
as those written.

A sweep file looks like this:

    defaults:
      design: {treat_time: 1989}
      inference: {replications: 100}
    specs:
    - label: SDID
    - label: Synthetic control
      estimator: sc

## Examples

``` r
specs <- list(
  sdid_spec(label = "SDID", design = list(treat_time = 1989)),
  sdid_spec(label = "SC", estimator = "sc", design = list(treat_time = 1989))
)
cat(write_specs_yaml(specs, NULL, title = "Prop 99 estimators"))
#> # Prop 99 estimators
#> # Placebo standard errors depend on the RNG stream: fit these in file
#> # order under one seed, or with run_specs(seed_each = TRUE).
#> defaults:
#>   design:
#>     treat_time: 1989.0
#> specs:
#> - label: SDID
#> - label: SC
#>   estimator: sc

f <- tempfile(fileext = ".yml")
write_specs_yaml(specs, f)
back <- read_specs(f)
identical(vapply(back, spec_hash, ""), vapply(specs, spec_hash, ""))
#> [1] TRUE
```
