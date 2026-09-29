# A specification's stable identifier

Twelve hex characters of a SHA-256 over the specification's canonical
JSON. Metadata (`id`, `label`, `note`, `tags`) and exclusion notes are
excluded, so relabelling a specification does not change its identity,
and any field still at its default is dropped before hashing, so adding
a new field with a default leaves every existing hash unchanged. Field
order never matters.

## Usage

``` r
spec_hash(spec)
```

## Arguments

- spec:

  An
  [`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md).

## Value

A 12-character string.

## Examples

``` r
a <- sdid_spec(design = list(treat_time = 1989), label = "first name")
b <- sdid_spec(label = "second name", design = list(treat_time = 1989L))
identical(spec_hash(a), spec_hash(b))
#> [1] TRUE
```
