# Default specification

Every field of a specification, at its default.
[`sdid_spec()`](https://james-hodges.github.io/sdidlab/reference/sdid_spec.md)
merges overrides onto this, and
[`spec_hash()`](https://james-hodges.github.io/sdidlab/reference/spec_hash.md)
drops fields still at their default before hashing, so adding a new
field with a default never changes an existing specification's identity.

## Usage

``` r
spec_defaults()
```

## Value

A named list.
