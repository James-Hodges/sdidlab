# Donors carrying almost no weight

Splits the donors at a weight threshold – useful for reporting how
concentrated a counterfactual is, or for a robustness refit on the
material donors only (decided in advance, not after seeing the result).

## Usage

``` r
low_weight_units(run, threshold = 0.01)
```

## Arguments

- run:

  An `sdid_run`.

- threshold:

  Weight below which a donor counts as negligible.

## Value

A list: `drop` and `keep` (unit identifiers), their counts, the share of
total weight below the threshold, and the full `table`.
