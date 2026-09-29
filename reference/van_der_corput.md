# A deterministic low-discrepancy point in \[0, 1)

The van der Corput sequence: write `i` in the given base and mirror its
digits about the radix point. Successive points fill the interval far
more evenly than independent draws, which is what makes it a good source
of restart positions for an optimiser – and, unlike draws, it consumes
no randomness and returns the same points on every machine.

## Usage

``` r
van_der_corput(i, base = 2L)
```

## Arguments

- i:

  Index, a non-negative integer.

- base:

  Radix, an integer of at least 2. Using a different base per coordinate
  is what keeps a multi-dimensional start from lying on a diagonal.

## Value

A single number in `[0, 1)`.

## Examples

``` r
vapply(1:4, van_der_corput, numeric(1), base = 2)
#> [1] 0.500 0.250 0.750 0.125
```
