# Simplex weights minimising a V-weighted predictor discrepancy

Solves the synthetic-control weight problem exactly: \$\$\min_w (Z_1 -
Z_0' w)' V (Z_1 - Z_0' w) + \lambda \sum_j w_j (Z_1 - Z\_{0j})' V (Z_1 -
Z\_{0j})\$\$ over the unit simplex (`w >= 0`, `sum(w) = 1`), where
`V = diag(v)` and `lambda` is the Abadie-L'Hour `penalty` (zero by
default, the classic synthetic control).

## Usage

``` r
sc_weights_simplex(Z0, Z1, v = rep(1, ncol(Z0)), penalty = 0)
```

## Arguments

- Z0:

  Donors x predictors matrix, already standardised.

- Z1:

  The treated unit's predictor vector.

- v:

  Non-negative predictor weights (the diagonal of `V`).

- penalty:

  The Abadie-L'Hour penalty `lambda`, non-negative. The scale is that of
  the standardised predictors. See
  [`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md)
  for when to use it.

## Value

A numeric vector of donor weights, non-negative and summing to 1.

## Method

A primal active-set method – Lawson and Hanson's non-negative least
squares algorithm, adapted to the simplex. It keeps a working set of
donors, solves the equality-constrained problem on that set exactly
(through a pseudo-inverse, so a singular set is handled), steps back to
feasibility when a weight would turn negative, and adds the donor whose
gradient most violates optimality. It stops at the exact KKT point:
every donor with positive weight has the same gradient and no donor has
a lower one. Starting from the nearest donor, it typically needs about
two steps per donor in the final support.

Two alternatives were tried and rejected. A first-order (Frank-Wolfe)
solve needs over 12,000 iterations on the highly collinear pre-period
outcomes of the Proposition 99 panel, and at any practical iteration cap
returns weights well short of the optimum. A dual QP solver needs a
ridge to make the singular quadratic term invertible, which with a
penalty makes the problem nearly linear and the solve numerically
unreliable.

## Uniqueness

When the treated unit lies inside the donors' convex hull, many weight
vectors can fit the predictors equally well, and they imply different
counterfactuals. Which of them is returned then depends on the order in
which donors enter the working set – deterministic, but arbitrary.
[`fit_scm()`](https://james-hodges.github.io/sdidlab/reference/fit_scm.md)
warns when the fit is exact. A positive `penalty` makes the solution
unique (Abadie and L'Hour 2021).

## References

Lawson, C. L. and Hanson, R. J. (1974). *Solving Least Squares
Problems*. Prentice-Hall, chapter 23.

## Examples

``` r
Z0 <- matrix(c(0, 1, 2, 0, 1, 2), nrow = 3)
sc_weights_simplex(Z0, c(1, 1), v = c(1, 1))
#> [1] 0 1 0
```
