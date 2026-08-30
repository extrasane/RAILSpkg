# Nested propensity score weights

Estimates non-probability-sample weights by maximizing the
pseudo-likelihood of the **nested propensity score** (NPS) model via
Newton-Raphson on aggregated covariate cells. The returned weights are
scaled to sum to `nsiz`.

## Usage

``` r
nps_weights(
  formula,
  cells_np,
  cells_ref,
  nsiz = sum(cells_ref$weight),
  theta_init = NULL,
  maxit = 1000,
  tol = 1e-10
)
```

## Arguments

- formula:

  A one-sided model formula, e.g. `~ agegroup + sex + agegroup:sex`.

- cells_np:

  Aggregated non-probability cell counts, with a column `weight` giving
  the cell size.

- cells_ref:

  Aggregated reference-sample cell counts, with a column `weight`.

- nsiz:

  Target population total the weights are scaled to. Defaults to the
  total reference-sample weight, `sum(cells_ref$weight)`.

- theta_init:

  Optional warm-start coefficient vector (e.g. the previous step's
  solution, zero-padded for new columns). The pseudo-likelihood is
  globally concave, so the warm start only reduces the iteration count.

- maxit:

  Maximum Newton-Raphson iterations.

- tol:

  Convergence tolerance on the squared Newton step.

## Value

A numeric vector of propensity weights (one per cell), scaled to sum to
`nsiz`, with the converged coefficient vector attached as attribute
`"theta"`.

## Details

This is the NPS step of the RAILS algorithm: the pseudo-likelihood
merges the non-probability cohort with the probability sample and the
resulting inverse fitted propensities serve as base weights for the
calibration stage.

## See also

[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
which calls this at the start and at every step of the stepwise raking
walk.

## Examples

``` r
pop <- rails_simulate(2000, seed = 1)
np  <- rails_cells(pop[pop$aou == 1, ], c("agegroup", "sex"))
ref <- rails_cells(pop[pop$s == 1, ],   c("agegroup", "sex"))
w   <- nps_weights(~ agegroup + sex, np, ref)
round(w, 2)
#> [1]  36.01  29.71  25.45 118.27  33.05 244.52
#> attr(,"theta")
#> [1] -2.1701939  0.8562832  3.3126288  2.9970826
```
