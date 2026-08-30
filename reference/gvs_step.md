# Greedy variable selection: test one candidate interaction

Fits the nested propensity score model with one candidate interaction
term added and returns a likelihood-ratio test of that term against the
current model.
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
calls this once per candidate at every round of greedy variable
selection (GVS), the forward step of the RAILS algorithm; the candidate
with the largest significant increment per degree of freedom is the one
added.

## Usage

``` r
gvs_step(
  x,
  terms_current,
  loglik,
  cells_ref,
  cells_np,
  m_ref,
  maxit = 1000,
  tol = 0.001
)
```

## Arguments

- x:

  Term or terms to add, e.g. `"agegroup:sex:income"`. A vector is tested
  as a single increment, which is how
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
  keeps a selected higher-order term together with any lower-order
  components it needs.

- terms_current:

  Character vector of terms already in the current model.

- loglik:

  Log-likelihood of the current model.

- cells_ref:

  Aggregated reference-sample cell counts (column `weight`).

- cells_np:

  Aggregated non-probability cell counts (column `weight`).

- m_ref:

  Model matrix for the current model on `cells_ref`, used only for its
  column count, to compute degrees of freedom.

- maxit:

  Maximum Newton-Raphson iterations.

- tol:

  Convergence tolerance on the squared Newton step. Looser than
  [`nps_weights()`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md)
  because only the log-likelihood is needed, not the weights.

## Value

A length-5 numeric vector `c(loglik, LRT, p-value, LRT/df, df)`, or
`rep(NA, 5)` if any cell count is zero or the solver fails.

## See also

[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
which runs the GVS loop, and
[`nps_weights()`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md)
for the model being tested.

## Examples

``` r
pop <- rails_simulate(20000, seed = 42)
vars <- c("agegroup", "sex", "income")
np  <- rails_cells(pop[pop$aou == 1, ], vars)
ref <- rails_cells(pop[pop$s == 1, ], vars, weights = 1 / pop$ps[pop$s == 1])

# The current model: main effects only.
terms_now <- vars
f <- stats::formula(paste("~", paste(terms_now, collapse = "+")))
d <- nps_weights(f, np, ref)
m_ref <- Matrix::sparse.model.matrix(f, ref)
loglik <- sum(as.numeric(Matrix::crossprod(
    Matrix::sparse.model.matrix(f, np), np$weight)) * attr(d, "theta")) -
  sum(ref$weight * log1p(exp(as.numeric(m_ref %*% attr(d, "theta")))))

# Test one candidate interaction against it.
gvs_step("sex:income", terms_now, loglik, ref, np, m_ref)
#> [1] -8.698596e+03  1.637203e+01  9.512354e-04  5.457343e+00  3.000000e+00
```
