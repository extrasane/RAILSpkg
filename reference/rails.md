# RAILS: raking with assisted nested propensity score and LIFO selection

Combines a non-probability sample with a probability reference sample to
produce population-representative weights. The procedure is:

## Usage

``` r
rails(
  np,
  ref,
  vars,
  weights_ref = NULL,
  weights_np = NULL,
  pop_totals = NULL,
  aggregated = FALSE,
  start = 2,
  scope = 3,
  alpha = 0.05,
  nsiz = NULL,
  benchmarks = FALSE,
  fallback = FALSE,
  lifo = c("descending", "ascending"),
  keep_data = TRUE,
  verbose = TRUE
)
```

## Arguments

- np:

  The non-probability sample (e.g. a volunteer biobank): a data frame of
  individual records, or a cell table if `aggregated = TRUE`.

- ref:

  The probability reference sample, same form as `np`.

- vars:

  Character vector of covariate names shared by both samples. These are
  the main effects; interactions are generated from them.

- weights_ref:

  Design weights for `ref`: a numeric vector or the name of a column.
  Almost always required, since a probability sample's design weights
  are what make its margins represent the population. Ignored when
  `aggregated = TRUE`.

- weights_np:

  Optional weights for `np`. `NULL` counts each record as 1, which is
  the usual choice for a biobank. Ignored when `aggregated = TRUE`.

- pop_totals:

  Named numeric vector of population margins. `NULL` (the default)
  computes them from `ref` via
  [`rails_totals()`](https://extrasane.github.io/RAILSpkg/reference/rails_totals.md).
  Supply your own when the margins come from an external source (a
  census table, say) rather than from the reference sample itself.

- aggregated:

  Set `TRUE` when `np` and `ref` are already cell tables with a `weight`
  column. The row-to-cell map is then unavailable, so
  [`weights()`](https://rdrr.io/r/stats/weights.html) returns one weight
  per cell rather than per record, and
  [`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
  will not run. See the note below.

- start:

  The starting model: every term in it is included without being tested.
  Give an interaction order (`2`, the default, means all main effects
  and all two-way interactions of `vars`), a one-sided formula, or a
  character vector of terms such as `c("age", "sex", "age:sex")`.

- scope:

  The largest model selection may reach, in the same three forms. `3`
  (the default) means up to three-way interactions. Everything in
  `scope` but not in `start` is a candidate for the forward test, so
  `start = 2, scope = 3` offers the three-way terms, while
  `start = 1, scope = 3` offers the two-way *and* three-way terms.
  `scope` must contain `start`; when the two are equal there is nothing
  to select and the starting model is simply raked. Term components may
  be written in any order – `"sex:age"` matches `"age:sex"`.

- alpha:

  Significance threshold for the forward likelihood-ratio test.

- nsiz:

  Population total the weights are scaled to. `NULL` uses the total
  reference weight, which is the right choice when `ref` carries design
  weights that already sum to the population.

- benchmarks:

  Also compute the comparison methods used in the papers (unweighted,
  raking-only, propensity-only, propensity-plus-raking), placed in
  `fit$benchmarks` as per-record weights, on the same scale as
  [`weights()`](https://rdrr.io/r/stats/weights.html). Off by default:
  this is simulation-study output, not something you need in order to
  produce weights.

- fallback:

  When selection picks nothing, or the first stepwise step fails to
  converge, fall back to the raked starting model instead of returning
  `NA` weights. `FALSE` is the default so that `rails()` reproduces the
  published estimator exactly; `TRUE` is often more useful in applied
  work. Either way the outcome is recorded in `converged`. This has no
  effect when `scope` equals `start`: there the starting model is what
  was asked for, not a fallback, and is always used.

- lifo:

  Direction of the stepwise raking stage. `"descending"` (the default)
  is the LIFO strategy of Chen et al. (2025): start at the full selected
  model and pop the last-added term until raking converges.
  `"ascending"` instead climbs from the starting model and halts at the
  first non-convergent step, which is what the application code did; use
  it to reproduce published application results exactly. The two agree
  whenever raking convergence is monotone in the number of added terms;
  where it is not, `"descending"` keeps more terms.

- keep_data:

  Retain the `vars` columns of both samples inside the fit, so
  [`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
  can rebuild individual-level model matrices later. Costs memory
  proportional to the input; set `FALSE` if you only need weights.

- verbose:

  Report selection and stepwise progress via
  [`message()`](https://rdrr.io/r/base/message.html).

## Value

An object of class `rails_fit`: a list with `weights` (one per input
record, or per cell when `aggregated = TRUE`), `cells`, `weights_cell`
(cell totals), `terms_base`, `terms_scope`, `terms_selected`,
`terms_used`, `formula_used`, `theta`, `pop_totals`, `n_used`,
`converged`, `benchmarks`, and `call`. Use
[`weights()`](https://rdrr.io/r/stats/weights.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on it, and
[`rails_design()`](https://extrasane.github.io/RAILSpkg/reference/rails_design.md)
to hand it to survey.

## Details

1.  Fit the **nested propensity score** (NPS) model on the starting
    model `start`, giving base weights and a baseline log-likelihood.

2.  **Greedy variable selection** (GVS): a forward likelihood-ratio
    search over the terms in `scope` but not in `start` (via
    [`gvs_step()`](https://extrasane.github.io/RAILSpkg/reference/gvs_step.md))
    at significance level `alpha`, adding whichever significant
    candidate has the largest average increment in the log
    pseudo-likelihood.

3.  **LIFO stepwise raking**: rake the full selected model to the
    population margins. If raking does not converge, pop the last-added
    term, refit NPS, and rake again – last in, first out – until it
    converges or no selected term is left.

Those three stages are what the method is named for in Chen et al.
(2025): **R**aking with **A**ssisted Nested Propens**i**ty Score and
**L**IFO **S**election.

All computation runs on aggregated covariate cells, so cost scales with
the number of distinct covariate patterns rather than the sample size.

## Benchmark names

With `benchmarks = TRUE`, `fit$benchmarks` carries the comparison
methods. They correspond to the estimators tabulated in Chen et al.
(2025):

|                |             |                                                 |
|----------------|-------------|-------------------------------------------------|
| package        | paper       | description                                     |
| `d_unweighted` | `naive`     | the unweighted cohort mean                      |
| `d_cal1`       | `cal-1`     | raking to univariate margins from equal weights |
| `d_cal2`       | `cal-2`     | raking on main effects and all two-way terms    |
| `d_nps1`       | `nps-1`     | NPS weights, main effects only                  |
| `d_nps2`       | `nps-2`     | NPS weights through the starting model          |
| `d_nps1_rake`  | `nps-cal-1` | NPS main effects, then univariate raking        |
| `d_nps2_rake`  |             | NPS through the starting model, then raking     |

`cal-2` is the one the paper reports as failing to converge in over 99
percent of simulations; here it returns `NA` with a warning when raking
fails, so a non-convergent benchmark is visible rather than silently
absent.

## Reproducing the published estimator

The defaults – `start = 2`, `scope = 3`, `fallback = FALSE` – are
exactly the procedure of
[`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md),
so `rails()` and the deprecated shim give identical weights on identical
input. Changing `start` or `scope` departs from the paper deliberately;
`fallback = TRUE` departs from it only in the degenerate case where the
published code would have returned `NA`.

## Aggregated input

`aggregated = TRUE` accepts cell tables directly, which is what the
original
[`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
required. It is the right choice when the microdata cannot be
materialized in one frame. The cost is that the fit no longer knows
which record belongs to which cell, so per-record weights and the
variance estimator are unavailable.

## See also

[`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md)
to run within strata,
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
for standard errors,
[`rails_cells()`](https://extrasane.github.io/RAILSpkg/reference/rails_cells.md)
and
[`rails_totals()`](https://extrasane.github.io/RAILSpkg/reference/rails_totals.md)
to build inputs by hand.

## Examples

``` r
pop <- rails_simulate(20000, seed = 42)
aou <- pop[pop$aou == 1, ]
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps

fit <- rails(aou, ref, vars = c("agegroup", "sex", "income"),
             weights_ref = "dweight", verbose = FALSE)
fit
#> RAILS fit
#>   covariates   : agegroup, sex, income
#>   start model  : 6 term(s)
#>   scope        : 7 term(s), so 1 candidate(s)
#>   selected     : 1 term(s) at alpha = 0.05
#>                  agegroup:sex:income
#>   raked        : 1 of 1 selected term(s) survived the stepwise walk
#>   cells        : 24
#>   weights      : n = 3799, sum = 19,875.78, range 1.73 to 28.7

# Weighted prevalence of the outcome, against the population truth:
w <- weights(fit)
sum(w * aou$y) / sum(w)
#> [1] 0.2040237
mean(pop$y)
#> [1] 0.1952
```
