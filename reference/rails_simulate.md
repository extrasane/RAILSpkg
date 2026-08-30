# Simulate a population with a biased volunteer sample

Generates the synthetic population used by the RAILS simulation studies:
five categorical covariates, a binary outcome, a probability sample, and
a non-probability ("volunteer") sample whose inclusion depends on the
same covariates. Because the true population is known, the bias of any
weighting method can be read off directly. This is what the examples and
vignettes run on – real biobank microdata cannot be redistributed.

## Usage

``` r
rails_simulate(
  n,
  alpha = c(-3.5, 0.05, 0.5, 0.2, 0.4, 0.6, 0.01, 0.1, 0.3, 0.1, 0.2, 0.3),
  beta = c(-2.2, 0.03, 0.2, 0.1, 0.2, 0.3, 0.005, 0.05, 0.1, 0.05, 0.1, 0.15),
  gamma = c(-4.8, 0.08, 0.6, 0.3, 0.6, 0.9, 0.02, 0.2, 0.4, 0.15, 0.3, 0.45),
  seed = NULL
)
```

## Arguments

- n:

  Population size.

- alpha:

  Length-12 coefficient vector for the outcome model.

- beta:

  Length-12 coefficient vector for probability-sample inclusion.

- gamma:

  Length-12 coefficient vector for volunteer-sample inclusion. Making
  `gamma` differ from `beta` is what creates the selection bias RAILS is
  meant to remove.

- seed:

  Optional integer seed. `NULL` leaves the RNG state alone.

## Value

A data frame with one row per population member:

- `id`:

  Row identifier.

- `agegroup`:

  `young`, `middle`, `older`.

- `sex`:

  `male`, `female`.

- `income`:

  `p0-10`, `p10-50`, `p50-90`, `p90-100` – cut at those percentiles of a
  truncated-Pareto draw.

- `homeown`:

  `rent`, `own`.

- `edu`:

  `<HS`, `HS`, `Some college`, `College+`.

- `y`,`py`:

  Binary outcome and its true probability.

- `s`,`ps`,`fpc`:

  Probability-sample indicator, its inclusion probability, and the
  population size.

- `aou`,`paou`:

  Volunteer-sample indicator and its inclusion probability.

Design weights for the probability sample are `1 / ps`.

## These are synthetic, not calibrated

The variable names and level labels are there so worked examples read
like the application rather than like `agegroup:sex`. The marginal
distributions behind them are arbitrary and resemble no real population:
`edu` in particular puts only 0.5 percent in its top category. Nothing
here should be read as a claim about any real group.

## Difference from the paper code

Income is drawn from a truncated Pareto by inverse-CDF sampling rather
than via `distributionsrd::rtruncpareto()`. Same distribution, one less
dependency; streams drawn under a given seed therefore differ from the
simulation archive.

## Examples

``` r
pop <- rails_simulate(1000, seed = 1)
head(pop[, c("agegroup", "sex", "income", "homeown", "edu", "y")])
#>   agegroup    sex  income homeown          edu y
#> 1    young   male  p50-90    rent Some college 0
#> 2   middle   male  p10-50    rent           HS 0
#> 3    young   male  p10-50    rent Some college 0
#> 4    older female  p10-50    rent           HS 0
#> 5   middle female  p50-90    rent          <HS 0
#> 6    young female p90-100    rent Some college 0
table(pop$aou, pop$s)
#>    
#>       0   1
#>   0 625 197
#>   1 123  55
mean(pop$y)                        # population prevalence
#> [1] 0.218
mean(pop$y[pop$aou == 1])          # what the volunteer sample would report
#> [1] 0.2865169
```
