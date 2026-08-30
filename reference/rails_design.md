# Hand RAILS weights to the survey package

Wraps a fit's weights in a
[`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html)
so the usual survey machinery – `svymean()`, `svytotal()`, `svyglm()`,
`svyby()` – can be used on the weighted sample.

## Usage

``` r
rails_design(fit, data, ...)
```

## Arguments

- fit:

  A `rails_fit` from
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
  or a `rails_subgroup_fit` from
  [`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md).

- data:

  The data frame the weights belong to: the same non-probability sample
  passed to
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
  with whatever outcome columns you want to analyse. Must have as many
  rows as the fit has weights.

- ...:

  Passed to
  [`survey::svydesign()`](https://rdrr.io/pkg/survey/man/svydesign.html).

## Value

A `survey.design` object with the RAILS weights attached and a
`rails_weights` column added to its data.

## What this design does and does not carry

The weights are exact: point estimates from this object – `svymean()`,
`svytotal()`, `svyglm()`, `svyby()` – are the RAILS estimates.

Standard errors are a different matter. The survey package computes its
own design-based variance from the weights it is handed, treating them
as fixed. That is close to, but not the same estimator as,
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
with `type = "simplified"`, which is computed independently; in practice
the two agree to about four significant figures, differing by little
more than an n/(n-1) factor. Neither accounts for the weights having
been estimated.

So: use this design freely for point estimates and for model fitting,
and use
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
when you want inference on a mean that reflects how the weights were
built. Do not expect the two standard errors to match to the last digit
– they are different estimators of the same fixed-weight quantity.

## See also

[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
for corrected inference on a weighted mean.

## Examples

``` r
pop <- rails_simulate(20000, seed = 42)
aou <- pop[pop$aou == 1, ]
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps

fit <- rails(aou, ref, vars = c("agegroup", "sex", "income"),
             weights_ref = "dweight", verbose = FALSE)
des <- rails_design(fit, aou)
survey::svymean(~y, des)
#>      mean     SE
#> y 0.20402 0.0076
```
