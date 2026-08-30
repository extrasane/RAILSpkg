# Summarize a RAILS fit

Summarize a RAILS fit

## Usage

``` r
# S3 method for class 'rails_fit'
summary(object, ...)

# S3 method for class 'summary.rails_fit'
print(x, ...)
```

## Arguments

- object:

  A `rails_fit`.

- ...:

  Ignored.

- x:

  A `summary.rails_fit`.

## Value

A `summary.rails_fit` object: a list with the selection path, the terms
actually raked to, weight diagnostics from
[`weight_summary()`](https://extrasane.github.io/RAILSpkg/reference/weight_summary.md),
the design effect, and the effective sample size.

## Examples

``` r
pop <- rails_simulate(2000, seed = 3)
ref <- pop[pop$s == 1, ]; ref$dweight <- 1 / ref$ps
fit <- rails(pop[pop$aou == 1, ], ref, c("agegroup", "sex"),
             weights_ref = "dweight", verbose = FALSE)
summary(fit)
#> RAILS fit summary
#> 
#> Covariates: agegroup, sex
#> Model     : 3 term(s) in the starting model, 0 candidate(s) in scope, alpha 0.05
#> 
#> Selected terms (0):
#>   none
#> 
#> Raked to 3 term(s); full walk completed.
#> 
#> Weight diagnostics:
#>                sum                var non_positive_ratio          below_one 
#>          2057.4737            18.1614             0.0000             0.0000 
#>                min                max 
#>             2.5325            21.5074 
#> 
#> Design effect (Kish): 1.602
#> Effective sample size: 234.126
```
