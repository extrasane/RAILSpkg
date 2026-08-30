# Extract RAILS weights

Extract RAILS weights

## Usage

``` r
# S3 method for class 'rails_fit'
weights(object, ...)

# S3 method for class 'rails_subgroup_fit'
weights(object, ...)
```

## Arguments

- object:

  A `rails_fit` or `rails_subgroup_fit`.

- ...:

  Ignored.

## Value

A numeric vector: one weight per record of the non-probability sample
that was fitted, or one per cell when the fit was built from cell
tables. Records dropped for missing covariates, and strata with no
reference records, are `NA`.

## Examples

``` r
pop <- rails_simulate(2000, seed = 3)
ref <- pop[pop$s == 1, ]; ref$dweight <- 1 / ref$ps
fit <- rails(pop[pop$aou == 1, ], ref, c("agegroup", "sex"),
             weights_ref = "dweight", verbose = FALSE)
summary(weights(fit))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   2.532   2.532   4.110   5.487   4.110  21.507 
```
