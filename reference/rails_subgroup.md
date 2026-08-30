# Subgroup RAILS

Runs
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
independently within each level of a stratifying variable. Each subgroup
gets its own NPS model, its own selection path, its own stepwise walk,
and its own population margins computed from that subgroup's reference
records, with weights scaled to the subgroup total.

## Usage

``` r
rails_subgroup(
  np,
  ref,
  vars,
  by,
  weights_ref = NULL,
  weights_np = NULL,
  aggregated = FALSE,
  start = 2,
  scope = 3,
  alpha = 0.05,
  benchmarks = FALSE,
  fallback = FALSE,
  lifo = c("descending", "ascending"),
  keep_data = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- np, ref, vars, weights_ref, weights_np:

  As in
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).

- by:

  Name of the stratifying variable. It must be a column of both samples,
  and is automatically excluded from `vars`.

- aggregated, start, scope, alpha, benchmarks, fallback, lifo,
  keep_data, verbose:

  As in
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
  applied within every subgroup.

- ...:

  Further arguments passed to
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).

## Value

An object of class `rails_subgroup_fit`: a list with `fits` (one
`rails_fit` per level, named by level), `weights` (a single vector
aligned to the rows of `np`, pooling every subgroup), `by`, `levels`,
and `call`.

## Details

This is the *subgroup calibration* strategy – "calibration-in-sub" –
that Chen et al. compare against *global calibration*, in which one set
of weights is built from the whole cohort by
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md) and
then simply applied within a subgroup. Which strategy does better is an
empirical question rather than a settled one; that comparison is what
the subgroup paper is about, and heterogeneity of the non-probability
selection mechanism across subgroups is the factor it identifies as most
influential. This function provides the subgroup arm;
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
provides the global arm.

## Empty strata

Levels with no reference records are skipped with a warning and their
non-probability records receive `NA` weights. Factor levels are never
dropped, so model matrices stay aligned between the two samples within
each stratum.

## See also

[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).

## Examples

``` r
pop <- rails_simulate(20000, seed = 7)
aou <- pop[pop$aou == 1, ]
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps

fit <- rails_subgroup(aou, ref, vars = c("agegroup", "income"),
                      by = "sex", weights_ref = "dweight",
                      verbose = FALSE)
fit
#> Subgroup RAILS fit
#>   stratifier : sex (2 of 2 level(s) fitted)
#>   covariates : agegroup, income
#> 
#>   sex = male: 0 selected, 0 raked, 12 cells
#>   sex = female: 0 selected, 0 raked, 12 cells
#> 
#>   weights: n = 3826, sum = 19,944.8
```
