# Variance of a RAILS-weighted mean

Standard error and confidence interval for a weighted mean computed with
RAILS weights, in either of two forms.

## Usage

``` r
rails_var(
  fit,
  y,
  type = c("simplified", "stacked", "both"),
  truth = NULL,
  level = 0.95
)
```

## Arguments

- fit:

  A `rails_fit` from
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
  or a `rails_subgroup_fit` from
  [`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md).
  Fits built with `aggregated = TRUE` are supported; see `y` and the
  section below.

- y:

  The outcome. For a fit built from microdata, a vector aligned to the
  rows of the non-probability sample passed to
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).
  For a fit built from cell tables, a data frame with columns `sum` and
  `sumsq` – the within-cell sum of the outcome and of its square – one
  row per cell, in the order of the cell table. For a 0/1 outcome the
  two columns are equal, so `sumsq` is just the case count again. Binary
  or continuous.

- type:

  Which variance to compute: `"simplified"` (the default), `"stacked"`,
  or `"both"`.

- truth:

  Optional known population value. When supplied, the returned vector
  includes a coverage indicator for each interval computed, which is
  what the simulation studies use. For a subgroup fit, either one value
  for every stratum or a vector named by stratum level.

- level:

  Confidence level.

## Value

For a single fit, a named numeric vector: `mu_hat`, the variance,
standard error and interval bounds for the requested form, a coverage
indicator, and the sample sizes. `"stacked"` and `"both"` add the six
sandwich components (`term_11` through `term_33`); `"both"` also adds
`ratio_full_naive`. Simplified quantities are named `*_naive`, stacked
ones `*_full`.

For a subgroup fit, a data frame with one row per stratum: the
stratifier in the first column, then those same quantities.

## Details

The **simplified** variance (the default) treats the weights as fixed.
It is what the simulation pipeline reports and what `svymean()` on a
plain design gives, and it is cheap: no model matrix is involved.

The **stacked** variance accounts for the weights having been estimated,
through the sandwich of three stacked estimating equations – (1) the
pseudo-likelihood propensity score, (2) the calibration/raking
constraints, and (3) the weighted mean itself. Because the simplified
form ignores the first two, the two forms differ; `type = "both"`
reports the pair along with `ratio_full_naive`, the ratio of stacked to
simplified.

The direction is not fixed. Ignoring the propensity stage tends to
understate, but raking to known population margins removes variance the
way a regression estimator does, so the stacked form is often the
*smaller* of the two – in the examples below the ratio is around 0.97.
Report the stacked figure because it is the one that accounts for how
the weights were built, not because it is conservative.

## Subgroup fits

Each stratum is a separate fit – its own propensity model, its own
selected terms, its own margins – so the variance is computed within
each and returned one row per stratum.

The subgroup paper reports the **simplified** variance by default,
because residual bias persists there and a tighter interval would
overstate what the estimator delivers. The stacked form is available and
applies within a stratum, but the per-stratum sandwich carries no
covariance between strata, so it warns: use those intervals stratum by
stratum, and not for anything that pools or differences across strata.

## Everything runs on cells

The variance is computed from the same aggregated cells the estimator
uses, not from individual records. Within a covariate cell the model
matrix row, the weight and the fitted propensity are all constant, so
each block of the sandwich needs at most three numbers per cell: the
record count, the sum of the outcome, and the sum of its square. Blocks
free of the outcome need only the count; blocks linear in it need the
sum; the one quadratic block needs the sum of squares. On the reference
side the meat is quadratic in the design weights, so
[`rails_cells()`](https://extrasane.github.io/RAILSpkg/reference/rails_cells.md)
records `weight_sq` alongside `weight`.

Cost therefore scales with the number of cells rather than the number of
records, which for a biobank is a difference of two or three orders of
magnitude. Given microdata, `rails_var()` reduces it to these statistics
itself; given cell tables, supply them through `y`. `keep_data = TRUE`
is no longer required for the stacked form.

## See also

[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
[`rails_design()`](https://extrasane.github.io/RAILSpkg/reference/rails_design.md).

## Examples

``` r
pop <- rails_simulate(20000, seed = 42)
aou <- pop[pop$aou == 1, ]
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps

fit <- rails(aou, ref, vars = c("agegroup", "sex", "income"),
             weights_ref = "dweight", verbose = FALSE)

# Simplified, the default:
rails_var(fit, aou$y, truth = mean(pop$y))
#>         mu_hat      var_naive       se_naive ci_naive_lower ci_naive_upper 
#>   2.040237e-01   5.769617e-05   7.595799e-03   1.891362e-01   2.189112e-01 
#>    cover_naive          N_hat           n_NP 
#>   1.000000e+00   1.987578e+04   3.799000e+03 

# Both, to compare the two:
v <- rails_var(fit, aou$y, type = "both", truth = mean(pop$y))
v[c("mu_hat", "se_naive", "se_full", "ratio_full_naive")]
#>           mu_hat         se_naive          se_full ratio_full_naive 
#>      0.204023738      0.007595799      0.007404475      0.974811841 

# One row per stratum for a subgroup fit:
sfit <- rails_subgroup(aou, ref, c("agegroup", "income"), by = "sex",
                       weights_ref = "dweight", start = 1, scope = 2,
                       verbose = FALSE)
rails_var(sfit, aou$y)
#>      sex    mu_hat    var_naive    se_naive ci_naive_lower ci_naive_upper
#> 1   male 0.1215393 2.136081e-04 0.014615339     0.09289375      0.1501848
#> 2 female 0.2494814 6.939087e-05 0.008330118     0.23315462      0.2658081
#>   cover_naive     N_hat n_NP
#> 1          NA  7061.835  571
#> 2          NA 12813.949 3228
```
