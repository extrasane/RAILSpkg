# Deprecated functions

The pre-package function names, kept so that existing analysis scripts
keep running. Each forwards to its replacement and warns once per
session.

## Usage

``` r
fun.rails.threeway(
  dt_agg_aou,
  dt_agg_pums,
  pop_totals,
  names_univar = c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
  alpha = 0.05,
  nsiz = sum(dt_agg_pums$weight)
)

fun.sub.rails.threeway(
  dt_agg_aou,
  dt_agg_pums,
  subgroup_var = "region",
  names_univar = setdiff(c("agegroup", "edu", "homeown", "income", "race_eth", "sex",
    "region"), subgroup_var),
  alpha = 0.05
)

fun.nps(cat_temp, dt_aou, dt_s, nsiz = sum(dt_s$weight), theta_init = NULL)

fun.lkd(x, names_univar, LKHD, dt_s, dt_aou, m_s)

fun.out(w, ...)
```

## Arguments

- dt_agg_aou:

  Aggregated non-probability cell counts (column `weight`).

- dt_agg_pums:

  Aggregated reference-sample cell counts (column `weight`).

- pop_totals:

  Named numeric vector of population margins.

- names_univar:

  Character vector of main-effect variables.

- alpha:

  Significance threshold for forward selection.

- nsiz:

  Population total the weights are scaled to.

- subgroup_var:

  Name of the stratifying variable.

- cat_temp:

  A one-sided model formula.

- dt_aou:

  Aggregated non-probability cell counts.

- dt_s:

  Aggregated reference-sample cell counts.

- theta_init:

  Optional warm-start coefficient vector.

- x:

  Candidate term.

- LKHD:

  Log-likelihood of the current model.

- m_s:

  Model matrix of the current model on the reference cells.

- w:

  Numeric vector of weights.

- ...:

  Ignored.

## Value

For `fun.rails.threeway()` and `fun.sub.rails.threeway()`, the wide data
frame described above.

## Details

- `fun.rails.threeway()`:

  use
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
  whose defaults (`start = 2`, `scope = 3`, `fallback = FALSE`)
  reproduce it exactly

- `fun.sub.rails.threeway()`:

  use
  [`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md)

- `fun.nps()`:

  use
  [`nps_weights()`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md)

- `fun.lkd()`:

  use
  [`gvs_step()`](https://extrasane.github.io/RAILSpkg/reference/gvs_step.md)

- `fun.out()`:

  use
  [`weight_summary()`](https://extrasane.github.io/RAILSpkg/reference/weight_summary.md)

The two estimator shims return the old wide data frame – `dt_agg_aou`
with `d_rails`, `d_cal1`, `d_cal2`, `d_nps1`, `d_nps2`, `d_nps1_rake`,
`d_nps2_rake`, `d_unweighted`, `selected_terms` and `calibrated_terms`
appended – not a `rails_fit`. The weights themselves are identical to
what
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
produces on the same input under its defaults; only the container
differs.
