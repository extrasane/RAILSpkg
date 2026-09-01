# Build a RAILS reference sample from recoded PUMS

Takes recoded PUMS microdata and produces the three things RAILS needs
from a reference sample: cleaned person-level records, the aggregated
cell table, and the population margins.

## Usage

``` r
pums_reference(
  data,
  vars = c("agegroup", "sex", "edu", "homeown", "income", "race_eth", "region"),
  order = 3,
  rescale = TRUE
)
```

## Arguments

- data:

  Recoded PUMS from
  [`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md),
  or any data frame with a `weight` column and the columns in `vars`.

- vars:

  Covariates to aggregate over. Defaults to the seven used in the
  application.

- order:

  Highest interaction order to compute margins for. Must be at least the
  `scope` you will pass to
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).

- rescale:

  Rescale the weights of the retained records to the population total
  before incomplete records were dropped. Leave `TRUE` unless you have a
  reason to want the reference sample to under-count.

## Value

An object of class `pums_reference`: a list with `micro` (the retained
person-level records), `cells` (the aggregated table, ready for
`rails(..., aggregated = TRUE)`), `totals` (the named margin vector for
`pop_totals`), `n_pop` (the population total the weights sum to),
`n_dropped` and `vars`.

## Details

Records with a missing value in any of `vars` cannot enter a cell, so
they are dropped – but dropping them would also lose their share of the
population. The weights are therefore rescaled so the retained records
still sum to the population total computed *before* the drop. This is
the step most easily got wrong by hand, and it is why the totals are
worth taking from here rather than from an ad-hoc aggregation.

## See also

[`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md),
[`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md),
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).

## Examples

``` r
# Standing in for real PUMS, with the same column names and levels:
pop <- rails_simulate(5000, seed = 1)
fake <- data.frame(weight = 1 / pop$ps, pop[, c("agegroup", "sex", "income")])

ref <- pums_reference(fake, vars = c("agegroup", "sex", "income"), order = 2)
ref
#> PUMS reference sample
#>   covariates : agegroup, sex, income
#>   records    : 5,000
#>   cells      : 24
#>   margins    : 18 up to order 2
#>   population : 21,198
head(ref$cells)
#>   agegroup    sex  income   weight   n weight_sq
#> 1    young female  p10-50 1855.646 399  8696.863
#> 2   middle   male p90-100  304.461  70  1329.445
#> 3    young   male  p50-90 1187.199 213  6667.814
#> 4    older female  p10-50 1033.802 308  3484.918
#> 5   middle female  p10-50 2265.910 578  8918.822
#> 6   middle   male  p50-90 1497.776 316  7126.524
```
