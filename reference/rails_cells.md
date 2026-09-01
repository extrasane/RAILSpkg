# Aggregate microdata into covariate cells

Collapses individual-level records into one row per unique combination
of `vars`, with a `weight` column giving the cell total. RAILS runs
entirely on these cells, so both the non-probability sample and the
reference sample must be aggregated the same way before estimation.
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
calls this for you when given microdata; call it directly when you want
to inspect or cache the cell tables.

## Usage

``` r
rails_cells(data, vars, weights = NULL, drop_na = TRUE)
```

## Arguments

- data:

  A data frame of individual records.

- vars:

  Character vector of covariate names to aggregate over.

- weights:

  Optional. Either a numeric vector of length `nrow(data)` or the name
  of a column of `data` holding design weights. `NULL` (the default)
  counts each row as 1, which is what you want for an unweighted biobank
  sample; pass the survey weight for a probability reference sample.

- drop_na:

  Drop rows with a missing value in any of `vars`. Model matrices would
  silently drop them anyway, so the default warns and removes them up
  front.

## Value

A data frame with one row per cell: the columns named in `vars`, plus
`weight`, `n` (records in the cell) and `weight_sq` (sum of squared
weights). The last two are what let
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
compute the stacked sandwich from cells alone. It carries class
`rails_cells` and an attribute `row_cell`, the integer index mapping
each row of `data` to its cell, which is what lets
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
return one weight per input record.

## Details

Factor levels are **not** dropped, so cells present in one sample but
empty in the other still line up in the model matrices.

## See also

[`rails_totals()`](https://extrasane.github.io/RAILSpkg/reference/rails_totals.md)
for the matching population margins.

## Examples

``` r
pop <- rails_simulate(500, seed = 1)
rails_cells(pop[pop$aou == 1, ], c("agegroup", "sex"))
#>   agegroup    sex weight  n weight_sq
#> 1    young female     15 15        15
#> 2   middle female     32 32        32
#> 3    older female     27 27        27
#> 4    older   male      7  7         7
#> 5   middle   male     10 10        10
#> 6    young   male      2  2         2

# A reference sample carrying design weights:
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps
rails_cells(ref, c("agegroup", "sex"), weights = "dweight")
#>   agegroup    sex    weight  n weight_sq
#> 1    older   male  30.11030  7  130.7037
#> 2    young   male  56.48730 10  321.5076
#> 3   middle female 167.50293 46  615.8940
#> 4    young female 140.59382 34  592.7221
#> 5   middle   male 119.02028 25  573.0375
#> 6    older female  70.20786 22  227.0849
```
