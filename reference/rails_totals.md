# Population margins matching a cell table

Builds the named vector of population totals that
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
rakes to: every one-way, two-way, ... up to `order`-way margin of
`vars`, named exactly as the calibration model matrix columns are named.
Assembling this by hand is the most common source of errors, because the
names have to match what
[`survey::calibrate()`](https://rdrr.io/pkg/survey/man/calibrate.html)
expects.

## Usage

``` r
rails_totals(
  cells,
  vars = attr(cells, "vars"),
  order = 3,
  weight_col = "weight"
)
```

## Arguments

- cells:

  A cell table from
  [`rails_cells()`](https://extrasane.github.io/RAILSpkg/reference/rails_cells.md),
  or any data frame with the columns in `vars` plus a `weight` column.

- vars:

  Character vector of covariate names. Defaults to the variables
  recorded on a `rails_cells` object.

- order:

  Highest interaction order to compute margins for. Must be at least as
  large as the highest order
  [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
  will select, otherwise raking has no target for the selected terms.

- weight_col:

  Name of the cell-total column.

## Value

A named numeric vector of margins, suitable as the `pop_totals` argument
of [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md).

## See also

[`rails_cells()`](https://extrasane.github.io/RAILSpkg/reference/rails_cells.md).

## Examples

``` r
pop <- rails_simulate(500, seed = 1)
ref <- rails_cells(pop[pop$s == 1, ], c("agegroup", "sex"))
rails_totals(ref, order = 2)
#>              (Intercept)           agegroupmiddle            agegroupolder 
#>                      144                       71                       29 
#>                sexfemale agegroupmiddle:sexfemale  agegroupolder:sexfemale 
#>                      102                       46                       22 
```
