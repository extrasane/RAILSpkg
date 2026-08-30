# Recode raw ACS PUMS into RAILS analysis variables

Applies the banding and labelling used in the RAILS application: age
groups, sex, race and ethnicity, income bands, educational attainment,
tenure, and Census region. Every result is a factor with fixed levels,
so cell tables built from different extracts line up.

## Usage

``` r
pums_recode(data, year = 2022, adults_only = TRUE)
```

## Arguments

- data:

  Raw PUMS records from
  [`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md),
  or a data frame with the same ACS column names.

- year:

  The ACS vintage `data` came from. Only used to check the codings below
  are the right ones; see the note.

- adults_only:

  Keep only records with `AGEP > 17`. `TRUE` for a reference sample of
  adults. Set `FALSE` when you need complete population counts, for
  instance to compute state population shares.

## Value

A data frame with one row per retained person and columns `weight`,
`age`, `agegroup`, `sex`, `race_eth`, `income`, `edu`, `homeown`,
`region`, and, when `ST` was fetched, `state` as a two-letter USPS code.
Demographic columns may be `NA` where a code fell outside the bands.

## Codings are vintage-specific

The `SCHL`, `TEN`, `SEX`, `HISP`, `RAC*` and `HINCP` bandings follow the
**2022** ACS 1-year PUMS data dictionary. Category codes do change
between releases, so a different `year` raises a warning: check each
mapping against that year's dictionary at
<https://www.census.gov/programs-surveys/acs/microdata/documentation.html>
before trusting the output.

Income bands are `<35k`, `35k-50k`, `50k-75k`, `75k-100k`, `>100k`; note
that `HINCP` can be negative, and values below -60000 are treated as
missing rather than as the lowest band.

## See also

[`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md),
[`pums_reference()`](https://extrasane.github.io/RAILSpkg/reference/pums_reference.md).

## Examples

``` r
# A handful of synthetic records in raw ACS coding:
raw <- data.frame(
  PWGTP = c(100, 250, 90), AGEP = c(30, 70, 16), SEX = c(1, 2, 2),
  HINCP = c(42000, 120000, 20000), SCHL = c(21, 16, 12), TEN = c(1, 3, 2),
  HISP = c(1, 5, 1), RACWHT = c(1, 0, 0), RACBLK = c(0, 0, 1),
  RACASN = c(0, 0, 0), REGION = c(1, 4, 3)
)
pums_recode(raw)
#>   weight age agegroup    sex race_eth  income                          edu
#> 1    100  30    25-44   Male NH White 35k-50k College graduate or advanced
#> 2    250  70    65-74 Female Hispanic   >100k          Highschool graduate
#>   homeown    region
#> 1     Own Northeast
#> 2    Rent      West
```
