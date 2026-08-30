# ACS PUMS variables RAILS uses

The Census API variable names fetched by
[`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md),
with what each is used for. Supplied as a function rather than a
constant so it can be printed, subsetted, or extended before being
passed back to
[`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md).

## Usage

``` r
pums_variables(geography = c("region", "state"))
```

## Arguments

- geography:

  `"region"` (the default) or `"state"`. The state form adds `ST`, the
  state FIPS code.

## Value

A named character vector: names are the ACS variable codes, values
describe what each contributes.

## See also

[`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md),
[`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md).

## Examples

``` r
pums_variables()
#>                                           PWGTP 
#>                                 "person weight" 
#>                                            AGEP 
#> "age, used for the adult filter and age groups" 
#>                                             SEX 
#>                                           "sex" 
#>                                           HINCP 
#>                      "household income, banded" 
#>                                            SCHL 
#>                "educational attainment, banded" 
#>                                             TEN 
#>            "tenure, banded into own/rent/other" 
#>                                            HISP 
#>                               "Hispanic origin" 
#>                                          RACWHT 
#>                                   "race: White" 
#>                                          RACBLK 
#>                                   "race: Black" 
#>                                          RACASN 
#>                                   "race: Asian" 
#>                                         RACAIAN 
#>        "race: American Indian or Alaska Native" 
#>                                           RACNH 
#>                         "race: Native Hawaiian" 
#>                                           RACPI 
#>                        "race: Pacific Islander" 
#>                                          RACSOR 
#>                         "race: some other race" 
#>                                          REGION 
#>                                 "Census region" 
```
