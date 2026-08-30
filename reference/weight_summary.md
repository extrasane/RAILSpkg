# Weight diagnostics

Summarizes a vector of weights for quality checking after estimation.

## Usage

``` r
weight_summary(w, ...)
```

## Arguments

- w:

  Numeric vector of weights.

- ...:

  Ignored; present for call compatibility.

## Value

A named numeric vector with elements `sum`, `var`, `non_positive_ratio`,
`below_one`, `min`, and `max`.

## Examples

``` r
weight_summary(c(1.2, 0.8, 3.4, 0.5))
#>                sum                var non_positive_ratio          below_one 
#>           5.900000           1.729167           0.000000           2.000000 
#>                min                max 
#>           0.500000           3.400000 
```
