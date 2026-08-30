# Plot the RAILS weight distribution

Draws a histogram of the fitted weights on a log scale, which is where a
heavy right tail – the usual failure mode of propensity weighting – is
visible.

## Usage

``` r
# S3 method for class 'rails_fit'
plot(x, breaks = 40, ...)
```

## Arguments

- x:

  A `rails_fit`.

- breaks:

  Passed to [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html).

- ...:

  Passed to [`graphics::hist()`](https://rdrr.io/r/graphics/hist.html).

## Value

`x`, invisibly. Called for the plot.

## Examples

``` r
pop <- rails_simulate(2000, seed = 3)
ref <- pop[pop$s == 1, ]; ref$dweight <- 1 / ref$ps
fit <- rails(pop[pop$aou == 1, ], ref, c("agegroup", "sex"),
             weights_ref = "dweight", verbose = FALSE)
plot(fit)

```
