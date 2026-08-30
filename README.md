# RAILS

<!-- badges: start -->
<!-- badges: end -->

**R**aking with **A**ssisted Nested Propens**i**ty Score and **L**IFO **S**election
(Chen et al., 2025). An R package for combining a non-probability sample — a volunteer biobank, a convenience cohort, an opt-in
panel — with a probability reference sample, producing weights that make the
non-probability sample represent the target population.

The estimator has three stages, which are what it is named for: a nested
propensity score (NPS) fitted by pseudo-likelihood against the reference sample;
greedy variable selection (GVS) of higher-order interactions by forward
likelihood-ratio test; and LIFO stepwise raking to the reference population
margins -- rake the full selected model, and pop the last-added term until
raking converges. Everything runs on aggregated covariate
cells, so cost scales with the number of distinct covariate patterns rather than
the sample size.

## Installation

```r
# install.packages("remotes")
remotes::install_github("extrasane/RAILSpkg")
```

## Usage

```r
library(RAILS)

# A synthetic population: `aou` is the biased volunteer sample,
# `s` is a probability sample with design weights 1 / ps.
pop <- rails_simulate(20000, seed = 1)
aou <- pop[pop$aou == 1, ]
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps

fit <- rails(aou, ref,
             vars = c("agegroup", "sex", "income"),
             weights_ref = "dweight")

summary(fit)
w <- weights(fit)
```

The weights correct the selection bias:

```r
mean(pop$y)                      # population truth
mean(aou$y)                      # unweighted volunteer sample
sum(w * aou$y) / sum(w)          # RAILS
```

Standard errors. `rails_var()` gives the simplified (fixed-weight) variance by
default, and the stacked variance that accounts for the weights having been
*estimated* on request:

```r
rails_var(fit, aou$y)                      # simplified, the default
rails_var(fit, aou$y, type = "both")       # and the stacked correction
```

And for everything else, hand the weights to **survey**:

```r
des <- rails_design(fit, aou)
survey::svyglm(y ~ agegroup, design = des, family = quasibinomial())
```

## Functions

| Function | Purpose |
|---|---|
| `rails()` | The estimator: propensity model, forward selection, stepwise raking |
| `rails_subgroup()` | Subgroup calibration: runs the estimator independently within each stratum |
| `rails_var()` | Variance of a weighted mean, global or per stratum: simplified by default, stacked on request |
| `rails_cells()` | Aggregate microdata into covariate cells |
| `rails_totals()` | Population margins, named to match the calibration model matrix |
| `rails_design()` | Wrap the weights in a `survey` design |
| `rails_simulate()` | Synthetic population for examples and testing |
| `pums_fetch()` | Download ACS PUMS microdata from the Census API |
| `pums_recode()` | Band raw ACS codes into the analysis variables |
| `pums_reference()` | Cleaned records, cell table and margins in one object |
| `nps_weights()` | Nested propensity score weights (one fit) |
| `gvs_step()` | Greedy variable selection: likelihood-ratio test for one candidate |
| `weight_summary()` | Weight diagnostics |

## Building the reference sample from ACS PUMS

The application uses ACS PUMS as its probability reference. Three functions
cover that end-to-end:

```r
Sys.setenv(CENSUS_API_KEY = "...")     # better: put it in ~/.Renviron

raw <- pums_fetch(2022, cache = "PUMS_2022.csv")   # cached after the first run
ref <- pums_reference(pums_recode(raw))

ref
#> PUMS reference sample
#>   covariates : agegroup, sex, edu, homeown, income, race_eth, region
#>   records    : ...
#>   cells      : ...
#>   margins    : ... 
#>   population : ...

fit <- rails(rails_cells(cohort, ref$vars), ref$cells, ref$vars,
             pop_totals = ref$totals, aggregated = TRUE)
```

`pums_recode()` applies the 2022 ACS bandings for age, sex, race and ethnicity,
income, education, tenure and region, as fixed-level factors so extracts taken
at different times line up. It warns if you hand it another vintage, because
ACS category codes change between releases.

`pums_reference()` handles the step that is easiest to get wrong by hand:
records missing any covariate cannot enter a cell, but dropping them would also
drop their share of the population, so the retained weights are rescaled to the
population total computed *before* the drop. It returns the cleaned records, the
aggregated cell table, and the margins together.

## Two notes on inference

`rails_var()` needs individual records, not cell tables -- it needs one `y` per
record. The stacked form needs them for a second reason: within a covariate cell
the model matrix row, the weight and the fitted propensity are all constant, but
the outcome is not. Fits built with `aggregated = TRUE` cannot be passed to it.

Standard errors taken off a `rails_design()` object treat the weights as fixed,
matching `rails_var()`'s simplified form. For the correction, use
`rails_var(type = "both")`; the `ratio_full_naive` it reports is the ratio of
the two, which can fall either side of 1 -- raking to known margins removes
variance, so the corrected figure is often the smaller one.

## Migrating from the script version

The pre-package names still work and forward to their replacements with a
deprecation warning:

| Old | New |
|---|---|
| `fun.rails.threeway()` | `rails()` (its defaults reproduce it exactly) |
| `fun.sub.rails.threeway()` | `rails_subgroup()` |
| `fun.nps()` | `nps_weights()` |
| `fun.lkd()` | `gvs_step()` |
| `fun.out()` | `weight_summary()` |

`fun.rails.threeway()` still returns exactly what it always did. `rails()`
matches it too, with one argument: the application code ran the LIFO stage
ascending (climb from the starting model, halt at the first failure) where the
manuscripts define it as descending, so `rails()` defaults to `lifo =
"descending"` and `lifo = "ascending"` reproduces the published application.
The two agree whenever raking convergence is monotone, which covers every case
in the package's own tests. What changed
is that the model is now specified rather than fixed: `start` is the model taken
as given, `scope` the largest model selection may reach. `start = 2, scope = 3`
is the published estimator; `scope = 2` gives a two-way-only fit with no
selection; `start = 1, scope = 3` offers the two-way terms for selection too.
Each accepts an interaction order, a formula, or an explicit vector of terms.

## Citation

```r
citation("RAILS")
```

## License

MIT © Huiding Chen

## The paper code lives elsewhere

This repository is the R package. The simulations, the *All of Us* / PUMS
application pipeline, and the results behind the two papers are at
[extrasane/RAILS](https://github.com/extrasane/RAILS).
