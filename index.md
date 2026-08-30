# RAILS

**R**aking with **A**ssisted Nested Propens**i**ty Score and **L**IFO
**S**election (Chen et al., 2025). An R package for combining a
non-probability sample — a volunteer biobank, a convenience cohort, an
opt-in panel — with a probability reference sample, producing weights
that make the non-probability sample represent the target population.

The estimator has three stages, which are what it is named for: a nested
propensity score (NPS) fitted by pseudo-likelihood against the reference
sample; greedy variable selection (GVS) of higher-order interactions by
forward likelihood-ratio test; and LIFO stepwise raking to the reference
population margins – rake the full selected model, and pop the
last-added term until raking converges. Everything runs on aggregated
covariate cells, so cost scales with the number of distinct covariate
patterns rather than the sample size.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("extrasane/RAILSpkg")
```

## Usage

### A worked example, end to end

[`rails_simulate()`](https://extrasane.github.io/RAILSpkg/reference/rails_simulate.md)
gives a population where the truth is known: `aou` is a biased volunteer
sample, `s` is a probability sample carrying design weights.

``` r

library(RAILS)

pop <- rails_simulate(20000, seed = 42)
aou <- pop[pop$aou == 1, ]        # the volunteer sample
ref <- pop[pop$s == 1, ]          # the probability sample
ref$dweight <- 1 / ref$ps         # its design weights
```

[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
needs the two samples, the covariates they share, and the reference
sample’s design weights. `start` and `scope` say which models the search
runs between — see below.

``` r

fit <- rails(aou, ref,
             vars        = c("agegroup", "sex", "income"),
             weights_ref = "dweight",
             start = 1, scope = 2)
fit
#> RAILS fit
#>   covariates   : agegroup, sex, income
#>   start model  : 3 term(s)
#>   scope        : 6 term(s), so 3 candidate(s)
#>   selected     : 2 term(s) at alpha = 0.05
#>                  agegroup:sex, sex:income
#>   raked        : 2 of 2 selected term(s) survived the stepwise walk
#>   cells        : 24
#>   weights      : n = 3799, sum = 19,875.78, range 1.65 to 25.9
```

The weights remove most of the selection bias:

``` r

w <- weights(fit)
c(truth      = mean(pop$y),
  unweighted = mean(aou$y),
  rails      = sum(w * aou$y) / sum(w))
#>      truth unweighted      rails
#>     0.1952     0.2545     0.2051
```

### Preparing the inputs

[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
takes individual records and aggregates them for you. It only ever works
on covariate cells, so cost scales with the number of distinct covariate
patterns, not the sample size.

When the microdata will not fit in one data frame — or when you want to
inspect what is being fitted — build the pieces yourself:

``` r

vars      <- c("agegroup", "sex", "income")
cells_np  <- rails_cells(aou, vars)                          # one row per cell
cells_ref <- rails_cells(ref, vars, weights = "dweight")
totals    <- rails_totals(cells_ref, order = 2)              # population margins

fit <- rails(cells_np, cells_ref, vars,
             pop_totals = totals, aggregated = TRUE)
```

[`rails_totals()`](https://extrasane.github.io/RAILSpkg/reference/rails_totals.md)
names the margins exactly as the calibration model matrix names them,
which is the step most easily got wrong by hand. Supply your own
`pop_totals` when the margins come from an external table rather than
from the reference sample.

Two consequences of `aggregated = TRUE`:
[`weights()`](https://rdrr.io/r/stats/weights.html) returns one weight
per cell rather than per record, and
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
will not run, since the variance needs one outcome per record.

### Choosing which models the search runs between

`start` is the model taken as given; `scope` is the largest model
selection may reach. Everything in `scope` but not in `start` is a
candidate for the forward test. Each accepts an interaction order, a
one-sided formula, or an explicit character vector of terms.

``` r

rails(aou, ref, vars, weights_ref = "dweight")                    # start = 2, scope = 3
rails(aou, ref, vars, weights_ref = "dweight", scope = 2)         # no selection: just rake
rails(aou, ref, vars, weights_ref = "dweight", start = 1, scope = 2)
rails(aou, ref, vars, weights_ref = "dweight",
      start = ~ (agegroup + sex + income)^2,
      scope = ~ (agegroup + sex + income)^3)                      # same as the default
```

The defaults `start = 2, scope = 3` are the published estimator. Term
components may be written in any order: `"sex:agegroup"` matches
`"agegroup:sex"`.

### Reading the fit

``` r

summary(fit)
#> Selected terms (2):
#>   1. agegroup:sex
#>   2. sex:income
#>
#> Raked to 5 term(s); full walk completed.
#>
#> Weight diagnostics:
#>                sum                var non_positive_ratio          below_one
#>         19875.7840            16.4808             0.0000             0.0000
#>                min                max
#>             1.6543            25.8569
#>
#> Design effect (Kish): 1.602
#> Effective sample size: 2371.5
```

The design effect is the price of weighting — the factor by which
unequal weights inflate the variance of a mean — so an effective sample
size well below the nominal one is expected. A very large one means a
handful of records carry the estimate. `plot(fit)` shows the weight
distribution on a log scale, which is where a heavy right tail is
visible.

Useful fields on the fit: `terms_selected`, `terms_used` (what was
actually raked to), `converged`, and `cells`.

### Standard errors

[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
reports the **simplified** variance by default, which treats the weights
as fixed. That is what the papers report.

``` r

rails_var(fit, aou$y)
#>         mu_hat      var_naive       se_naive ci_naive_lower ci_naive_upper
#>        0.20512        0.00006        0.00762        0.19019        0.22006
```

The **stacked** variance additionally accounts for the weights having
been estimated, through a sandwich over the propensity, calibration and
mean stages:

``` r

rails_var(fit, aou$y, type = "stacked")
rails_var(fit, aou$y, type = "both")     # both, plus ratio_full_naive
```

`ratio_full_naive` can fall either side of 1. Ignoring the propensity
stage tends to understate, but raking to *known* margins removes
variance the way a regression estimator does, and the second effect
often wins.

Pass `truth =` to get coverage indicators, which is what the simulation
studies use.

### Downstream analysis

``` r

des <- rails_design(fit, aou)
survey::svymean(~y, des)
survey::svyglm(y ~ agegroup, design = des, family = quasibinomial())
```

Point estimates from this design are exact. Its standard errors are
survey’s own design-based ones, close to but not identical with
`rails_var(type = "simplified")`, and neither accounts for the weights
having been estimated. Use
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
for inference on a mean.

### Subgroup calibration

[`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md)
runs the whole procedure independently within each level of a
stratifier: its own propensity model, selection path, stepwise walk and
margins, with weights scaled to that stratum’s total.

``` r

sfit <- rails_subgroup(aou, ref,
                       vars = c("sex", "income"),
                       by   = "agegroup",
                       weights_ref = "dweight",
                       start = 1, scope = 2)

weights(sfit)                       # pooled, aligned to the rows of `aou`
sfit$fits[["young"]]                # the individual stratum fits
rails_var(sfit, aou$y)              # one row per stratum
```

`by` must not appear in `vars` — within a stratum it is constant. Which
variable you stratify by decides what is left for selection to find: an
interaction involving the stratifier flattens into a main effect inside
each stratum. See
[`vignette("subgroup-rails")`](https://extrasane.github.io/RAILSpkg/articles/subgroup-rails.md).

### When raking does not converge

Two arguments control what happens then.

`lifo` sets the direction of the stepwise stage. `"descending"` (the
default) is the method as published: rake the full selected model, and
pop the last-added term until raking converges. `"ascending"` climbs
from the starting model and halts at the first failure, which is what
the application code did — use it to reproduce published application
results exactly.

`fallback` decides what comes back when selection picks nothing, or when
no selected model rakes. `FALSE` (the default) returns `NA` weights,
which is what the published code does and makes a non-result impossible
to miss. `TRUE` returns the raked starting model instead. Either way
`fit$converged` records which happened.

### Reproducing the published estimator

The defaults `start = 2`, `scope = 3`, `fallback = FALSE`, plus
`lifo = "ascending"`, are exactly the procedure of the application code.
The deprecated
[`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
is pinned to that and returns the old wide data frame.
`tests/testthat/test-reproduces-original.R` asserts the package still
matches the original implementation’s output.

## Function reference

Every function has a help page with arguments, return value and a
runnable example:
[`?rails`](https://extrasane.github.io/RAILSpkg/reference/rails.md),
[`?rails_var`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md),
and so on. [`help(package = "RAILS")`](https://rdrr.io/pkg/RAILS/man)
lists them all.

**Fitting**

| Function | Purpose |
|----|----|
| [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md) | The estimator: NPS, greedy variable selection, LIFO stepwise raking |
| [`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md) | Subgroup calibration: the estimator run independently within each stratum |

**Inspecting a fit** — S3 methods on `rails_fit`

| Function | Purpose |
|----|----|
| [`weights()`](https://rdrr.io/r/stats/weights.html) | The fitted weights, one per record |
| [`summary()`](https://rdrr.io/r/base/summary.html) | Selection path, weight diagnostics, design effect, effective sample size |
| [`print()`](https://rdrr.io/r/base/print.html) | One-screen overview of the fit |
| [`plot()`](https://rdrr.io/r/graphics/plot.default.html) | Weight distribution on a log scale |
| [`weight_summary()`](https://extrasane.github.io/RAILSpkg/reference/weight_summary.md) | Weight diagnostics for any weight vector |

**Inference**

| Function | Purpose |
|----|----|
| [`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md) | Variance of a weighted mean, global or per stratum: simplified by default, stacked on request |
| [`rails_design()`](https://extrasane.github.io/RAILSpkg/reference/rails_design.md) | Wrap the weights in a `survey` design for `svymean()`, `svyglm()`, `svyby()` |

**Preparing inputs**

| Function | Purpose |
|----|----|
| [`rails_cells()`](https://extrasane.github.io/RAILSpkg/reference/rails_cells.md) | Aggregate microdata into covariate cells |
| [`rails_totals()`](https://extrasane.github.io/RAILSpkg/reference/rails_totals.md) | Population margins, named to match the calibration model matrix |
| [`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md) | Download ACS PUMS microdata from the Census API |
| [`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md) | Band raw ACS codes into the analysis variables |
| [`pums_reference()`](https://extrasane.github.io/RAILSpkg/reference/pums_reference.md) | Cleaned records, cell table and margins in one object |
| [`pums_variables()`](https://extrasane.github.io/RAILSpkg/reference/pums_variables.md) | The ACS variables the fetch requests |
| [`rails_simulate()`](https://extrasane.github.io/RAILSpkg/reference/rails_simulate.md) | Synthetic population for examples and testing |

**Algorithm stages** — exported for inspection;
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
calls them for you

| Function | Purpose |
|----|----|
| [`nps_weights()`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md) | Nested propensity score weights for one model ([`?nps`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md)) |
| [`gvs_step()`](https://extrasane.github.io/RAILSpkg/reference/gvs_step.md) | Greedy variable selection: likelihood-ratio test for one candidate ([`?gvs`](https://extrasane.github.io/RAILSpkg/reference/gvs_step.md)) |

**Deprecated** — the pre-package names, kept working:
[`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md),
[`fun.sub.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md),
[`fun.nps()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md),
[`fun.lkd()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md),
[`fun.out()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md).
See
[`?"RAILS-deprecated"`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md).

## Longer documentation

Full documentation, including both vignettes rendered with their output,
is at <https://extrasane.github.io/RAILSpkg/>.

| Where | What |
|----|----|
| [`vignette("rails")`](https://extrasane.github.io/RAILSpkg/articles/rails.md) | Getting started: the problem, a fit, diagnostics, variance, `survey` |
| [`vignette("subgroup-rails")`](https://extrasane.github.io/RAILSpkg/articles/subgroup-rails.md) | Subgroup calibration, choosing a stratifier, per-stratum variance |
| [`?rails`](https://extrasane.github.io/RAILSpkg/reference/rails.md) | Every argument, the benchmark-name mapping, how to reproduce the paper |

## Building the reference sample from ACS PUMS

The application uses ACS PUMS as its probability reference. Three
functions cover that end-to-end:

``` r

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

[`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md)
applies the 2022 ACS bandings for age, sex, race and ethnicity, income,
education, tenure and region, as fixed-level factors so extracts taken
at different times line up. It warns if you hand it another vintage,
because ACS category codes change between releases.

[`pums_reference()`](https://extrasane.github.io/RAILSpkg/reference/pums_reference.md)
handles the step that is easiest to get wrong by hand: records missing
any covariate cannot enter a cell, but dropping them would also drop
their share of the population, so the retained weights are rescaled to
the population total computed *before* the drop. It returns the cleaned
records, the aggregated cell table, and the margins together.

## Two notes on inference

[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
needs individual records, not cell tables – it needs one `y` per record.
The stacked form needs them for a second reason: within a covariate cell
the model matrix row, the weight and the fitted propensity are all
constant, but the outcome is not. Fits built with `aggregated = TRUE`
cannot be passed to it.

Standard errors taken off a
[`rails_design()`](https://extrasane.github.io/RAILSpkg/reference/rails_design.md)
object treat the weights as fixed, matching
[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)’s
simplified form. For the correction, use `rails_var(type = "both")`; the
`ratio_full_naive` it reports is the ratio of the two, which can fall
either side of 1 – raking to known margins removes variance, so the
corrected figure is often the smaller one.

## Migrating from the script version

The pre-package names still work and forward to their replacements with
a deprecation warning:

| Old | New |
|----|----|
| [`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md) | [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md) (its defaults reproduce it exactly) |
| [`fun.sub.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md) | [`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md) |
| [`fun.nps()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md) | [`nps_weights()`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md) |
| [`fun.lkd()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md) | [`gvs_step()`](https://extrasane.github.io/RAILSpkg/reference/gvs_step.md) |
| [`fun.out()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md) | [`weight_summary()`](https://extrasane.github.io/RAILSpkg/reference/weight_summary.md) |

[`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
still returns exactly what it always did.
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
matches it too, with one argument: the application code ran the LIFO
stage ascending (climb from the starting model, halt at the first
failure) where the manuscripts define it as descending, so
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
defaults to `lifo = "descending"` and `lifo = "ascending"` reproduces
the published application. The two agree whenever raking convergence is
monotone, which covers every case in the package’s own tests. What
changed is that the model is now specified rather than fixed: `start` is
the model taken as given, `scope` the largest model selection may reach.
`start = 2, scope = 3` is the published estimator; `scope = 2` gives a
two-way-only fit with no selection; `start = 1, scope = 3` offers the
two-way terms for selection too. Each accepts an interaction order, a
formula, or an explicit vector of terms.

## Citation

``` r

citation("RAILS")
```

## License

MIT © Huiding Chen

## The paper code lives elsewhere

This repository is the R package. The simulations, the *All of Us* /
PUMS application pipeline, and the results behind the two papers are at
[extrasane/RAILS](https://github.com/extrasane/RAILS).
