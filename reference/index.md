# Package index

## Fitting

The estimator.
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
calibrates the whole cohort;
[`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md)
runs the same procedure independently within each stratum.

- [`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md) :
  RAILS: raking with assisted nested propensity score and LIFO selection
- [`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md)
  : Subgroup RAILS

## Inspecting a fit

S3 methods on `rails_fit`, plus a diagnostic that works on any weight
vector.

- [`weights(`*`<rails_fit>`*`)`](https://extrasane.github.io/RAILSpkg/reference/weights.rails_fit.md)
  [`weights(`*`<rails_subgroup_fit>`*`)`](https://extrasane.github.io/RAILSpkg/reference/weights.rails_fit.md)
  : Extract RAILS weights
- [`summary(`*`<rails_fit>`*`)`](https://extrasane.github.io/RAILSpkg/reference/summary.rails_fit.md)
  [`print(`*`<summary.rails_fit>`*`)`](https://extrasane.github.io/RAILSpkg/reference/summary.rails_fit.md)
  : Summarize a RAILS fit
- [`plot(`*`<rails_fit>`*`)`](https://extrasane.github.io/RAILSpkg/reference/plot.rails_fit.md)
  : Plot the RAILS weight distribution
- [`weight_summary()`](https://extrasane.github.io/RAILSpkg/reference/weight_summary.md)
  : Weight diagnostics

## Inference

[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
is the reference for inference on a weighted mean.
[`rails_design()`](https://extrasane.github.io/RAILSpkg/reference/rails_design.md)
hands the weights to the survey package for everything else.

- [`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
  : Variance of a RAILS-weighted mean
- [`rails_design()`](https://extrasane.github.io/RAILSpkg/reference/rails_design.md)
  : Hand RAILS weights to the survey package

## Preparing inputs

Cell tables and population margins, from your own data or from ACS PUMS.

- [`rails_cells()`](https://extrasane.github.io/RAILSpkg/reference/rails_cells.md)
  : Aggregate microdata into covariate cells
- [`rails_totals()`](https://extrasane.github.io/RAILSpkg/reference/rails_totals.md)
  : Population margins matching a cell table
- [`pums_fetch()`](https://extrasane.github.io/RAILSpkg/reference/pums_fetch.md)
  : Download ACS PUMS microdata from the Census API
- [`pums_recode()`](https://extrasane.github.io/RAILSpkg/reference/pums_recode.md)
  : Recode raw ACS PUMS into RAILS analysis variables
- [`pums_reference()`](https://extrasane.github.io/RAILSpkg/reference/pums_reference.md)
  : Build a RAILS reference sample from recoded PUMS
- [`pums_variables()`](https://extrasane.github.io/RAILSpkg/reference/pums_variables.md)
  : ACS PUMS variables RAILS uses

## Simulation

A synthetic population with a known truth, used by every example and
vignette here.

- [`rails_simulate()`](https://extrasane.github.io/RAILSpkg/reference/rails_simulate.md)
  : Simulate a population with a biased volunteer sample

## Algorithm stages

Exported so the individual steps can be inspected or reused.
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
calls them for you.

- [`nps_weights()`](https://extrasane.github.io/RAILSpkg/reference/nps_weights.md)
  : Nested propensity score weights
- [`gvs_step()`](https://extrasane.github.io/RAILSpkg/reference/gvs_step.md)
  : Greedy variable selection: test one candidate interaction

## Deprecated

The pre-package function names. They still work and forward to their
replacements.

- [`fun.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
  [`fun.sub.rails.threeway()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
  [`fun.nps()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
  [`fun.lkd()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
  [`fun.out()`](https://extrasane.github.io/RAILSpkg/reference/RAILS-deprecated.md)
  : Deprecated functions

## Package

- [`RAILS`](https://extrasane.github.io/RAILSpkg/reference/RAILS-package.md)
  [`RAILS-package`](https://extrasane.github.io/RAILSpkg/reference/RAILS-package.md)
  : RAILS: Raking with Assisted Nested Propensity Score and LIFO
  Selection
