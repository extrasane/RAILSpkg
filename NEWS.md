# RAILS 0.2.0

## The variance estimators now run on cells

Every other part of the package already computed on aggregated covariate cells;
`rails_var()` was the exception, rebuilding individual-level model matrices and
refusing any fit made with `aggregated = TRUE`. That is fixed.

* Both variance forms are computed from cell-level sufficient statistics. Within
  a cell the model matrix row, the weight and the fitted propensity are
  constant, so each block of the sandwich needs at most the record count, the
  cell sum of the outcome, and the cell sum of its square. The reduction is
  exact; only the summation order changes.
* On a simulated cohort of 75,270 records collapsing to 164 cells, the stacked
  variance is 43x faster than before and agrees with the old record-level
  result to 4e-14 relative. The regression fixture against the original
  implementation passes unchanged.
* `rails_var()` accepts fits built with `aggregated = TRUE`. Pass the outcome
  per cell as a data frame with columns `sum` and `sumsq`; for a 0/1 outcome
  the two are the same case count.
* The stacked form no longer requires `keep_data = TRUE`. Neither form touches
  the retained records.
* `rails_cells()` gains two columns, `n` (records in the cell) and `weight_sq`
  (sum of squared weights). The sandwich meat is quadratic in the design
  weights, so a cell table carrying only `weight` cannot reconstruct it. Cell
  tables built by hand still work for the estimator and for the simplified
  variance; the stacked form needs `weight_sq` and says so.
* `rails()` retains the reference cell table on the fit as `cells_ref`.

Fits saved by 0.1.0 carry no `cells_ref` and no `weight_sq`, so the stacked
variance errors on them with an explicit message. Refit to use it.

# RAILS 0.1.0

First public release.

* `rails()` and `rails_subgroup()` replace `fun.rails.threeway()` and
  `fun.sub.rails.threeway()`. The model is now specified through `start` (the
  model taken as given) and `scope` (the largest model selection may reach)
  rather than baked into the function name. Each accepts an interaction order,
  a one-sided formula, or an explicit character vector of terms. The defaults
  reproduce the published estimator exactly.
* `rails_var()` brings the variance estimator into the package; it was
  previously only available in the simulation sources. `type` selects the
  simplified (fixed-weight) form, the stacked estimating-equation form, or both.
  The default is simplified, matching the simulation pipeline.
* `rails_cells()` and `rails_totals()` build the aggregated cell tables and
  name-matched population margins that the estimator needs, so callers no
  longer assemble them by hand.
* `rails_design()` hands RAILS weights to `survey`.
* `rails_simulate()` generates a synthetic population for examples and
  vignettes.
* `fun.rails.threeway()`, `fun.sub.rails.threeway()`, `fun.nps()`,
  `fun.lkd()` and `fun.out()` are deprecated but still work.

## Verified against the published implementation

Every estimator and the variance estimator were compared directly against the
pre-package sources on identical inputs: 35 of 35 quantities are bit-identical,
including all seven benchmark methods and all six sandwich variance components.
`tests/testthat/test-reproduces-original.R` pins this, so any future change that
alters a published number fails the test suite.

## ACS PUMS reference sample

* `pums_fetch()`, `pums_recode()` and `pums_reference()` build a reference
  sample from ACS PUMS: download, band the raw ACS codes into fixed-level
  analysis factors, then return cleaned records, the aggregated cell table and
  the population margins as one object. The recoding grammar was previously
  duplicated between two application scripts.
* `pums_recode()` warns when handed a vintage other than 2022, since ACS
  category codes change between releases.
* `pums_fetch()` reads its key from `CENSUS_API_KEY` and never takes one inline.
* `rails_simulate()` now returns `agegroup`, `sex`, `income`, `homeown` and
  `edu` with labelled levels instead of `catx1`-`catx5`, so worked examples read
  as `agegroup:sex` rather than `catx1:catx2`. The data-generating process is
  unchanged.

## Subgroup variance

* `rails_var()` now accepts a `rails_subgroup_fit` and returns one row per
  stratum. `truth` may be one value for all strata or a vector named by level.
* The default stays `type = "simplified"`, which is what the subgroup paper
  reports. `type = "stacked"` applies within each stratum and warns that the
  per-stratum sandwich carries no between-stratum covariance, so those intervals
  are usable stratum by stratum but not for pooling or differencing.
* A stratum with no converged model yields an NA row and a warning instead of
  aborting the call.
* `rails_design()` is documented as what it is: exact point estimates, and
  survey's own design-based SEs, which are close to but not identical with
  `rails_var(type = "simplified")`.
* Corrected a documentation error: the simplified variance does not reliably
  understate. Raking to known margins removes variance, so the stacked SE is
  often the smaller of the two.

## Checked against the manuscripts

* Corrected the expansion of the acronym. RAILS is **R**aking with **A**ssisted
  Nested Propens**i**ty Score and **L**IFO **S**election, not "Raking-Assisted
  Integration of Linked Surveys" -- the latter was inherited from an early
  scaffold and had reached the package Title, README and citation.
* Documentation now uses the manuscripts' terminology (NPS, GVS, LIFO), and
  `?rails` maps the benchmark columns onto the paper's estimator names.
* `rails_subgroup()` is documented as the subgroup-calibration arm of the
  global-versus-subgroup comparison, which is what the subgroup paper studies,
  rather than as a separate estimator.
* `?rails` records that LIFO here ascends and halts at the first non-convergent
  step, whereas the manuscripts describe it as descending from the full selected
  model. The implementation matches the code that produced the published
  results; the two agree when raking convergence is monotone.

## LIFO direction

* `rails()` and `rails_subgroup()` take `lifo = c("descending", "ascending")`.
  Descending is the default and the method as defined in Chen et al. (2025):
  rake the full selected model, and pop the last-added term until raking
  converges. Ascending -- climb from the starting model, halt at the first
  failure -- is what the application code did, and reproduces its results.
* `fun.rails.threeway()` and `fun.sub.rails.threeway()` are pinned to ascending,
  so existing scripts are unaffected.
* `survey`'s "Failed to converge" warning is muffled and folded into the
  package's own warning, instead of firing once per LIFO attempt.

## Names aligned with the manuscripts

* `ps_weights()` is now `nps_weights()` and `lrt_step()` is now `gvs_step()`,
  matching the nested propensity score (NPS) and greedy variable selection (GVS)
  stages the papers name. `?nps` and `?gvs` reach them as aliases.
* The deprecated `fun.nps()` and `fun.lkd()` forward to the new names, so
  existing scripts are unaffected.
