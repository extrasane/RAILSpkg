# RAILS: Raking with Assisted Nested Propensity Score and LIFO Selection

Combines a non-probability sample (e.g. a volunteer biobank) with a
probability reference sample to produce population-representative
weights. The core estimator fits a pseudo-likelihood propensity model on
all interactions up to a base order, selects higher-order interactions
by a forward likelihood-ratio test, and applies last-in-first-out
stepwise raking to reference population margins. All computation runs on
aggregated covariate cells. Variance is available from a stacked
estimating-equation sandwich that accounts for uncertainty in the
estimated weights, and a subgroup wrapper runs the procedure within
levels of a stratifying variable.

## See also

Useful links:

- <https://github.com/extrasane/RAILSpkg>

- <https://extrasane.github.io/RAILSpkg/>

- Report bugs at <https://github.com/extrasane/RAILSpkg/issues>

## Author

**Maintainer**: Huiding Chen <extrasane@gmail.com>
([ORCID](https://orcid.org/0000-0002-5102-434X))
