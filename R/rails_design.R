#' Hand RAILS weights to the survey package
#'
#' Wraps a fit's weights in a [survey::svydesign()] so the usual survey
#' machinery -- `svymean()`, `svytotal()`, `svyglm()`, `svyby()` -- can be used
#' on the weighted sample.
#'
#' @param fit A `rails_fit` from [rails()], or a `rails_subgroup_fit` from
#'   [rails_subgroup()].
#' @param data The data frame the weights belong to: the same non-probability
#'   sample passed to [rails()], with whatever outcome columns you want to
#'   analyse. Must have as many rows as the fit has weights.
#' @param ... Passed to [survey::svydesign()].
#'
#' @return A `survey.design` object with the RAILS weights attached and a
#'   `rails_weights` column added to its data.
#'
#' @section What this design does and does not carry:
#' The weights are exact: point estimates from this object -- `svymean()`,
#' `svytotal()`, `svyglm()`, `svyby()` -- are the RAILS estimates.
#'
#' Standard errors are a different matter. The survey package computes its own
#' design-based variance from the weights it is handed, treating them as fixed.
#' That is close to, but not the same estimator as, [rails_var()] with
#' `type = "simplified"`, which is computed independently; in practice the two
#' agree to about four significant figures, differing by little more than an
#' n/(n-1) factor. Neither accounts for the weights having been estimated.
#'
#' So: use this design freely for point estimates and for model fitting, and use
#' [rails_var()] when you want inference on a mean that reflects how the weights
#' were built. Do not expect the two standard errors to match to the last digit --
#' they are different estimators of the same fixed-weight quantity.
#'
#' @seealso [rails_var()] for corrected inference on a weighted mean.
#'
#' @examples
#' pop <- rails_simulate(20000, seed = 42)
#' aou <- pop[pop$aou == 1, ]
#' ref <- pop[pop$s == 1, ]
#' ref$dweight <- 1 / ref$ps
#'
#' fit <- rails(aou, ref, vars = c("agegroup", "sex", "income"),
#'              weights_ref = "dweight", verbose = FALSE)
#' des <- rails_design(fit, aou)
#' survey::svymean(~y, des)
#'
#' @export
rails_design <- function(fit, data, ...) {
  w <- if (inherits(fit, "rails_subgroup_fit")) fit$weights else fit$weights

  if (is.null(w)) {
    stop("This fit has no per-record weights (it was built from cell tables). ",
         "Refit from microdata to build a survey design.", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (nrow(data) != length(w)) {
    stop("`data` has ", nrow(data), " rows but the fit carries ", length(w),
         " weights. Pass the same non-probability sample you fitted on.",
         call. = FALSE)
  }
  if (anyNA(w)) {
    warning("rails_design(): ", sum(is.na(w)), " record(s) have NA weights and ",
            "will contribute nothing to survey estimates.", call. = FALSE)
  }

  data <- as.data.frame(data)
  data$rails_weights <- w
  survey::svydesign(id = ~1, weights = ~rails_weights, data = data, ...)
}
