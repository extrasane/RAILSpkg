#' Extract RAILS weights
#'
#' @param object A `rails_fit` or `rails_subgroup_fit`.
#' @param ... Ignored.
#'
#' @return A numeric vector: one weight per record of the non-probability sample
#'   that was fitted, or one per cell when the fit was built from cell tables.
#'   Records dropped for missing covariates, and strata with no reference
#'   records, are `NA`.
#'
#' @examples
#' pop <- rails_simulate(2000, seed = 3)
#' ref <- pop[pop$s == 1, ]; ref$dweight <- 1 / ref$ps
#' fit <- rails(pop[pop$aou == 1, ], ref, c("agegroup", "sex"),
#'              weights_ref = "dweight", verbose = FALSE)
#' summary(weights(fit))
#'
#' @export
weights.rails_fit <- function(object, ...) object$weights

#' @rdname weights.rails_fit
#' @export
weights.rails_subgroup_fit <- function(object, ...) object$weights

#' @export
print.rails_fit <- function(x, ...) {
  cat("RAILS fit\n")
  cat("  covariates   : ", paste(x$vars, collapse = ", "), "\n", sep = "")
  cat("  start model  : ", length(x$terms_base), " term(s)\n", sep = "")
  cat("  scope        : ", length(x$terms_scope), " term(s), so ",
      length(x$terms_scope) - length(x$terms_base), " candidate(s)\n", sep = "")
  if (isTRUE(x$selectable)) {
    cat("  selected     : ", length(x$terms_selected), " term(s) at alpha = ",
        x$alpha, "\n", sep = "")
    if (length(x$terms_selected)) {
      cat("                 ", paste(x$terms_selected, collapse = ", "), "\n",
          sep = "")
    }
  } else {
    cat("  selected     : nothing to select (scope equals start)\n")
  }
  n_extra <- max(0, x$n_used - x$n_base)
  cat("  raked        : ", n_extra, " of ", length(x$terms_selected),
      " selected term(s) survived the stepwise walk\n", sep = "")
  cat("  cells        : ", nrow(x$cells), "\n", sep = "")
  w <- x$weights
  if (all(is.na(w))) {
    cat("  weights      : all NA (no model converged)\n")
  } else {
    cat("  weights      : n = ", sum(!is.na(w)), ", sum = ",
        format(sum(w, na.rm = TRUE), big.mark = ",", digits = 7),
        ", range ", format(min(w, na.rm = TRUE), digits = 3), " to ",
        format(max(w, na.rm = TRUE), digits = 3), "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a RAILS fit
#'
#' @param object A `rails_fit`.
#' @param x A `summary.rails_fit`.
#' @param ... Ignored.
#'
#' @return A `summary.rails_fit` object: a list with the selection path, the
#'   terms actually raked to, weight diagnostics from [weight_summary()], the
#'   design effect, and the effective sample size.
#'
#' @examples
#' pop <- rails_simulate(2000, seed = 3)
#' ref <- pop[pop$s == 1, ]; ref$dweight <- 1 / ref$ps
#' fit <- rails(pop[pop$aou == 1, ], ref, c("agegroup", "sex"),
#'              weights_ref = "dweight", verbose = FALSE)
#' summary(fit)
#'
#' @export
summary.rails_fit <- function(object, ...) {
  w  <- object$weights[!is.na(object$weights)]
  ## Kish's approximation: the variance inflation from unequal weights.
  deff <- if (length(w)) length(w) * sum(w^2) / sum(w)^2 else NA_real_
  out <- list(
    vars           = object$vars,
    terms_base     = object$terms_base,
    terms_scope    = object$terms_scope,
    selectable     = object$selectable,
    alpha          = object$alpha,
    terms_selected = object$terms_selected,
    terms_used     = object$terms_used,
    n_cells        = nrow(object$cells),
    converged      = object$converged,
    diagnostics    = weight_summary(w),
    deff           = deff,
    n_eff          = if (is.na(deff)) NA_real_ else length(w) / deff
  )
  class(out) <- "summary.rails_fit"
  out
}

#' @rdname summary.rails_fit
#' @export
print.summary.rails_fit <- function(x, ...) {
  cat("RAILS fit summary\n\n")
  cat("Covariates: ", paste(x$vars, collapse = ", "), "\n", sep = "")
  cat("Model     : ", length(x$terms_base), " term(s) in the starting model, ",
      length(x$terms_scope) - length(x$terms_base), " candidate(s) in scope, ",
      "alpha ", x$alpha, "\n\n", sep = "")

  cat("Selected terms (", length(x$terms_selected), "):\n", sep = "")
  if (length(x$terms_selected)) {
    cat(paste0("  ", seq_along(x$terms_selected), ". ", x$terms_selected),
        sep = "\n")
  } else {
    cat("  none\n")
  }
  cat("\nRaked to ", if (identical(x$terms_used, NA_character_)) 0L else
      length(x$terms_used), " term(s); full walk ",
      if (isTRUE(x$converged)) "completed" else "stopped early", ".\n\n", sep = "")

  cat("Weight diagnostics:\n")
  print(round(x$diagnostics, 4))
  cat("\nDesign effect (Kish): ", format(x$deff, digits = 4),
      "\nEffective sample size: ", format(x$n_eff, digits = 6), "\n", sep = "")
  invisible(x)
}

#' Plot the RAILS weight distribution
#'
#' Draws a histogram of the fitted weights on a log scale, which is where a
#' heavy right tail -- the usual failure mode of propensity weighting -- is
#' visible.
#'
#' @param x A `rails_fit`.
#' @param breaks Passed to [graphics::hist()].
#' @param ... Passed to [graphics::hist()].
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' pop <- rails_simulate(2000, seed = 3)
#' ref <- pop[pop$s == 1, ]; ref$dweight <- 1 / ref$ps
#' fit <- rails(pop[pop$aou == 1, ], ref, c("agegroup", "sex"),
#'              weights_ref = "dweight", verbose = FALSE)
#' plot(fit)
#'
#' @export
plot.rails_fit <- function(x, breaks = 40, ...) {
  w <- x$weights[!is.na(x$weights) & x$weights > 0]
  if (!length(w)) stop("No positive weights to plot.", call. = FALSE)
  graphics::hist(log10(w), breaks = breaks,
                 main = "RAILS weights", xlab = "log10(weight)", ...)
  graphics::abline(v = log10(stats::median(w)), lty = 2)
  invisible(x)
}

#' @export
print.rails_subgroup_fit <- function(x, ...) {
  cat("Subgroup RAILS fit\n")
  cat("  stratifier : ", x$by, " (", length(x$fits), " of ",
      length(x$levels), " level(s) fitted)\n", sep = "")
  cat("  covariates : ", paste(x$vars, collapse = ", "), "\n\n", sep = "")
  for (nm in names(x$fits)) {
    f <- x$fits[[nm]]
    cat("  ", x$by, " = ", nm, ": ", length(f$terms_selected),
        " selected, ", max(0, f$n_used - f$n_base), " raked, ",
        nrow(f$cells), " cells\n", sep = "")
  }
  if (!is.null(x$weights)) {
    w <- x$weights
    cat("\n  weights: n = ", sum(!is.na(w)), ", sum = ",
        format(sum(w, na.rm = TRUE), big.mark = ",", digits = 7), "\n", sep = "")
  }
  invisible(x)
}
