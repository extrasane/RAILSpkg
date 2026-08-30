#' Subgroup RAILS
#'
#' Runs [rails()] independently within each level of a stratifying variable.
#' Each subgroup gets its own NPS model, its own selection path, its own stepwise
#' walk, and its own population margins computed from that subgroup's reference
#' records, with weights scaled to the subgroup total.
#'
#' This is the *subgroup calibration* strategy -- "calibration-in-sub" --
#' that Chen et al. compare against *global calibration*, in which one set of
#' weights is built from the whole cohort by [rails()] and then simply applied
#' within a subgroup. Which strategy does better is an empirical question rather
#' than a settled one; that comparison is what the subgroup paper is about, and
#' heterogeneity of the non-probability selection mechanism across subgroups is
#' the factor it identifies as most influential. This function provides the
#' subgroup arm; [rails()] provides the global arm.
#'
#' @param np,ref,vars,weights_ref,weights_np As in [rails()].
#' @param by Name of the stratifying variable. It must be a column of both
#'   samples, and is automatically excluded from `vars`.
#' @param aggregated,start,scope,alpha,benchmarks,fallback,lifo,keep_data,verbose
#'   As in [rails()], applied within every subgroup.
#' @param ... Further arguments passed to [rails()].
#'
#' @return An object of class `rails_subgroup_fit`: a list with `fits` (one
#'   `rails_fit` per level, named by level), `weights` (a single vector aligned
#'   to the rows of `np`, pooling every subgroup), `by`, `levels`, and `call`.
#'
#' @section Empty strata:
#' Levels with no reference records are skipped with a warning and their
#' non-probability records receive `NA` weights. Factor levels are never
#' dropped, so model matrices stay aligned between the two samples within each
#' stratum.
#'
#' @seealso [rails()].
#'
#' @examples
#' pop <- rails_simulate(20000, seed = 7)
#' aou <- pop[pop$aou == 1, ]
#' ref <- pop[pop$s == 1, ]
#' ref$dweight <- 1 / ref$ps
#'
#' fit <- rails_subgroup(aou, ref, vars = c("agegroup", "income"),
#'                       by = "sex", weights_ref = "dweight",
#'                       verbose = FALSE)
#' fit
#'
#' @export
rails_subgroup <- function(np, ref, vars, by,
                           weights_ref = NULL,
                           weights_np  = NULL,
                           aggregated  = FALSE,
                           start        = 2,
                           scope        = 3,
                           alpha        = 0.05,
                           benchmarks   = FALSE,
                           fallback     = FALSE,
                           lifo         = c("descending", "ascending"),
                           keep_data    = TRUE,
                           verbose      = TRUE,
                           ...) {

  lifo <- match.arg(lifo)

  cl <- match.call()

  if (by %in% vars) {
    stop("`by` must not also appear in `vars`; it is the stratifier, not a ",
         "covariate within strata.", call. = FALSE)
  }
  if (!by %in% names(np) || !by %in% names(ref)) {
    stop("`by` must be a column of both `np` and `ref`.", call. = FALSE)
  }

  np  <- as.data.frame(np)
  ref <- as.data.frame(ref)

  ## Reference design weights, resolved once so subsets stay aligned.
  d_ref <- if (isTRUE(aggregated) || is.null(weights_ref)) {
    NULL
  } else if (is.character(weights_ref) && length(weights_ref) == 1L) {
    as.numeric(ref[[weights_ref]])
  } else {
    as.numeric(weights_ref)
  }
  d_np <- if (isTRUE(aggregated) || is.null(weights_np)) {
    NULL
  } else if (is.character(weights_np) && length(weights_np) == 1L) {
    as.numeric(np[[weights_np]])
  } else {
    as.numeric(weights_np)
  }

  levs <- sort(unique(np[[by]]))
  w_all <- rep(NA_real_, nrow(np))
  fits  <- list()

  for (lev in levs) {
    say(verbose, "subgroup ", by, " = ", lev)

    i_np  <- which(np[[by]]  == lev)
    i_ref <- which(ref[[by]] == lev)

    if (length(i_ref) == 0) {
      warning("rails_subgroup(): no reference records for ", by, " = ", lev,
              " -- that stratum is skipped and its weights are NA.",
              call. = FALSE)
      next
    }

    fit_lev <- rails(
      np  = np[i_np, , drop = FALSE],
      ref = ref[i_ref, , drop = FALSE],
      vars = vars,
      weights_ref = if (is.null(d_ref)) weights_ref else d_ref[i_ref],
      weights_np  = if (is.null(d_np))  weights_np  else d_np[i_np],
      pop_totals  = NULL,
      aggregated  = aggregated,
      start       = start,
      scope       = scope,
      alpha        = alpha,
      nsiz         = NULL,
      benchmarks   = benchmarks,
      fallback     = fallback,
      lifo         = lifo,
      keep_data    = keep_data,
      verbose      = verbose,
      ...
    )

    ## Which rows of `np` this stratum owns. rails_var() needs it to line an
    ## outcome vector up with each stratum's fit.
    fit_lev$subgroup_rows <- i_np

    fits[[as.character(lev)]] <- fit_lev
    if (!isTRUE(aggregated)) w_all[i_np] <- fit_lev$weights
  }

  out <- list(
    fits    = fits,
    weights = if (isTRUE(aggregated)) NULL else w_all,
    by      = by,
    levels  = levs,
    vars    = vars,
    call    = cl
  )
  class(out) <- "rails_subgroup_fit"
  out
}
