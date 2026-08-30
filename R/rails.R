#' RAILS: raking with assisted nested propensity score and LIFO selection
#'
#' Combines a non-probability sample with a probability reference sample to
#' produce population-representative weights. The procedure is:
#'
#' \enumerate{
#'   \item Fit the \strong{nested propensity score} (NPS) model on the starting
#'     model `start`, giving base weights and a baseline log-likelihood.
#'   \item \strong{Greedy variable selection} (GVS): a forward likelihood-ratio
#'     search over the terms in `scope` but not in `start` (via [gvs_step()]) at
#'     significance level `alpha`, adding whichever significant candidate has the
#'     largest average increment in the log pseudo-likelihood.
#'   \item \strong{LIFO stepwise raking}: rake the full selected model to the
#'     population margins. If raking does not converge, pop the last-added term,
#'     refit NPS, and rake again -- last in, first out -- until it converges or
#'     no selected term is left.
#' }
#'
#' Those three stages are what the method is named for in Chen et al. (2025):
#' \strong{R}aking with \strong{A}ssisted Nested Propens\strong{i}ty Score and
#' \strong{L}IFO \strong{S}election.
#'
#'
#' All computation runs on aggregated covariate cells, so cost scales with the
#' number of distinct covariate patterns rather than the sample size.
#'
#' @param np The non-probability sample (e.g. a volunteer biobank): a data frame
#'   of individual records, or a cell table if `aggregated = TRUE`.
#' @param ref The probability reference sample, same form as `np`.
#' @param vars Character vector of covariate names shared by both samples. These
#'   are the main effects; interactions are generated from them.
#' @param weights_ref Design weights for `ref`: a numeric vector or the name of
#'   a column. Almost always required, since a probability sample's design
#'   weights are what make its margins represent the population. Ignored when
#'   `aggregated = TRUE`.
#' @param weights_np Optional weights for `np`. `NULL` counts each record as 1,
#'   which is the usual choice for a biobank. Ignored when `aggregated = TRUE`.
#' @param pop_totals Named numeric vector of population margins. `NULL` (the
#'   default) computes them from `ref` via [rails_totals()]. Supply your own
#'   when the margins come from an external source (a census table, say) rather
#'   than from the reference sample itself.
#' @param aggregated Set `TRUE` when `np` and `ref` are already cell tables with
#'   a `weight` column. The row-to-cell map is then unavailable, so `weights()`
#'   returns one weight per cell rather than per record, and [rails_var()] will
#'   not run. See the note below.
#' @param start The starting model: every term in it is included without being
#'   tested. Give an interaction order (`2`, the default, means all main effects
#'   and all two-way interactions of `vars`), a one-sided formula, or a character
#'   vector of terms such as `c("age", "sex", "age:sex")`.
#' @param scope The largest model selection may reach, in the same three forms.
#'   `3` (the default) means up to three-way interactions. Everything in `scope`
#'   but not in `start` is a candidate for the forward test, so `start = 2,
#'   scope = 3` offers the three-way terms, while `start = 1, scope = 3` offers
#'   the two-way *and* three-way terms. `scope` must contain `start`; when the
#'   two are equal there is nothing to select and the starting model is simply
#'   raked. Term components may be written in any order -- `"sex:age"` matches
#'   `"age:sex"`.
#' @param alpha Significance threshold for the forward likelihood-ratio test.
#' @param nsiz Population total the weights are scaled to. `NULL` uses the total
#'   reference weight, which is the right choice when `ref` carries design
#'   weights that already sum to the population.
#' @param benchmarks Also compute the comparison methods used in the papers
#'   (unweighted, raking-only, propensity-only, propensity-plus-raking), placed
#'   in `fit$benchmarks` as per-record weights, on the same scale as
#'   [weights()]. Off by default: this is simulation-study output, not something
#'   you need in order to produce weights.
#' @param fallback When selection picks nothing, or the first stepwise step
#'   fails to converge, fall back to the raked starting model instead of
#'   returning `NA` weights. `FALSE` is the default so that `rails()` reproduces
#'   the published estimator exactly; `TRUE` is often more useful in applied
#'   work. Either way the outcome is recorded in `converged`. This has no effect
#'   when `scope` equals `start`: there the starting model is what was asked for,
#'   not a fallback, and is always used.
#' @param lifo Direction of the stepwise raking stage. `"descending"` (the
#'   default) is the LIFO strategy of Chen et al. (2025): start at the full
#'   selected model and pop the last-added term until raking converges.
#'   `"ascending"` instead climbs from the starting model and halts at the first
#'   non-convergent step, which is what the application code did; use it to
#'   reproduce published application results exactly. The two agree whenever
#'   raking convergence is monotone in the number of added terms; where it is
#'   not, `"descending"` keeps more terms.
#' @param keep_data Retain the `vars` columns of both samples inside the fit, so
#'   [rails_var()] can rebuild individual-level model matrices later. Costs
#'   memory proportional to the input; set `FALSE` if you only need weights.
#' @param verbose Report selection and stepwise progress via [message()].
#'
#' @return An object of class `rails_fit`: a list with `weights` (one per input
#'   record, or per cell when `aggregated = TRUE`), `cells`, `weights_cell` (cell
#'   totals), `terms_base`, `terms_scope`, `terms_selected`, `terms_used`,
#'   `formula_used`, `theta`, `pop_totals`, `n_used`, `converged`, `benchmarks`,
#'   and `call`. Use [weights()], [summary()] and [plot()] on it, and
#'   [rails_design()] to hand it to \pkg{survey}.
#'
#' @section Benchmark names:
#' With `benchmarks = TRUE`, `fit$benchmarks` carries the comparison methods.
#' They correspond to the estimators tabulated in Chen et al. (2025):
#'
#' | package               | paper       | description                                     |
#' |-----------------------|-------------|-------------------------------------------------|
#' | `d_unweighted`        | `naive`     | the unweighted cohort mean                      |
#' | `d_cal1`              | `cal-1`     | raking to univariate margins from equal weights  |
#' | `d_cal2`              | `cal-2`     | raking on main effects and all two-way terms     |
#' | `d_nps1`              | `nps-1`     | NPS weights, main effects only                   |
#' | `d_nps2`              | `nps-2`     | NPS weights through the starting model           |
#' | `d_nps1_rake`         | `nps-cal-1` | NPS main effects, then univariate raking         |
#' | `d_nps2_rake`         |             | NPS through the starting model, then raking      |
#'
#' `cal-2` is the one the paper reports as failing to converge in over 99 percent
#' of simulations; here it returns `NA` with a warning when raking fails, so a
#' non-convergent benchmark is visible rather than silently absent.
#'
#' @section Reproducing the published estimator:
#' The defaults -- `start = 2`, `scope = 3`, `fallback = FALSE` -- are exactly
#' the procedure of `fun.rails.threeway()`, so `rails()` and the deprecated shim
#' give identical weights on identical input. Changing `start` or `scope` departs
#' from the paper deliberately; `fallback = TRUE` departs from it only in the
#' degenerate case where the published code would have returned `NA`.
#'
#' @section Aggregated input:
#' `aggregated = TRUE` accepts cell tables directly, which is what the original
#' `fun.rails.threeway()` required. It is the right choice when the microdata
#' cannot be materialized in one frame. The cost is that the fit no longer knows
#' which record belongs to which cell, so per-record weights and the variance
#' estimator are unavailable.
#'
#' @seealso [rails_subgroup()] to run within strata, [rails_var()] for standard
#'   errors, [rails_cells()] and [rails_totals()] to build inputs by hand.
#'
#' @examples
#' pop <- rails_simulate(20000, seed = 42)
#' aou <- pop[pop$aou == 1, ]
#' ref <- pop[pop$s == 1, ]
#' ref$dweight <- 1 / ref$ps
#'
#' fit <- rails(aou, ref, vars = c("agegroup", "sex", "income"),
#'              weights_ref = "dweight", verbose = FALSE)
#' fit
#'
#' # Weighted prevalence of the outcome, against the population truth:
#' w <- weights(fit)
#' sum(w * aou$y) / sum(w)
#' mean(pop$y)
#'
#' @export
rails <- function(np, ref, vars,
                  weights_ref = NULL,
                  weights_np  = NULL,
                  pop_totals  = NULL,
                  aggregated  = FALSE,
                  start       = 2,
                  scope       = 3,
                  alpha       = 0.05,
                  nsiz        = NULL,
                  benchmarks  = FALSE,
                  fallback    = FALSE,
                  lifo        = c("descending", "ascending"),
                  keep_data   = TRUE,
                  verbose     = TRUE) {

  lifo <- match.arg(lifo)

  cl <- match.call()

  if (length(vars) < 2) {
    stop("`vars` needs at least two variables for interactions to exist.",
         call. = FALSE)
  }

  if (isTRUE(aggregated)) {
    check_cells(np,  vars, "`np`")
    check_cells(ref, vars, "`ref`")
    cells_np  <- as.data.frame(np)
    cells_ref <- as.data.frame(ref)
    row_cell  <- NULL
    kept      <- NULL
  } else {
    if (is.null(weights_ref)) {
      warning("`weights_ref` is NULL: every reference record counts as 1. ",
              "Pass the design weights unless `ref` is genuinely unweighted.",
              call. = FALSE)
    }
    cells_np  <- rails_cells(np,  vars, weights = weights_np)
    cells_ref <- rails_cells(ref, vars, weights = weights_ref)
    row_cell  <- attr(cells_np, "row_cell")
    kept <- if (isTRUE(keep_data)) {
      d_ref <- if (is.null(weights_ref)) {
        rep(1, nrow(ref))
      } else if (is.character(weights_ref)) {
        as.numeric(ref[[weights_ref]])
      } else {
        as.numeric(weights_ref)
      }
      list(np    = as.data.frame(np)[, vars, drop = FALSE],
           ref   = as.data.frame(ref)[, vars, drop = FALSE],
           d_ref = d_ref)
    } else NULL
    cells_np  <- as.data.frame(cells_np)
    cells_ref <- as.data.frame(cells_ref)
  }

  if (is.null(nsiz)) nsiz <- sum(cells_ref$weight)
  if (is.null(pop_totals)) {
    ## Margins are needed up to the highest order the walk can reach.
    need <- max(max_term_order(resolve_model(start, vars, "start")),
                max_term_order(resolve_model(scope, vars, "scope")))
    pop_totals <- rails_totals(cells_ref, vars = vars, order = need)
  }

  eng <- rails_engine(
    cells_np = cells_np, cells_ref = cells_ref, pop_totals = pop_totals,
    vars = vars, start = start, scope = scope,
    alpha = alpha, nsiz = nsiz, benchmarks = benchmarks,
    fallback = fallback, lifo = lifo, verbose = verbose
  )

  ## Cell totals -> per-individual weights, then spread back over records.
  w_cell_ind <- eng$weights_cell / cells_np$weight
  w_out <- if (is.null(row_cell)) w_cell_ind else w_cell_ind[row_cell]

  fit <- c(eng, list(
    weights      = w_out,
    weights_unit = w_cell_ind,
    aggregated   = isTRUE(aggregated),
    row_cell     = row_cell,
    data         = kept,
    n_np         = if (is.null(row_cell)) NA_integer_ else length(row_cell),
    call         = cl
  ))
  class(fit) <- "rails_fit"
  fit
}
