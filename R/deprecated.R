#' Deprecated functions
#'
#' The pre-package function names, kept so that existing analysis scripts keep
#' running. Each forwards to its replacement and warns once per session.
#'
#' @details
#' \describe{
#'   \item{`fun.rails.threeway()`}{use [rails()], whose defaults (`start = 2`,
#'     `scope = 3`, `fallback = FALSE`) reproduce it exactly}
#'   \item{`fun.sub.rails.threeway()`}{use [rails_subgroup()]}
#'   \item{`fun.nps()`}{use [nps_weights()]}
#'   \item{`fun.lkd()`}{use [gvs_step()]}
#'   \item{`fun.out()`}{use [weight_summary()]}
#' }
#'
#' The two estimator shims return the old wide data frame -- `dt_agg_aou` with
#' `d_rails`, `d_cal1`, `d_cal2`, `d_nps1`, `d_nps2`, `d_nps1_rake`,
#' `d_nps2_rake`, `d_unweighted`, `selected_terms` and `calibrated_terms`
#' appended -- not a `rails_fit`. The weights themselves are identical to what
#' [rails()] produces on the same input under its defaults; only the container
#' differs.
#'
#' @name RAILS-deprecated
#' @keywords internal
NULL

#' @rdname RAILS-deprecated
#'
#' @param dt_agg_aou Aggregated non-probability cell counts (column `weight`).
#' @param dt_agg_pums Aggregated reference-sample cell counts (column `weight`).
#' @param pop_totals Named numeric vector of population margins.
#' @param names_univar Character vector of main-effect variables.
#' @param alpha Significance threshold for forward selection.
#' @param nsiz Population total the weights are scaled to.
#'
#' @return For `fun.rails.threeway()` and `fun.sub.rails.threeway()`, the wide
#'   data frame described above.
#'
#' @export
fun.rails.threeway <- function(
    dt_agg_aou,
    dt_agg_pums,
    pop_totals,
    names_univar = c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
    alpha        = 0.05,
    nsiz         = sum(dt_agg_pums$weight)) {

  .Deprecated("rails", package = "RAILS",
              msg = paste("fun.rails.threeway() is deprecated; use",
                          "rails() -- its defaults reproduce this exactly."))

  eng <- rails_engine(
    cells_np = as.data.frame(dt_agg_aou), cells_ref = as.data.frame(dt_agg_pums),
    pop_totals = pop_totals, vars = names_univar,
    start = 2, scope = 3, alpha = alpha, nsiz = nsiz,
    benchmarks = TRUE, fallback = FALSE, lifo = "ascending", verbose = TRUE
  )
  engine_to_wide(eng, dt_agg_aou)
}

#' Reshape an engine result into the pre-package wide data frame
#'
#' @keywords internal
#' @noRd
engine_to_wide <- function(eng, dt_agg_aou) {
  b  <- eng$benchmarks
  wt <- eng$cells$weight   # only d_rails still needs dividing down
  terms_all <- c(eng$terms_base, eng$terms_selected)

  dt_agg_aou %>%
    mutate(
      d_unweighted     = b$d_unweighted,
      d_cal1           = b$d_cal1,
      d_cal2           = b$d_cal2,
      d_nps1           = b$d_nps1,
      d_nps2           = b$d_nps2,
      d_nps1_rake      = b$d_nps1_rake,
      d_nps2_rake      = b$d_nps2_rake,
      d_rails          = eng$weights_cell / wt,
      selected_terms   = paste(terms_all, collapse = " + "),
      calibrated_terms = if (eng$n_used > 0) {
        paste(terms_all[seq_len(eng$n_used)], collapse = " + ")
      } else {
        NA_character_
      }
    )
}

#' @rdname RAILS-deprecated
#'
#' @param subgroup_var Name of the stratifying variable.
#'
#' @export
fun.sub.rails.threeway <- function(
    dt_agg_aou,
    dt_agg_pums,
    subgroup_var = "region",
    names_univar = setdiff(
      c("agegroup", "edu", "homeown", "income", "race_eth", "sex", "region"),
      subgroup_var
    ),
    alpha = 0.05) {

  .Deprecated("rails_subgroup", package = "RAILS",
              msg = "fun.sub.rails.threeway() is deprecated; use rails_subgroup().")

  if (subgroup_var %in% names_univar) {
    stop("subgroup_var must not appear in names_univar", call. = FALSE)
  }
  if (!subgroup_var %in% names(dt_agg_aou) || !subgroup_var %in% names(dt_agg_pums)) {
    stop("subgroup_var must be a column of both cell tables", call. = FALSE)
  }

  subgroup_levels <- sort(unique(dt_agg_aou[[subgroup_var]]))

  out <- lapply(subgroup_levels, function(lev) {
    message("=== Subgroup ", subgroup_var, " = ", lev, " ===")

    agg_aou_sub  <- dt_agg_aou  %>% filter(.data[[subgroup_var]] == lev)
    agg_pums_sub <- dt_agg_pums %>% filter(.data[[subgroup_var]] == lev)

    if (nrow(agg_pums_sub) == 0) {
      warning("no reference cells for ", subgroup_var, " = ", lev, " -- skipped.",
              call. = FALSE)
      return(NULL)
    }

    eng <- rails_engine(
      cells_np   = as.data.frame(agg_aou_sub),
      cells_ref  = as.data.frame(agg_pums_sub),
      pop_totals = rails_totals(as.data.frame(agg_pums_sub), vars = names_univar,
                                order = 3),
      vars = names_univar, start = 2, scope = 3, alpha = alpha,
      nsiz = sum(agg_pums_sub$weight),
      benchmarks = TRUE, fallback = FALSE, lifo = "ascending", verbose = TRUE
    )

    engine_to_wide(eng, agg_aou_sub) %>% mutate(subgroup_run = lev)
  })

  bind_rows(out)
}

#' @rdname RAILS-deprecated
#'
#' @param cat_temp A one-sided model formula.
#' @param dt_aou Aggregated non-probability cell counts.
#' @param dt_s Aggregated reference-sample cell counts.
#' @param theta_init Optional warm-start coefficient vector.
#'
#' @export
fun.nps <- function(cat_temp, dt_aou, dt_s, nsiz = sum(dt_s$weight),
                    theta_init = NULL) {
  .Deprecated("nps_weights", package = "RAILS")
  nps_weights(cat_temp, dt_aou, dt_s, nsiz = nsiz, theta_init = theta_init)
}

#' @rdname RAILS-deprecated
#'
#' @param x Candidate term.
#' @param LKHD Log-likelihood of the current model.
#' @param m_s Model matrix of the current model on the reference cells.
#'
#' @export
fun.lkd <- function(x, names_univar, LKHD, dt_s, dt_aou, m_s) {
  .Deprecated("gvs_step", package = "RAILS")
  gvs_step(x, terms_current = names_univar, loglik = LKHD,
           cells_ref = dt_s, cells_np = dt_aou, m_ref = m_s)
}

#' @rdname RAILS-deprecated
#'
#' @param w Numeric vector of weights.
#' @param ... Ignored.
#'
#' @export
fun.out <- function(w, ...) {
  .Deprecated("weight_summary", package = "RAILS")
  weight_summary(w, ...)
}
