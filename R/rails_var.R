#' Variance of a RAILS-weighted mean
#'
#' Standard error and confidence interval for a weighted mean computed with
#' RAILS weights, in either of two forms.
#'
#' The **simplified** variance (the default) treats the weights as fixed. It is
#' what the simulation pipeline reports and what `svymean()` on a plain design
#' gives, and it is cheap: no model matrix is involved.
#'
#' The **stacked** variance accounts for the weights having been estimated,
#' through the sandwich of three stacked estimating equations -- (1) the
#' pseudo-likelihood propensity score, (2) the calibration/raking constraints,
#' and (3) the weighted mean itself. Because the simplified form ignores the
#' first two, the two forms differ; `type = "both"` reports the pair along with
#' `ratio_full_naive`, the ratio of stacked to simplified.
#'
#' The direction is not fixed. Ignoring the propensity stage tends to understate,
#' but raking to known population margins removes variance the way a regression
#' estimator does, so the stacked form is often the *smaller* of the two -- in the
#' examples below the ratio is around 0.97. Report the stacked figure because it
#' is the one that accounts for how the weights were built, not because it is
#' conservative.
#'
#' @param fit A `rails_fit` from [rails()], or a `rails_subgroup_fit` from
#'   [rails_subgroup()]. Fits built with `aggregated = TRUE` are supported; see
#'   `y` and the section below.
#' @param y The outcome. For a fit built from microdata, a vector aligned to the
#'   rows of the non-probability sample passed to [rails()]. For a fit built
#'   from cell tables, a data frame with columns `sum` and `sumsq` -- the
#'   within-cell sum of the outcome and of its square -- one row per cell, in
#'   the order of the cell table. For a 0/1 outcome the two columns are equal,
#'   so `sumsq` is just the case count again. Binary or continuous.
#' @param type Which variance to compute: `"simplified"` (the default),
#'   `"stacked"`, or `"both"`.
#' @param truth Optional known population value. When supplied, the returned
#'   vector includes a coverage indicator for each interval computed, which is
#'   what the simulation studies use. For a subgroup fit, either one value for
#'   every stratum or a vector named by stratum level.
#' @param level Confidence level.
#'
#' @return For a single fit, a named numeric vector: `mu_hat`, the variance,
#'   standard error and interval bounds for the requested form, a coverage
#'   indicator, and the sample sizes. `"stacked"` and `"both"` add the six
#'   sandwich components (`term_11` through `term_33`); `"both"` also adds
#'   `ratio_full_naive`. Simplified quantities are named `*_naive`, stacked ones
#'   `*_full`.
#'
#'   For a subgroup fit, a data frame with one row per stratum: the stratifier
#'   in the first column, then those same quantities.
#'
#' @section Subgroup fits:
#' Each stratum is a separate fit -- its own propensity model, its own selected
#' terms, its own margins -- so the variance is computed within each and returned
#' one row per stratum.
#'
#' The subgroup paper reports the **simplified** variance by default, because
#' residual bias persists there and a tighter interval would overstate what the
#' estimator delivers. The stacked form is available and applies within a stratum,
#' but the per-stratum sandwich carries no covariance between strata, so it warns:
#' use those intervals stratum by stratum, and not for anything that pools or
#' differences across strata.
#'
#' @section Everything runs on cells:
#' The variance is computed from the same aggregated cells the estimator uses,
#' not from individual records. Within a covariate cell the model matrix row,
#' the weight and the fitted propensity are all constant, so each block of the
#' sandwich needs at most three numbers per cell: the record count, the sum of
#' the outcome, and the sum of its square. Blocks free of the outcome need only
#' the count; blocks linear in it need the sum; the one quadratic block needs
#' the sum of squares. On the reference side the meat is quadratic in the design
#' weights, so [rails_cells()] records `weight_sq` alongside `weight`.
#'
#' Cost therefore scales with the number of cells rather than the number of
#' records, which for a biobank is a difference of two or three orders of
#' magnitude. Given microdata, `rails_var()` reduces it to these statistics
#' itself; given cell tables, supply them through `y`. `keep_data = TRUE` is no
#' longer required for the stacked form.
#'
#' @seealso [rails()], [rails_design()].
#'
#' @examples
#' pop <- rails_simulate(20000, seed = 42)
#' aou <- pop[pop$aou == 1, ]
#' ref <- pop[pop$s == 1, ]
#' ref$dweight <- 1 / ref$ps
#'
#' fit <- rails(aou, ref, vars = c("agegroup", "sex", "income"),
#'              weights_ref = "dweight", verbose = FALSE)
#'
#' # Simplified, the default:
#' rails_var(fit, aou$y, truth = mean(pop$y))
#'
#' # Both, to compare the two:
#' v <- rails_var(fit, aou$y, type = "both", truth = mean(pop$y))
#' v[c("mu_hat", "se_naive", "se_full", "ratio_full_naive")]
#'
#' # One row per stratum for a subgroup fit:
#' sfit <- rails_subgroup(aou, ref, c("agegroup", "income"), by = "sex",
#'                        weights_ref = "dweight", start = 1, scope = 2,
#'                        verbose = FALSE)
#' rails_var(sfit, aou$y)
#'
#' @export
rails_var <- function(fit, y, type = c("simplified", "stacked", "both"),
                      truth = NULL, level = 0.95) {

  type <- match.arg(type)

  if (inherits(fit, "rails_subgroup_fit")) {
    return(rails_var_subgroup(fit, y, type, truth, level))
  }
  rails_var_one(fit, y, type, truth, level)
}

#' Variance for a single (non-subgroup) fit
#'
#' @keywords internal
#' @noRd
rails_var_one <- function(fit, y, type, truth, level) {

  want_stacked <- type %in% c("stacked", "both")

  if (!inherits(fit, "rails_fit")) {
    stop("`fit` must be a rails_fit object from rails().", call. = FALSE)
  }
  if (want_stacked && (is.null(fit$formula_used) || is.null(fit$theta))) {
    stop("`fit` has no converged model: no stepwise step succeeded, so there ",
         "are no weights to compute a variance for.", call. = FALSE)
  }

  cells  <- fit$cells
  w_cell <- as.vector(fit$weights_unit)
  n_c    <- nrow(cells)

  ## Records per cell. rails_cells() records this; a cell table built by hand
  ## may not, in which case `weight` is the count -- true whenever the
  ## non-probability sample is unweighted, which is the biobank case.
  m_cell <- if (!is.null(cells$n)) as.numeric(cells$n) else as.numeric(cells$weight)

  ## ---- Outcome, reduced to per-cell sufficient statistics ----------------
  if (is.data.frame(y) || is.matrix(y)) {
    yc <- as.data.frame(y)
    if (!all(c("sum", "sumsq") %in% names(yc))) {
      stop("Cell-level `y` needs columns `sum` and `sumsq`: the within-cell ",
           "sum of the outcome and of its square. For a 0/1 outcome the two ",
           "are equal.", call. = FALSE)
    }
    if (nrow(yc) != n_c) {
      stop("Cell-level `y` has ", nrow(yc), " rows but the fit has ", n_c,
           " cells.", call. = FALSE)
    }
    S_cell <- as.numeric(yc$sum)
    Q_cell <- as.numeric(yc$sumsq)
  } else {
    y  <- as.vector(y)
    rc <- fit$row_cell
    if (is.null(rc)) {
      stop("`fit` was built from cell tables, so `y` must be per cell too: ",
           "pass a data frame with columns `sum` and `sumsq`, the within-cell ",
           "sum of the outcome and of its square, one row per cell in the ",
           "order of `fit$cells`. For a 0/1 outcome the two are equal.",
           call. = FALSE)
    }
    if (length(y) != length(rc)) {
      stop("`y` has length ", length(y), " but the non-probability sample has ",
           length(rc), " rows.", call. = FALSE)
    }
    ok <- !is.na(rc) & !is.na(y)
    if (any(!ok)) {
      warning("rails_var(): dropping ", sum(!ok),
              " record(s) with a missing weight or outcome.", call. = FALSE)
    }
    into_cells <- function(v, g) {
      out <- numeric(n_c)
      tab <- rowsum(v, g)
      out[as.integer(rownames(tab))] <- tab[, 1]
      out
    }
    S_cell <- into_cells(y[ok], rc[ok])
    Q_cell <- into_cells(y[ok]^2, rc[ok])
    m_cell <- into_cells(rep(1, sum(ok)), rc[ok])
  }

  ## A cell with an unusable weight carries no information and drops out whole.
  keep <- !is.na(w_cell) & m_cell > 0
  if (!any(keep)) {
    stop("`fit` has no usable weights, so there is no variance to compute.",
         call. = FALSE)
  }

  m <- m_cell[keep]
  w <- w_cell[keep]
  S <- S_cell[keep]
  Q <- Q_cell[keep]

  n_NP   <- sum(m)
  N_hat  <- sum(m * w)
  mu_hat <- sum(w * S) / N_hat
  z      <- stats::qnorm(1 - (1 - level) / 2)

  ## ---- Simplified variance ----------------------------------------------
  ## EE3 alone, with the weights treated as fixed. Over records this is
  ## sum_i (w_i (y_i - mu))^2; within a cell w is constant, so that cell's
  ## contribution is w^2 (sum y^2 - 2 mu sum y + m mu^2). No model matrix is
  ## involved, which is why this needs nothing but the cell outcome totals.
  var_naive <- sum((w / N_hat)^2 * (Q - 2 * mu_hat * S + m * mu_hat^2))
  se_naive  <- sqrt(var_naive)
  ci_naive  <- mu_hat + c(-1, 1) * z * se_naive

  cover <- function(ci) {
    if (is.null(truth)) NA_real_ else
      as.numeric(truth >= ci[1] && truth <= ci[2])
  }

  if (!want_stacked) {
    return(c(mu_hat         = mu_hat,
             var_naive      = var_naive,
             se_naive       = se_naive,
             ci_naive_lower = ci_naive[1],
             ci_naive_upper = ci_naive[2],
             cover_naive    = cover(ci_naive),
             N_hat          = N_hat,
             n_NP           = n_NP))
  }

  ## ---- Stacked variance -------------------------------------------------
  ## Every block below is the record-level sandwich collapsed onto cells.
  ## Within a cell the model matrix row, the weight and the fitted propensity
  ## are constant, so a block linear in the outcome needs the cell sum of y, a
  ## block quadratic in it needs the cell sum of y^2, and a block free of the
  ## outcome needs only the record count. On the reference side the meat is
  ## quadratic in the design weights, which is why the reference cell table has
  ## to carry sum(d^2) as well as sum(d).
  cat_temp <- fit$formula_used
  cref     <- fit$cells_ref
  if (is.null(cref)) {
    stop("`fit` carries no reference cell table, so the stacked variance ",
         "cannot be formed. Refit with this version of rails().", call. = FALSE)
  }
  if (is.null(cref$weight_sq)) {
    stop("The reference cell table has no `weight_sq` column, the within-cell ",
         "sum of squared design weights. The stacked variance is quadratic in ",
         "those weights, so cell totals alone are not enough. Rebuild the ",
         "reference cells with rails_cells(), which records it.", call. = FALSE)
  }

  D  <- as.numeric(cref$weight)
  D2 <- as.numeric(cref$weight_sq)

  ## Built on the full cell tables and then subset, so that factor levels
  ## absent from the kept rows still contribute their columns.
  X <- stats::model.matrix(cat_temp, as.data.frame(cells))[keep, , drop = FALSE]
  Z <- stats::model.matrix(cat_temp, as.data.frame(cref))

  theta <- fit$theta
  if (ncol(X) != length(theta)) {
    stop("Model matrix has ", ncol(X), " columns but the fitted ",
         "coefficient vector has ", length(theta), ". The cell table and the ",
         "fitted model disagree, most likely because factor levels differ.",
         call. = FALSE)
  }
  if (nrow(Z) != length(D)) {
    stop("Reference model matrix has ", nrow(Z), " rows but ", length(D),
         " cells; reference cells with missing covariates were dropped by ",
         "model.matrix().", call. = FALSE)
  }

  pi_NP <- stats::plogis(as.numeric(X %*% theta))
  pi_P  <- stats::plogis(as.numeric(Z %*% theta))
  T_pop <- create_v3(fit$pop_totals, colnames(X))
  if (length(T_pop) != ncol(X)) {
    stop("Population margins are missing for ",
         paste(setdiff(colnames(X), names(T_pop)), collapse = ", "), ".",
         call. = FALSE)
  }

  q   <- ncol(X)
  n_P <- if (!is.null(cref$n)) sum(as.numeric(cref$n)) else NA_real_

  ## psi2 is constant within a cell; psi3 sums to t_cell over one and its
  ## squares sum to u_cell.
  A      <- (X * w - matrix(T_pop / n_NP, nrow = nrow(X), ncol = q,
                            byrow = TRUE)) / N_hat
  t_cell <- w * (S - m * mu_hat) / N_hat
  u_cell <- (w / N_hat)^2 * (Q - 2 * mu_hat * S + m * mu_hat^2)

  ## ---- Meat -------------------------------------------------------------
  M_11 <- crossprod(Z * (pi_P * sqrt(D2)) / N_hat) +
          crossprod(X * sqrt(m) / N_hat)
  M_12 <- crossprod(X * (m / N_hat), A)
  M_13 <- crossprod(X / N_hat, t_cell)
  M_22 <- crossprod(A * sqrt(m))
  M_23 <- crossprod(A, t_cell)
  M_33 <- sum(u_cell)

  ## ---- Bread ------------------------------------------------------------
  B_11 <-  crossprod(Z * sqrt(D * pi_P * (1 - pi_P))) / N_hat
  B_22 <- -crossprod(X * sqrt(m * w))                 / N_hat
  B_33 <- 1
  B_21 <- crossprod(X, X * (m * w * (1 - pi_NP)))     / N_hat

  B_31 <- matrix(colSums(X * (w * (1 - pi_NP) * (S - m * mu_hat))) / N_hat,
                 nrow = 1)
  B_32 <- matrix(-colSums(X * (w * (S - m * mu_hat))) / N_hat, nrow = 1)

  ## ---- Bread inverse ----------------------------------------------------
  B_11_inv <- solve(B_11)
  B_22_inv <- solve(B_22)

  C_33 <- 1 / B_33
  C_32 <- -C_33 * B_32 %*% B_22_inv
  C_31 <- -C_33 * (B_31 - B_32 %*% B_22_inv %*% B_21) %*% B_11_inv

  ## ---- Components -------------------------------------------------------
  term_33 <- C_33^2 * M_33
  term_11 <- as.numeric(C_31 %*% M_11 %*% t(C_31))
  term_22 <- as.numeric(C_32 %*% M_22 %*% t(C_32))
  term_12 <- as.numeric(2 * C_31 %*% M_12 %*% t(C_32))
  term_13 <- as.numeric(2 * C_33 * C_31 %*% M_13)
  term_23 <- as.numeric(2 * C_33 * C_32 %*% M_23)

  ## term_33 is the simplified variance, recomputed here from the same psi3.
  var_full <- term_33 + term_11 + term_22 + term_12 + term_13 + term_23

  if (var_full < 0) {
    warning("rails_var(): the stacked variance is negative; using its absolute ",
            "value. This usually signals a near-singular calibration model.",
            call. = FALSE)
    var_full <- abs(var_full)
  }

  se_full <- sqrt(var_full)
  ci_full <- mu_hat + c(-1, 1) * z * se_full

  cover_naive <- cover(ci_naive)
  cover_full  <- cover(ci_full)

  stacked_only <- c(
    mu_hat           = mu_hat,
    var_full         = var_full,
    se_full          = se_full,
    ci_full_lower    = ci_full[1],
    ci_full_upper    = ci_full[2],
    cover_full       = cover_full,
    term_33 = term_33,
    term_11 = term_11,
    term_22 = term_22,
    term_12 = term_12,
    term_13 = term_13,
    term_23 = term_23,
    N_hat   = N_hat,
    n_NP    = n_NP,
    n_P     = n_P
  )
  if (type == "stacked") return(stacked_only)

  c(mu_hat           = mu_hat,
    var_naive        = var_naive,
    se_naive         = se_naive,
    ci_naive_lower   = ci_naive[1],
    ci_naive_upper   = ci_naive[2],
    var_full         = var_full,
    se_full          = se_full,
    ci_full_lower    = ci_full[1],
    ci_full_upper    = ci_full[2],
    ratio_full_naive = se_full / se_naive,
    cover_naive      = cover_naive,
    cover_full       = cover_full,
    term_33 = term_33,
    term_11 = term_11,
    term_22 = term_22,
    term_12 = term_12,
    term_13 = term_13,
    term_23 = term_23,
    N_hat   = N_hat,
    n_NP    = n_NP,
    n_P     = n_P)
}

#' Per-stratum variance for a subgroup fit
#'
#' Runs the single-fit variance within each stratum. Each stratum has its own
#' propensity model, its own selected terms and its own margins, so the
#' estimating equations are genuinely separate; what the per-stratum sandwich
#' does not carry is the covariance between strata.
#'
#' @keywords internal
#' @noRd
rails_var_subgroup <- function(fit, y, type, truth, level) {

  ## `y` is per record for a microdata fit and per cell for an aggregated one.
  ## Either way `subgroup_rows` indexes the object the stratum was cut from, so
  ## the same slice works for both.
  cell_y <- is.data.frame(y) || is.matrix(y)
  if (cell_y) {
    y   <- as.data.frame(y)
    n_y <- nrow(y)
  } else {
    y   <- as.vector(y)
    n_y <- length(y)
  }

  if (is.null(fit$weights) && !cell_y) {
    stop("This subgroup fit was built from cell tables, so `y` must be per ",
         "cell too: a data frame with columns `sum` and `sumsq`, one row per ",
         "row of the cell table given to rails_subgroup().", call. = FALSE)
  }
  if (!is.null(fit$weights) && !cell_y && n_y != length(fit$weights)) {
    stop("`y` has length ", n_y, " but the non-probability sample has ",
         length(fit$weights), " rows.", call. = FALSE)
  }
  if (!length(fit$fits)) {
    stop("No stratum was fitted, so there is no variance to compute.",
         call. = FALSE)
  }

  levs <- names(fit$fits)

  ## `truth` may be one value for every stratum, or one per stratum named by
  ## level. Simulations normally want the latter.
  truth_for <- function(lev) {
    if (is.null(truth)) return(NULL)
    if (length(truth) == 1L && is.null(names(truth))) return(unname(truth))
    if (!lev %in% names(truth)) {
      stop("`truth` has no entry for ", fit$by, " = ", lev,
           ". Give one value for all strata, or a vector named by level.",
           call. = FALSE)
    }
    unname(truth[[lev]])
  }

  if (type %in% c("stacked", "both")) {
    warning("rails_var(): the stacked variance is computed within each stratum, ",
            "so it ignores the covariance between strata. Per-stratum intervals ",
            "are usable; anything pooling or differencing across strata is not.",
            call. = FALSE)
  }

  ## A stratum whose stepwise walk never converged has NA weights. Report NA for
  ## it and carry on, the way rails_subgroup() itself does, rather than losing
  ## the strata that did fit.
  usable <- vapply(levs, function(lev) {
    w <- fit$fits[[lev]]$weights
    !is.null(w) && !all(is.na(w))
  }, logical(1))

  if (!any(usable)) {
    stop("No stratum has usable weights, so there is no variance to compute.",
         call. = FALSE)
  }
  if (any(!usable)) {
    warning("rails_var(): no converged model for ", fit$by, " = ",
            paste(levs[!usable], collapse = ", "),
            "; those rows are NA.", call. = FALSE)
  }

  rows <- lapply(levs[usable], function(lev) {
    f <- fit$fits[[lev]]
    i <- f$subgroup_rows
    if (is.null(i)) {
      stop("This subgroup fit predates per-stratum row tracking; refit with ",
           "rails_subgroup().", call. = FALSE)
    }
    y_i <- if (cell_y) y[i, , drop = FALSE] else y[i]
    v <- rails_var_one(f, y_i, type, truth_for(lev), level)
    data.frame(subgroup = lev, as.list(v), check.names = FALSE,
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)

  if (any(!usable)) {
    blank <- out[1, , drop = FALSE]
    blank[, -1] <- NA_real_
    filler <- blank[rep(1, sum(!usable)), , drop = FALSE]
    filler$subgroup <- levs[!usable]
    out <- rbind(out, filler)
    out <- out[match(levs, out$subgroup), , drop = FALSE]
  }

  names(out)[1] <- fit$by
  rownames(out) <- NULL
  attr(out, "between_stratum_covariance") <- FALSE
  out
}
