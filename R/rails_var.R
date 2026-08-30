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
#'   [rails_subgroup()]. Either must have been fitted with `aggregated = FALSE`;
#'   the stacked variance additionally needs `keep_data = TRUE` (the default).
#' @param y Outcome vector, aligned to the rows of the non-probability sample
#'   passed to [rails()]. Binary or continuous.
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
#' @section Why individual records are required:
#' The estimator itself runs on aggregated cells, but the variance does not: it
#' needs one `y` per record. For the stacked form there is a further reason.
#' Within a covariate cell the model matrix row, the weight and the fitted
#' propensity are all constant, yet the outcome is not, so the meat matrix wants
#' within-cell sums of `y` rather than cell means. `rails_var()` therefore
#' rebuilds individual-level model matrices from the data retained by
#' `keep_data = TRUE`, and errors on a fit built from cell tables.
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
  if (isTRUE(fit$aggregated)) {
    stop("rails_var() needs individual records, but `fit` was built from cell ",
         "tables (aggregated = TRUE). Refit from microdata.", call. = FALSE)
  }
  if (want_stacked && is.null(fit$data)) {
    stop("The stacked variance needs the retained covariate columns, but `fit` ",
         "was built with keep_data = FALSE. Refit with keep_data = TRUE, or use ",
         'type = "simplified".', call. = FALSE)
  }
  if (want_stacked && (is.null(fit$formula_used) || is.null(fit$theta))) {
    stop("`fit` has no converged model: no stepwise step succeeded, so there ",
         "are no weights to compute a variance for.", call. = FALSE)
  }

  y    <- as.vector(y)
  w_NP <- as.vector(fit$weights)

  if (length(y) != length(w_NP)) {
    stop("`y` has length ", length(y), " but the non-probability sample has ",
         length(w_NP), " rows.", call. = FALSE)
  }

  ## Rows dropped at aggregation (missing covariates) carry NA weights.
  ok <- !is.na(w_NP) & !is.na(y)
  if (any(!ok)) {
    warning("rails_var(): dropping ", sum(!ok),
            " record(s) with a missing weight or outcome.", call. = FALSE)
  }

  y    <- y[ok]
  w_NP <- w_NP[ok]

  if (!length(w_NP) || all(is.na(w_NP))) {
    stop("`fit` has no usable weights, so there is no variance to compute.",
         call. = FALSE)
  }

  n_NP   <- length(w_NP)
  N_hat  <- sum(w_NP)
  mu_hat <- sum(w_NP * y) / N_hat
  z      <- stats::qnorm(1 - (1 - level) / 2)

  ## ---- Simplified variance ----------------------------------------------
  ## The weights are treated as fixed, so this is EE3 alone. No model matrix is
  ## involved, which is why it is available even without keep_data.
  psi3_NP   <- w_NP * (y - mu_hat) / N_hat
  var_naive <- sum(psi3_NP^2)
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
  cat_temp <- fit$formula_used
  dt_np    <- fit$data$np[ok, , drop = FALSE]
  dt_ref   <- fit$data$ref
  d_P      <- as.vector(fit$data$d_ref)

  X_ps_NP <- stats::model.matrix(cat_temp, dt_np)
  X_ps_P  <- stats::model.matrix(cat_temp, dt_ref)
  X_cal   <- X_ps_NP

  theta <- fit$theta
  if (ncol(X_ps_NP) != length(theta)) {
    stop("Model matrix has ", ncol(X_ps_NP), " columns but the fitted ",
         "coefficient vector has ", length(theta), ". The retained data and ",
         "the fitted model disagree, most likely because factor levels differ.",
         call. = FALSE)
  }
  if (nrow(X_ps_P) != length(d_P)) {
    stop("Reference model matrix has ", nrow(X_ps_P), " rows but ", length(d_P),
         " design weights; reference records with missing covariates were ",
         "dropped by model.matrix().", call. = FALSE)
  }

  pi_NP <- stats::plogis(as.numeric(X_ps_NP %*% theta))
  pi_P  <- stats::plogis(as.numeric(X_ps_P  %*% theta))
  T_pop <- create_v3(fit$pop_totals, colnames(X_cal))
  if (length(T_pop) != ncol(X_cal)) {
    stop("Population margins are missing for ",
         paste(setdiff(colnames(X_cal), names(T_pop)), collapse = ", "), ".",
         call. = FALSE)
  }

  n_P <- length(d_P)
  q   <- ncol(X_cal)

  ## ---- Individual psi vectors -------------------------------------------
  ## EE1 and EE2 are scaled by 1/N_hat for numerical stability; EE3 is scaled
  ## by 1/N_hat by definition. The sandwich is invariant to this scaling.
  ## psi3_NP is the same vector the simplified variance used above.
  psi1_P  <- -X_ps_P * d_P * pi_P / N_hat
  psi1_NP <- X_ps_NP / N_hat
  psi2_NP <- (X_cal * w_NP -
                matrix(T_pop / n_NP, nrow = n_NP, ncol = q, byrow = TRUE)) / N_hat

  ## ---- Meat -------------------------------------------------------------
  M_11 <- crossprod(psi1_P) + crossprod(psi1_NP)
  M_12 <- crossprod(psi1_NP, psi2_NP)
  M_13 <- crossprod(psi1_NP, psi3_NP)
  M_22 <- crossprod(psi2_NP)
  M_23 <- crossprod(psi2_NP, psi3_NP)
  M_33 <- sum(psi3_NP^2)

  ## ---- Bread ------------------------------------------------------------
  ## Rows 1 and 2 inherit the 1/N_hat scaling from EE1 and EE2; row 3 does not.
  B_11 <-  crossprod(X_ps_P * sqrt(d_P * pi_P * (1 - pi_P))) / N_hat
  B_22 <- -crossprod(X_cal  * sqrt(w_NP))                    / N_hat
  B_33 <- 1
  B_21 <- t(X_cal) %*% (X_ps_NP * w_NP * (1 - pi_NP))        / N_hat

  multiplier_31 <- w_NP * (1 - pi_NP) * (y - mu_hat) / N_hat
  B_31 <- matrix(colSums(sweep(X_ps_NP, 1, multiplier_31, "*")), nrow = 1)

  multiplier_32 <- w_NP * (y - mu_hat) / N_hat
  B_32 <- matrix(-colSums(sweep(X_cal, 1, multiplier_32, "*")), nrow = 1)

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

  if (is.null(fit$weights)) {
    stop("rails_var() needs individual records, but this subgroup fit was ",
         "built from cell tables (aggregated = TRUE). Refit from microdata.",
         call. = FALSE)
  }
  y <- as.vector(y)
  if (length(y) != length(fit$weights)) {
    stop("`y` has length ", length(y), " but the non-probability sample has ",
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
    v <- rails_var_one(f, y[i], type, truth_for(lev), level)
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
