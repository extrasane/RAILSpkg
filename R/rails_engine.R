## The RAILS algorithm, unwrapped. Internal.
##
## rails() wraps the return value in a `rails_fit` object; the deprecated
## fun.rails.threeway() reshapes it into the old wide data frame. Neither
## duplicates the algorithm.

#' @keywords internal
#' @noRd
rails_engine <- function(cells_np, cells_ref, pop_totals, vars,
                         start = 2, scope = 3, alpha = 0.05,
                         nsiz = sum(cells_ref$weight),
                         benchmarks = FALSE, fallback = FALSE,
                         lifo = c("descending", "ascending"), verbose = TRUE) {

  lifo <- match.arg(lifo)

  ## ---- Starting model and its likelihood --------------------------------
  terms_base  <- resolve_model(start, vars, "start")
  terms_scope <- resolve_model(scope, vars, "scope")

  ## Terms compare equal regardless of the order their components are written
  ## in, so "sex:region" in `start` covers "region:sex" in `scope`.
  std <- function(x) vapply(x, standardize_string, character(1), USE.NAMES = FALSE)
  missing_from_scope <- setdiff(std(terms_base), std(terms_scope))
  if (length(missing_from_scope)) {
    stop("`scope` must contain every term in `start`; it is the largest model ",
         "selection may reach. Missing: ",
         paste(missing_from_scope, collapse = ", "), ".", call. = FALSE)
  }

  cat_base <- terms_formula(terms_base)

  d_base <- nps_weights(cat_base, cells_np, cells_ref, nsiz = nsiz)
  m_ref  <- Matrix::sparse.model.matrix(cat_base, cells_ref)
  m_np   <- Matrix::sparse.model.matrix(cat_base, cells_np)
  theta_base <- attr(d_base, "theta")
  loglik_base <- sum(as.numeric(Matrix::crossprod(m_np, cells_np$weight)) * theta_base) -
    sum(cells_ref$weight * log1p(exp(as.numeric(m_ref %*% theta_base))))

  ## ---- Forward selection over everything in scope but not in start ------
  candidates   <- terms_scope[!std(terms_scope) %in% std(terms_base)]
  selectable   <- length(candidates) > 0
  terms_model  <- terms_base
  loglik_cur   <- loglik_base
  m_ref_cur    <- m_ref

  keep_going <- length(candidates) > 0
  while (keep_going) {
    if (length(candidates) == 0) break

    ## Each candidate is tested as a hierarchically complete increment: the term
    ## plus any lower-order components still missing. Under the published
    ## settings nothing is ever missing, so the increment is the term alone.
    increments <- lapply(candidates, complete_term, present_std = std(terms_model))

    temp_p <- vapply(increments, gvs_step,
                     FUN.VALUE = numeric(5),
                     terms_current = terms_model, loglik = loglik_cur,
                     cells_ref = cells_ref, cells_np = cells_np, m_ref = m_ref_cur)
    colnames(temp_p) <- candidates
    temp_sign <- temp_p[, !is.na(temp_p[3, ]) & temp_p[3, ] < alpha, drop = FALSE]
    if (ncol(temp_sign) == 0) {
      keep_going <- FALSE
    } else {
      selected_var <- colnames(temp_sign)[which.max(temp_sign[4, ])]
      added <- increments[[match(selected_var, candidates)]]
      say(verbose, "selection: added '", paste(added, collapse = "' + '"), "'")
      terms_model <- c(terms_model, added)
      loglik_cur  <- temp_sign[1, selected_var]
      m_ref_cur   <- Matrix::sparse.model.matrix(terms_formula(terms_model), cells_ref)
      ## Drop the selected term and anything its increment pulled in with it.
      candidates  <- candidates[!std(candidates) %in% std(terms_model)]
    }
  }

  terms_selected <- setdiff(terms_model, terms_base)
  cat_full <- terms_formula(terms_model)

  ## Subset pop_totals ONCE to the full selected model's margins; each stepwise
  ## step then subsets this small vector.
  pop_totals <- create_v3(
    pop_totals, colnames(Matrix::sparse.model.matrix(cat_full, cells_np))
  )

  ## ---- LIFO stepwise raking (direction set by `lifo`) -------------------
  n_base <- length(terms_base)
  n_max  <- length(terms_model)

  ## Rake a cell design with the given starting cell totals to the margins of
  ## `cat_f`. Returns NULL on non-convergence.
  rake_cells <- function(cat_f, start_cell) {
    temp_c <- as.data.frame(cells_np)
    temp_c$count <- start_cell
    clus_c <- survey::svydesign(id = ~1, weights = ~count, data = temp_c)
    vars_c <- create_v3(pop_totals, survey::cal_names(cat_f, clus_c))
    ## survey::grake() warns "Failed to converge" and sets attr(g, "failed") --
    ## it does not stop. calibrate.survey.design2() is what stops, with
    ## "Calibration failed", and only when `force = FALSE` and the flag survives;
    ## a `trim =` argument would clear the flag first. We pass neither `force`
    ## nor `trim`, so non-convergence always reaches the error branch here and
    ## unconverged weights cannot be mistaken for converged ones. A test pins
    ## that, so a future change to survey -- or to this call -- fails loudly.
    ##
    ## The warning is muffled because it would fire once per LIFO attempt; its
    ## text is folded into the error message instead, so the two always describe
    ## the same failed call. Any other warning passes through untouched.
    detail <- NA_character_
    withCallingHandlers(
      tryCatch(
        ## as.numeric() drops the names and the `eta` attribute survey attaches;
        ## they are its internals, not part of a weight vector.
        as.numeric(stats::weights(survey::calibrate(
          clus_c, cat_f, vars_c,
          calfun = "raking", maxit = 1e3,
          epsilon = rep(1e-7, length(vars_c))
        ))),
        error = function(e) {
          ## Pair the error with the muffled warning from the same call.
          rake_error <<- paste(c(conditionMessage(e), detail), collapse = "; ")
          NULL
        }
      ),
      warning = function(w) {
        if (grepl("Failed to converge", conditionMessage(w), fixed = TRUE)) {
          detail <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  rake_error <- NA_character_
  rake_reason <- function() {
    if (is.na(rake_error)) "" else paste0(" (", rake_error, ")")
  }

  scale_cells <- function(d) {
    ## as.numeric() drops the "theta" attribute nps_weights() carries; it belongs
    ## to the solver's return value, not to a weight vector.
    ct <- as.numeric(d) * cells_np$weight
    ct / sum(ct) * nsiz
  }

  weights_cell <- rep(NA_real_, nrow(cells_np))
  n_used       <- 0L
  theta_used   <- NULL
  terms_used   <- NA_character_
  theta_prev   <- theta_base

  ## Step 0: the base model itself, propensity-weighted then raked. This is the
  ## fallback when nothing is selected or the first higher-order step fails.
  ## When nothing is selectable, the starting model IS the requested model, not
  ## a fallback, so it is always kept. Otherwise `fallback` decides whether it
  ## stands in for a walk that never got off the ground.
  keep_base <- isTRUE(fallback) || !selectable
  base_cell <- rake_cells(cat_base, scale_cells(d_base))
  if (keep_base) {
    if (!is.null(base_cell)) {
      weights_cell <- base_cell
      n_used       <- n_base
      theta_used   <- theta_base
      terms_used   <- terms_base
    } else {
      warning("rails(): the starting model failed to converge under raking",
              rake_reason(), ".", call. = FALSE)
    }
  }

  if (selectable && length(terms_selected) == 0 && !isTRUE(fallback)) {
    warning("rails(): no term from `scope` was selected -- weights are NA. ",
            "Set fallback = TRUE to use the raked starting model instead.",
            call. = FALSE)
  }

  ## Fit the model on the first `k` terms and rake it. NULL if raking fails.
  try_model <- function(k, theta_init = NULL) {
    cat_temp <- terms_formula(terms_model[seq_len(k)])
    temp_d   <- nps_weights(cat_temp, cells_np, cells_ref, nsiz = nsiz,
                           theta_init = theta_init)
    step <- rake_cells(cat_temp, scale_cells(temp_d))
    if (is.null(step)) NULL else
      list(weights = step, theta = attr(temp_d, "theta"),
           terms = terms_model[seq_len(k)], n = k)
  }

  take <- function(fit) {
    weights_cell <<- fit$weights
    n_used       <<- fit$n
    theta_used   <<- fit$theta
    terms_used   <<- fit$terms
  }

  if (identical(lifo, "descending")) {
    ## The LIFO strategy of Chen et al. (2025): begin at the full selected model
    ## and pop the last-added term until raking converges.
    for (k in rev(seq.int(n_base + 1L, length.out = max(0L, n_max - n_base)))) {
      say(verbose, "LIFO ", n_max - k + 1L, "/", n_max - n_base,
          ": trying ", k - n_base, " selected term(s)")
      fit_k <- try_model(k)
      if (!is.null(fit_k)) {
        take(fit_k)
        if (k < n_max) {
          warning("rails(): raking did not converge with all ", n_max - n_base,
                  " selected term(s)", rake_reason(), "; LIFO dropped back to ",
                  k - n_base, ".", call. = FALSE)
        }
        break
      }
    }
    if (n_used <= n_base && n_max > n_base) {
      warning("rails(): raking did not converge for any selected term",
              rake_reason(), ".", call. = FALSE)
    }
  } else {
    ## Ascending: climb from the starting model and halt at the first failure.
    ## This is what the application code did; kept so its results reproduce.
    n_uni <- n_base
    while (n_uni < n_max) {
      n_uni <- n_uni + 1L
      say(verbose, "stepwise ", n_uni - n_base, "/", n_max - n_base,
          ": adding '", terms_model[n_uni], "'")
      fit_k <- try_model(n_uni, theta_init = theta_prev)
      if (is.null(fit_k)) {
        warning("rails(): covariate combination '", terms_model[n_uni],
                "' fails to pile up to higher order", rake_reason(),
                " -- keeping the previous model.", call. = FALSE)
        break
      }
      theta_prev <- fit_k$theta
      take(fit_k)
    }
  }

  ## ---- Benchmarks (all on the same aggregated cells) --------------------
  bench <- NULL
  if (isTRUE(benchmarks)) {
    cat_oneway  <- terms_formula(vars)
    d_nps1_cell <- scale_cells(nps_weights(cat_oneway, cells_np, cells_ref, nsiz = nsiz))
    d_nps2_cell <- scale_cells(d_base)
    d_eq_cell   <- cells_np$weight / sum(cells_np$weight) * nsiz

    na_cell <- rep(NA_real_, nrow(cells_np))
    or_na   <- function(x, what) {
      if (is.null(x)) {
        warning("rails(): benchmark '", what, "' failed to converge",
                rake_reason(), ".", call. = FALSE)
        na_cell
      } else x
    }

    ## Benchmarks are reported as per-record weights, like `weights()`. The
    ## unweighted one is the scalar nsiz / n computed directly rather than as a
    ## cell total divided back down, so it matches the published code bit for
    ## bit instead of to within rounding.
    per_unit <- function(x) x / cells_np$weight

    bench <- list(
      d_unweighted = rep(nsiz / sum(cells_np$weight), nrow(cells_np)),
      d_cal1       = per_unit(or_na(rake_cells(cat_oneway, d_eq_cell),   "d_cal1")),
      d_cal2       = per_unit(or_na(rake_cells(cat_base,   d_eq_cell),   "d_cal2")),
      d_nps1       = per_unit(d_nps1_cell),
      d_nps2       = per_unit(d_nps2_cell),
      d_nps1_rake  = per_unit(or_na(rake_cells(cat_oneway, d_nps1_cell), "d_nps1_rake")),
      d_nps2_rake  = per_unit(or_na(rake_cells(cat_base,   d_nps2_cell), "d_nps2_rake"))
    )
  }

  list(
    cells          = cells_np,
    weights_cell   = weights_cell,
    terms_base     = terms_base,
    terms_selected = terms_selected,
    terms_used     = terms_used,
    formula_used   = if (n_used > 0) terms_formula(terms_used) else NULL,
    theta          = theta_used,
    pop_totals     = pop_totals,
    nsiz           = nsiz,
    vars           = vars,
    terms_scope    = terms_scope,
    start          = start,
    scope          = scope,
    alpha          = alpha,
    n_used         = n_used,
    n_base         = n_base,
    selectable     = selectable,
    converged      = n_used == n_max && n_used > 0,
    benchmarks     = bench
  )
}
