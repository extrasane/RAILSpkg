## Internal helpers -- not exported.

#' Standardize an interaction term name
#'
#' Sorts the colon-separated components of a term alphabetically so that,
#' e.g., `"sex:region"` and `"region:sex"` compare equal.
#'
#' @param x A single character string, possibly containing `":"`.
#' @return The term with components sorted, rejoined by `":"`.
#' @keywords internal
#' @noRd
standardize_string <- function(x) {
  paste(sort(unlist(strsplit(x, ":"))), collapse = ":")
}

#' Match calibration targets to design columns
#'
#' Subsets a named vector of population margins (`v1`) to the calibration
#' column names of a survey design (`v2`), matching by standardized term name
#' and returning the values in `v2` order (which is the order
#' [survey::calibrate()] expects).
#'
#' @param v1 Named numeric vector of population margins.
#' @param v2 Character vector of calibration column names.
#' @return Numeric vector named by `v2`, containing the matched `v1` values.
#' @keywords internal
#' @noRd
create_v3 <- function(v1, v2) {
  v1_std <- vapply(names(v1), standardize_string, character(1), USE.NAMES = FALSE)
  v2_std <- vapply(v2,        standardize_string, character(1), USE.NAMES = FALSE)
  idx  <- match(v2_std, v1_std)
  keep <- !is.na(idx)
  stats::setNames(as.numeric(v1[idx[keep]]), v2[keep])
}

#' All interaction terms of `vars` up to (or of exactly) a given order
#'
#' @param vars Character vector of main-effect variable names.
#' @param order Integer interaction order.
#' @param exactly If `TRUE`, return only terms of size `order`; otherwise all
#'   terms of size `1:order`.
#' @return Character vector of terms, colon-separated.
#' @keywords internal
#' @noRd
interaction_terms <- function(vars, order, exactly = FALSE) {
  if (order < 1) return(character(0))
  ## An order above the number of variables has no terms of that size at all;
  ## capping it would silently return lower-order terms instead.
  if (exactly && order > length(vars)) return(character(0))
  order <- min(order, length(vars))
  sizes <- if (exactly) order else seq_len(order)
  unlist(lapply(sizes, function(k) {
    utils::combn(vars, k, FUN = function(x) paste(x, collapse = ":"))
  }), use.names = FALSE)
}

#' Resolve a model specification into a character vector of terms
#'
#' `start` and `scope` may each be given as an interaction order (an integer),
#' a one-sided formula, or an explicit character vector of terms. This turns any
#' of those into terms, and validates that every component names a variable in
#' `vars`.
#'
#' @param spec The user's specification.
#' @param vars Character vector of main-effect variable names.
#' @param what Argument name, for error messages.
#' @return Character vector of terms.
#' @keywords internal
#' @noRd
resolve_model <- function(spec, vars, what) {
  if (inherits(spec, "formula")) {
    spec <- attr(stats::terms(spec, simplify = TRUE), "term.labels")
  }
  if (is.numeric(spec) && length(spec) == 1L) {
    if (spec < 1 || spec != round(spec)) {
      stop("`", what, "` given as an order must be a positive whole number, not ",
           spec, ".", call. = FALSE)
    }
    return(interaction_terms(vars, spec))
  }
  if (!is.character(spec) || length(spec) == 0L) {
    stop("`", what, "` must be an interaction order, a one-sided formula, or a ",
         "character vector of terms.", call. = FALSE)
  }
  parts <- unique(unlist(strsplit(spec, ":", fixed = TRUE)))
  unknown <- setdiff(parts, vars)
  if (length(unknown)) {
    stop("`", what, "` names variable(s) not in `vars`: ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  unique(spec)
}

#' Complete a candidate term to a hierarchically well-formed increment
#'
#' Returns the term preceded by any of its lower-order components that are not
#' already in the model, lowest order first. Adding `a:b:c` to a model lacking
#' `a:b` would give a model matrix whose columns are coded differently from the
#' full model's, so the calibration margins would not line up; completing the
#' increment keeps every model along the walk well formed.
#'
#' With the published settings (`start = 2`, `scope = 3`) every two-way
#' component is already present, so this returns the term unchanged.
#'
#' @param term A single candidate term.
#' @param present_std Standardized names of the terms already in the model.
#' @return Character vector of terms to add, in order.
#' @keywords internal
#' @noRd
complete_term <- function(term, present_std) {
  comp <- strsplit(term, ":", fixed = TRUE)[[1]]
  k <- length(comp)
  if (k <= 1) return(term)
  subs <- unlist(lapply(seq_len(k - 1), function(m) {
    utils::combn(comp, m, FUN = function(z) paste(z, collapse = ":"))
  }), use.names = FALSE)
  subs_std <- vapply(subs, standardize_string, character(1), USE.NAMES = FALSE)
  c(subs[!subs_std %in% present_std], term)
}

#' Highest interaction order appearing in a set of terms
#'
#' @keywords internal
#' @noRd
max_term_order <- function(terms) {
  if (!length(terms)) return(1L)
  max(vapply(strsplit(terms, ":", fixed = TRUE), length, integer(1)))
}

#' Build a one-sided formula from a character vector of terms
#'
#' @keywords internal
#' @noRd
terms_formula <- function(terms) {
  stats::formula(paste0("~", paste(terms, collapse = "+")))
}

#' Emit a timestamped progress message when `verbose` is TRUE
#'
#' @keywords internal
#' @noRd
say <- function(verbose, ...) {
  if (isTRUE(verbose)) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
}

#' Validate that a cell table looks like one
#'
#' @keywords internal
#' @noRd
check_cells <- function(x, vars, what, weight_col = "weight") {
  if (!is.data.frame(x)) {
    stop(what, " must be a data frame.", call. = FALSE)
  }
  if (!weight_col %in% names(x)) {
    stop(what, " has no `", weight_col, "` column. Aggregated input must carry ",
         "one row per covariate cell with a `", weight_col, "` column giving ",
         "the cell total. Use rails_cells() to build one.", call. = FALSE)
  }
  missing_vars <- setdiff(vars, names(x))
  if (length(missing_vars)) {
    stop(what, " is missing the variable(s): ",
         paste(missing_vars, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}
