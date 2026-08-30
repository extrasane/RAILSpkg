#' Aggregate microdata into covariate cells
#'
#' Collapses individual-level records into one row per unique combination of
#' `vars`, with a `weight` column giving the cell total. RAILS runs entirely on
#' these cells, so both the non-probability sample and the reference sample must
#' be aggregated the same way before estimation. [rails()] calls this for you
#' when given microdata; call it directly when you want to inspect or cache the
#' cell tables.
#'
#' Factor levels are **not** dropped, so cells present in one sample but empty in
#' the other still line up in the model matrices.
#'
#' @param data A data frame of individual records.
#' @param vars Character vector of covariate names to aggregate over.
#' @param weights Optional. Either a numeric vector of length `nrow(data)` or
#'   the name of a column of `data` holding design weights. `NULL` (the default)
#'   counts each row as 1, which is what you want for an unweighted biobank
#'   sample; pass the survey weight for a probability reference sample.
#' @param drop_na Drop rows with a missing value in any of `vars`. Model
#'   matrices would silently drop them anyway, so the default warns and removes
#'   them up front.
#'
#' @return A data frame with one row per cell: the columns named in `vars`, plus
#'   `weight`. It carries class `rails_cells` and an attribute `row_cell`, the
#'   integer index mapping each row of `data` to its cell, which is what lets
#'   [rails()] return one weight per input record.
#'
#' @seealso [rails_totals()] for the matching population margins.
#'
#' @examples
#' pop <- rails_simulate(500, seed = 1)
#' rails_cells(pop[pop$aou == 1, ], c("agegroup", "sex"))
#'
#' # A reference sample carrying design weights:
#' ref <- pop[pop$s == 1, ]
#' ref$dweight <- 1 / ref$ps
#' rails_cells(ref, c("agegroup", "sex"), weights = "dweight")
#'
#' @export
rails_cells <- function(data, vars, weights = NULL, drop_na = TRUE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop("`data` is missing the variable(s): ",
         paste(missing_vars, collapse = ", "), ".", call. = FALSE)
  }

  w <- if (is.null(weights)) {
    rep(1, nrow(data))
  } else if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) {
      stop("`weights` names a column not found in `data`: ", weights, ".",
           call. = FALSE)
    }
    as.numeric(data[[weights]])
  } else {
    if (length(weights) != nrow(data)) {
      stop("`weights` must have length nrow(data) (", nrow(data), "), not ",
           length(weights), ".", call. = FALSE)
    }
    as.numeric(weights)
  }

  keep <- rep(TRUE, nrow(data))
  if (isTRUE(drop_na)) {
    keep <- stats::complete.cases(data[, vars, drop = FALSE]) & !is.na(w)
    if (any(!keep)) {
      warning("rails_cells(): dropped ", sum(!keep),
              " row(s) with missing covariate or weight values.", call. = FALSE)
    }
  }

  sub <- data[keep, vars, drop = FALSE]
  key <- do.call(paste, c(unname(lapply(sub, as.character)), sep = "\r"))
  first <- !duplicated(key)
  idx   <- match(key, key[first])

  cells <- sub[first, , drop = FALSE]
  cells$weight <- as.numeric(rowsum(w[keep], idx))
  rownames(cells) <- NULL

  ## row_cell is indexed against the *original* rows; dropped rows get NA.
  row_cell <- rep(NA_integer_, nrow(data))
  row_cell[keep] <- idx

  attr(cells, "row_cell") <- row_cell
  attr(cells, "vars")     <- vars
  class(cells) <- c("rails_cells", class(cells))
  cells
}

#' Population margins matching a cell table
#'
#' Builds the named vector of population totals that [rails()] rakes to: every
#' one-way, two-way, ... up to `order`-way margin of `vars`, named exactly as the
#' calibration model matrix columns are named. Assembling this by hand is the
#' most common source of errors, because the names have to match what
#' [survey::calibrate()] expects.
#'
#' @param cells A cell table from [rails_cells()], or any data frame with the
#'   columns in `vars` plus a `weight` column.
#' @param vars Character vector of covariate names. Defaults to the variables
#'   recorded on a `rails_cells` object.
#' @param order Highest interaction order to compute margins for. Must be at
#'   least as large as the highest order [rails()] will select, otherwise raking
#'   has no target for the selected terms.
#' @param weight_col Name of the cell-total column.
#'
#' @return A named numeric vector of margins, suitable as the `pop_totals`
#'   argument of [rails()].
#'
#' @seealso [rails_cells()].
#'
#' @examples
#' pop <- rails_simulate(500, seed = 1)
#' ref <- rails_cells(pop[pop$s == 1, ], c("agegroup", "sex"))
#' rails_totals(ref, order = 2)
#'
#' @export
rails_totals <- function(cells, vars = attr(cells, "vars"), order = 3,
                         weight_col = "weight") {
  if (is.null(vars)) {
    stop("`vars` could not be inferred; pass it explicitly.", call. = FALSE)
  }
  check_cells(cells, vars, "`cells`", weight_col)

  terms <- interaction_terms(vars, order)
  mat   <- Matrix::sparse.model.matrix(
    terms_formula(terms), data = as.data.frame(cells), keep.order = TRUE
  )
  stats::setNames(
    as.numeric(Matrix::crossprod(mat, cells[[weight_col]])),
    colnames(mat)
  )
}
