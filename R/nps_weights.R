#' Nested propensity score weights
#'
#' Estimates non-probability-sample weights by maximizing the pseudo-likelihood
#' of the **nested propensity score** (NPS) model via Newton-Raphson on aggregated
#' covariate cells. The returned weights are scaled to sum to `nsiz`.
#'
#' This is the NPS step of the RAILS algorithm: the pseudo-likelihood merges the
#' non-probability cohort with the probability sample and the resulting
#' inverse fitted propensities serve as base weights for the calibration stage.
#'
#' @param formula A one-sided model formula, e.g. `~ agegroup + sex +
#'   agegroup:sex`.
#' @param cells_np Aggregated non-probability cell counts, with a column
#'   `weight` giving the cell size.
#' @param cells_ref Aggregated reference-sample cell counts, with a column
#'   `weight`.
#' @param nsiz Target population total the weights are scaled to. Defaults to
#'   the total reference-sample weight, `sum(cells_ref$weight)`.
#' @param theta_init Optional warm-start coefficient vector (e.g. the previous
#'   step's solution, zero-padded for new columns). The pseudo-likelihood is
#'   globally concave, so the warm start only reduces the iteration count.
#' @param maxit Maximum Newton-Raphson iterations.
#' @param tol Convergence tolerance on the squared Newton step.
#'
#' @return A numeric vector of propensity weights (one per cell), scaled to sum
#'   to `nsiz`, with the converged coefficient vector attached as attribute
#'   `"theta"`.
#'
#' @aliases nps
#' @seealso [rails()], which calls this at the start and at every step of the
#'   stepwise raking walk.
#'
#' @examples
#' pop <- rails_simulate(2000, seed = 1)
#' np  <- rails_cells(pop[pop$aou == 1, ], c("agegroup", "sex"))
#' ref <- rails_cells(pop[pop$s == 1, ],   c("agegroup", "sex"))
#' w   <- nps_weights(~ agegroup + sex, np, ref)
#' round(w, 2)
#'
#' @export
nps_weights <- function(formula, cells_np, cells_ref,
                       nsiz = sum(cells_ref$weight), theta_init = NULL,
                       maxit = 1000, tol = 1e-10) {
  m_np  <- Matrix::sparse.model.matrix(formula, cells_np)
  m_ref <- Matrix::sparse.model.matrix(formula, cells_ref)
  theta <- if (is.null(theta_init)) rep(0, ncol(m_np)) else
    c(theta_init, rep(0, ncol(m_np) - length(theta_init)))
  pia   <- stats::plogis(as.numeric(m_ref %*% theta))
  W_ref <- cells_ref$weight
  W_np  <- cells_np$weight
  U1    <- as.numeric(Matrix::crossprod(m_np, W_np))

  res  <- 1
  iter <- 0
  while (res >= tol) {
    iter <- iter + 1
    if (iter > maxit) {
      stop("nps_weights(): no convergence within ", maxit, " iterations.",
           call. = FALSE)
    }
    W1     <- W_ref * pia
    Uscore <- U1 - as.numeric(Matrix::crossprod(m_ref, W1))
    Hsov   <- tryCatch(
      solve(as.matrix(Matrix::crossprod(m_ref, m_ref * (W1 * (1 - pia))))),
      error = function(e) {
        stop("nps_weights(): the propensity model is not identified from the ",
             "reference sample -- its information matrix is singular. The model ",
             "has ", ncol(m_ref), " parameters but the reference sample occupies ",
             "only ", nrow(m_ref), " covariate cells. Use a smaller `start`, drop a ",
             "variable, or collapse sparse categories.", call. = FALSE)
      }
    )
    theta0 <- as.numeric(theta + Hsov %*% Uscore)
    pia    <- stats::plogis(as.numeric(m_ref %*% theta0))
    res    <- sum((Hsov %*% Uscore)^2)
    theta  <- theta0
  }

  d <- 1 / stats::plogis(as.numeric(m_np %*% theta))
  d <- d / sum(d) * nsiz
  attr(d, "theta") <- theta
  d
}
