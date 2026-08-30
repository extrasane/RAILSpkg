#' Greedy variable selection: test one candidate interaction
#'
#' Fits the nested propensity score model with one candidate interaction term
#' added and returns a likelihood-ratio test of that term against the current
#' model. [rails()] calls this once per candidate at every round of greedy
#' variable selection (GVS), the forward step of the RAILS algorithm; the
#' candidate with the largest significant increment per degree of freedom is the
#' one added.
#'
#' @param x Term or terms to add, e.g. `"agegroup:sex:income"`. A vector is
#'   tested as a single increment, which is how [rails()] keeps a selected
#'   higher-order term together with any lower-order components it needs.
#' @param terms_current Character vector of terms already in the current model.
#' @param loglik Log-likelihood of the current model.
#' @param cells_ref Aggregated reference-sample cell counts (column `weight`).
#' @param cells_np Aggregated non-probability cell counts (column `weight`).
#' @param m_ref Model matrix for the current model on `cells_ref`, used only for
#'   its column count, to compute degrees of freedom.
#' @param maxit Maximum Newton-Raphson iterations.
#' @param tol Convergence tolerance on the squared Newton step. Looser than
#'   [nps_weights()] because only the log-likelihood is needed, not the weights.
#'
#' @return A length-5 numeric vector `c(loglik, LRT, p-value, LRT/df, df)`, or
#'   `rep(NA, 5)` if any cell count is zero or the solver fails.
#'
#' @aliases gvs
#' @seealso [rails()], which runs the GVS loop, and [nps_weights()] for the
#'   model being tested.
#'
#' @export
gvs_step <- function(x, terms_current, loglik, cells_ref, cells_np, m_ref,
                     maxit = 1000, tol = 1e-3) {
  cat_formula <- stats::formula(
    paste0("~", paste0(paste0(terms_current, collapse = "+"), "+",
                       paste0(x, collapse = "+")))
  )
  t_ref <- Matrix::sparse.model.matrix(cat_formula, cells_ref)
  t_np  <- Matrix::sparse.model.matrix(cat_formula, cells_np)
  theta <- rep(0, ncol(t_np))
  pia   <- rep(1 / 2, nrow(t_ref))
  W_ref <- cells_ref$weight
  W_np  <- cells_np$weight
  U1    <- as.numeric(Matrix::crossprod(t_np, W_np))

  if (any(U1 == 0)) return(rep(NA, 5))

  res       <- 1
  temp_LKHD <- NA
  iter      <- 0

  tryCatch({
    while (res >= tol) {
      iter <- iter + 1
      if (iter > maxit) stop("non-convergence")
      W1     <- W_ref * pia
      Uscore <- U1 - as.numeric(Matrix::crossprod(t_ref, W1))
      Hsov   <- solve(as.matrix(Matrix::crossprod(t_ref, t_ref * (W1 * (1 - pia)))))
      theta0 <- as.numeric(theta + Hsov %*% Uscore)
      pia    <- stats::plogis(as.numeric(t_ref %*% theta0))
      res    <- sum((Hsov %*% Uscore)^2)
      theta  <- theta0
      temp_LKHD <- sum(U1 * theta) -
        sum(W_ref * log1p(exp(as.numeric(t_ref %*% theta))))
    }
    df   <- ncol(t_ref) - ncol(m_ref)
    LRT  <- 2 * (temp_LKHD - loglik)
    pval <- 1 - stats::pchisq(LRT, df = df)
    c(temp_LKHD, LRT, pval, LRT / df, df)
  }, error = function(e) rep(NA, 5))
}
