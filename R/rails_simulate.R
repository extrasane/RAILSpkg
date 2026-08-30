#' Simulate a population with a biased volunteer sample
#'
#' Generates the synthetic population used by the RAILS simulation studies: five
#' categorical covariates, a binary outcome, a probability sample, and a
#' non-probability ("volunteer") sample whose inclusion depends on the same
#' covariates. Because the true population is known, the bias of any weighting
#' method can be read off directly. This is what the examples and vignettes run
#' on -- real biobank microdata cannot be redistributed.
#'
#' @param n Population size.
#' @param alpha Length-12 coefficient vector for the outcome model.
#' @param beta Length-12 coefficient vector for probability-sample inclusion.
#' @param gamma Length-12 coefficient vector for volunteer-sample inclusion.
#'   Making `gamma` differ from `beta` is what creates the selection bias RAILS
#'   is meant to remove.
#' @param seed Optional integer seed. `NULL` leaves the RNG state alone.
#'
#' @return A data frame with one row per population member:
#'   \describe{
#'     \item{`id`}{Row identifier.}
#'     \item{`agegroup`}{`young`, `middle`, `older`.}
#'     \item{`sex`}{`male`, `female`.}
#'     \item{`income`}{`p0-10`, `p10-50`, `p50-90`, `p90-100` -- cut at those
#'       percentiles of a truncated-Pareto draw.}
#'     \item{`homeown`}{`rent`, `own`.}
#'     \item{`edu`}{`<HS`, `HS`, `Some college`, `College+`.}
#'     \item{`y`,`py`}{Binary outcome and its true probability.}
#'     \item{`s`,`ps`,`fpc`}{Probability-sample indicator, its inclusion
#'       probability, and the population size.}
#'     \item{`aou`,`paou`}{Volunteer-sample indicator and its inclusion
#'       probability.}
#'   }
#'   Design weights for the probability sample are `1 / ps`.
#'
#' @section These are synthetic, not calibrated:
#' The variable names and level labels are there so worked examples read like the
#' application rather than like `agegroup:sex`. The marginal distributions behind
#' them are arbitrary and resemble no real population: `edu` in particular puts
#' only 0.5 percent in its top category. Nothing here should be read as a claim
#' about any real group.
#'
#' @section Difference from the paper code:
#' Income is drawn from a truncated Pareto by inverse-CDF sampling rather than
#' via `distributionsrd::rtruncpareto()`. Same distribution, one less
#' dependency; streams drawn under a given seed therefore differ from the
#' simulation archive.
#'
#' @examples
#' pop <- rails_simulate(1000, seed = 1)
#' head(pop[, c("agegroup", "sex", "income", "homeown", "edu", "y")])
#' table(pop$aou, pop$s)
#' mean(pop$y)                        # population prevalence
#' mean(pop$y[pop$aou == 1])          # what the volunteer sample would report
#'
#' @export
rails_simulate <- function(
    n,
    alpha = c(-3.5, 0.05, 0.50, 0.20, 0.40, 0.60, 0.010, 0.10, 0.30, 0.10, 0.20, 0.30),
    beta  = c(-2.2, 0.03, 0.20, 0.10, 0.20, 0.30, 0.005, 0.05, 0.10, 0.05, 0.10, 0.15),
    gamma = c(-4.8, 0.08, 0.60, 0.30, 0.60, 0.90, 0.020, 0.20, 0.40, 0.15, 0.30, 0.45),
    seed  = NULL) {

  if (!is.null(seed)) set.seed(seed)
  for (nm in c("alpha", "beta", "gamma")) {
    v <- get(nm)
    if (length(v) != 12) {
      stop("`", nm, "` must have length 12, not ", length(v), ".", call. = FALSE)
    }
  }

  ## Age
  x1    <- stats::rnorm(n, mean = 20, sd = 5)
  agegroup <- factor((x1 >= stats::quantile(x1, 0.30)) +
                     (x1 >= stats::quantile(x1, 0.75)),
                     labels = c("young", "middle", "older"))

  ## Sex (1 = female)
  x2  <- stats::rbinom(n, size = 1, prob = 0.65)
  sex <- factor(x2, levels = c(0, 1), labels = c("male", "female"))

  ## Income, truncated Pareto on [8, 300] with shape 1, by inverse CDF
  x3    <- rtrunc_pareto(n, lower = 8, upper = 300, shape = 1)
  qx3   <- stats::quantile(x3, probs = c(0.1, 0.5, 0.9))
  x3    <- (x3 >= qx3[1] & x3 < qx3[2]) * 1 +
           (x3 >= qx3[2] & x3 < qx3[3]) * 2 +
           (x3 >= qx3[3]) * 3
  income <- factor(x3, levels = 0:3,
                   labels = c("p0-10", "p10-50", "p50-90", "p90-100"))

  ## Two further covariates
  x5  <- sample(c(0, 1, 2, 3), n, replace = TRUE, prob = c(0.29, 0.30, 0.45, 0.005))
  edu <- factor(x5, levels = 0:3,
                labels = c("<HS", "HS", "Some college", "College+"))
  x4      <- stats::rbinom(n, size = 1, prob = 0.1 * x5)
  homeown <- factor(x4, levels = c(0, 1), labels = c("rent", "own"))

  ## Shared linear predictor, parameterized by a length-12 coefficient vector
  lp <- function(k) {
    k[1] + k[2] * x1 + k[3] * x2 + k[4] * (x3 == 1) + k[5] * (x3 == 2) +
      k[6] * (x3 == 3) + k[7] * x1 * x2 + k[8] * x2 * x3 + k[9] * x4 +
      k[10] * (x5 == 1) + k[11] * (x5 == 2) + k[12] * (x5 == 3)
  }

  py   <- stats::plogis(lp(alpha))
  ps   <- stats::plogis(lp(beta))
  paou <- stats::plogis(lp(gamma))

  data.frame(
    id       = seq_len(n),
    agegroup = agegroup, sex = sex, income = income,
    homeown  = homeown,  edu = edu,
    y     = stats::rbinom(n, size = 1, prob = py),
    py    = py,
    s     = stats::rbinom(n, size = 1, prob = ps),
    ps    = ps,
    fpc   = n,
    aou   = stats::rbinom(n, size = 1, prob = paou),
    paou  = paou,
    row.names = NULL
  )
}

#' Truncated Pareto by inverse CDF
#'
#' @keywords internal
#' @noRd
rtrunc_pareto <- function(n, lower, upper, shape) {
  u <- stats::runif(n)
  ratio <- (lower / upper)^shape
  lower / (1 - u * (1 - ratio))^(1 / shape)
}
