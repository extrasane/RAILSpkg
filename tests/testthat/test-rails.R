make_samples <- function(n = 20000, seed = 42) {
  pop <- rails_simulate(n, seed = seed)
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps
  list(pop = pop, aou = pop[pop$aou == 1, ], ref = ref)
}

test_that("rails() returns one weight per record, summing to the population", {
  d <- make_samples()
  fit <- rails(d$aou, d$ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", verbose = FALSE)

  expect_s3_class(fit, "rails_fit")
  expect_equal(length(weights(fit)), nrow(d$aou))
  expect_false(anyNA(weights(fit)))
  expect_equal(sum(weights(fit)), sum(d$ref$dweight), tolerance = 1e-6)
  expect_true(all(weights(fit) > 0))
})

test_that("raking reproduces the reference one-way margins", {
  d <- make_samples()
  fit <- rails(d$aou, d$ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", verbose = FALSE)
  w <- weights(fit)

  got  <- tapply(w, d$aou$agegroup, sum)
  want <- tapply(d$ref$dweight, d$ref$agegroup, sum)
  expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-4)
})

test_that("RAILS reduces selection bias relative to the unweighted sample", {
  d <- make_samples()
  fit <- rails(d$aou, d$ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", verbose = FALSE)
  w <- weights(fit)

  truth      <- mean(d$pop$y)
  unweighted <- mean(d$aou$y)
  weighted   <- sum(w * d$aou$y) / sum(w)

  expect_lt(abs(weighted - truth), abs(unweighted - truth))
})

test_that("scope equal to start leaves nothing to select", {
  d <- make_samples()
  fit <- rails(d$aou, d$ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", scope = 2, verbose = FALSE)
  expect_length(fit$terms_selected, 0L)
  expect_false(anyNA(weights(fit)))
})

test_that("aggregated input reproduces the microdata fit", {
  d <- make_samples()
  vars <- c("agegroup", "sex", "income")

  cells_np  <- rails_cells(d$aou, vars)
  cells_ref <- rails_cells(d$ref, vars, weights = "dweight")

  fit_micro <- rails(d$aou, d$ref, vars, weights_ref = "dweight", verbose = FALSE)
  fit_cells <- rails(cells_np, cells_ref, vars, aggregated = TRUE, verbose = FALSE)

  expect_equal(fit_micro$weights_unit, fit_cells$weights, tolerance = 1e-8)
  expect_equal(fit_micro$terms_selected, fit_cells$terms_selected)
})

test_that("benchmarks are opt-in", {
  d <- make_samples()
  fit <- rails(d$aou, d$ref, c("agegroup", "sex"),
               weights_ref = "dweight", verbose = FALSE)
  expect_null(fit$benchmarks)

  fit_b <- rails(d$aou, d$ref, c("agegroup", "sex"), weights_ref = "dweight",
                 benchmarks = TRUE, verbose = FALSE)
  expect_named(fit_b$benchmarks,
               c("d_unweighted", "d_cal1", "d_cal2", "d_nps1", "d_nps2",
                 "d_nps1_rake", "d_nps2_rake"))
})

test_that("print and summary methods run", {
  d <- make_samples(2000, seed = 5)
  fit <- rails(d$aou, d$ref, c("agegroup", "sex"),
               weights_ref = "dweight", verbose = FALSE)
  expect_output(print(fit), "RAILS fit")
  expect_output(print(summary(fit)), "Design effect")
  expect_s3_class(summary(fit), "summary.rails_fit")
  expect_gte(summary(fit)$deff, 1)
})

test_that("input validation fires before any fitting", {
  d <- make_samples(5000, seed = 2)
  expect_error(rails(d$aou, d$ref, "agegroup", weights_ref = "dweight"),
               "at least two variables")
  ## Omitting the reference design weights makes the "population" the reference
  ## sample itself, which is smaller than the biobank in some cells; the
  ## propensity model then has no solution. The warning must fire regardless of
  ## what happens downstream, hence the try().
  expect_warning(
    try(rails(d$aou, d$ref, c("agegroup", "sex", "income"), scope = 2,
              verbose = FALSE), silent = TRUE),
    "weights_ref` is NULL"
  )
})

test_that("a non-identified propensity model fails with a useful message", {
  d <- make_samples(5000, seed = 2)
  expect_error(
    suppressWarnings(
      rails(d$aou, d$ref, c("agegroup", "sex", "income"), scope = 2,
            verbose = FALSE)
    ),
    "not identified from the reference sample"
  )
})

test_that("start and scope accept orders, formulas and explicit terms", {
  d <- make_samples()
  vars <- c("agegroup", "sex", "income")

  by_order <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                    start = 2, scope = 3, verbose = FALSE)
  by_terms <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                    start = c("agegroup", "sex", "income",
                              "agegroup:sex", "agegroup:income", "sex:income"),
                    scope = c("agegroup", "sex", "income",
                              "agegroup:sex", "agegroup:income", "sex:income",
                              "agegroup:sex:income"),
                    verbose = FALSE)
  by_formula <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                      start = ~ (agegroup + sex + income)^2,
                      scope = ~ (agegroup + sex + income)^3,
                      verbose = FALSE)

  expect_equal(weights(by_order), weights(by_terms))
  expect_equal(weights(by_order), weights(by_formula))
})

test_that("term components may be written in any order", {
  d <- make_samples()
  vars <- c("agegroup", "sex", "income")
  flipped <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                   start = c("agegroup", "sex", "income",
                             "sex:agegroup", "income:agegroup", "income:sex"),
                   scope = 3, verbose = FALSE)
  canonical <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                     start = 2, scope = 3, verbose = FALSE)
  ## "sex:agegroup" in start must cancel "agegroup:sex" in scope, so only the
  ## three-way term is a candidate and the fit is identical.
  expect_equal(weights(flipped), weights(canonical))
  expect_equal(flipped$terms_selected, canonical$terms_selected)
})

test_that("start = 1 offers two-way and three-way terms alike", {
  d <- make_samples()
  fit <- rails(d$aou, d$ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", start = 1, scope = 3, verbose = FALSE)
  expect_length(fit$terms_base, 3L)     # main effects only
  expect_length(fit$terms_scope, 7L)    # 3 main + 3 two-way + 1 three-way
})

test_that("scope must contain start", {
  d <- make_samples(5000, seed = 1)
  expect_error(
    rails(d$aou, d$ref, c("agegroup", "sex", "income"), weights_ref = "dweight",
          start = 3, scope = 2, verbose = FALSE),
    "must contain every term in `start`"
  )
  expect_error(
    rails(d$aou, d$ref, c("agegroup", "sex"), weights_ref = "dweight",
          start = c("agegroup", "nosuchvar"), verbose = FALSE),
    "names variable"
  )
})

test_that("the default is NA rather than a fallback, matching the paper", {
  ## seed 3 at this size selects no three-way term, which is exactly the case
  ## where fun.rails.threeway() returned NA.
  d <- make_samples(4000, seed = 3)
  vars <- c("agegroup", "sex", "income")

  expect_warning(
    strict <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                    verbose = FALSE),
    "no term from `scope` was selected"
  )
  expect_true(all(is.na(weights(strict))))

  relaxed <- rails(d$aou, d$ref, vars, weights_ref = "dweight",
                   fallback = TRUE, verbose = FALSE)
  expect_false(anyNA(weights(relaxed)))
  expect_equal(relaxed$terms_used, relaxed$terms_base)
})

## ---- LIFO direction --------------------------------------------------------
## Chen et al. (2025) define LIFO as descending: rake the full selected model,
## and pop the last-added term until raking converges. The package defaults to
## that; "ascending" reproduces the application code, which climbed from the
## starting model and halted at the first failure.

lifo_fixture <- function() {
  v   <- c("agegroup", "sex", "income")
  pop <- rails_simulate(20000, seed = 42)
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps
  list(
    vars = v,
    np   = as.data.frame(rails_cells(pop[pop$aou == 1, ], v)),
    ref  = as.data.frame(rails_cells(ref, v, weights = "dweight")),
    n    = sum(1 / ref$ps)
  )
}

test_that("descending is the default direction", {
  f <- lifo_fixture()
  tot <- rails_totals(f$ref, vars = f$vars, order = 2)

  default <- suppressWarnings(rails(f$np, f$ref, f$vars, pop_totals = tot,
                                    aggregated = TRUE, start = 1, scope = 2,
                                    verbose = FALSE))
  desc <- suppressWarnings(rails(f$np, f$ref, f$vars, pop_totals = tot,
                                 aggregated = TRUE, start = 1, scope = 2,
                                 lifo = "descending", verbose = FALSE))
  expect_equal(weights(default), weights(desc))
})

test_that("both directions agree when the full model rakes cleanly", {
  f <- lifo_fixture()
  tot <- rails_totals(f$ref, vars = f$vars, order = 2)

  asc <- suppressWarnings(rails(f$np, f$ref, f$vars, pop_totals = tot,
                                aggregated = TRUE, start = 1, scope = 2,
                                lifo = "ascending", verbose = FALSE))
  desc <- suppressWarnings(rails(f$np, f$ref, f$vars, pop_totals = tot,
                                 aggregated = TRUE, start = 1, scope = 2,
                                 lifo = "descending", verbose = FALSE))

  ## Nothing fails here, so the whole selected set survives either way.
  expect_equal(asc$n_used, asc$n_base + length(asc$terms_selected))
  expect_equal(weights(asc), weights(desc))
})

test_that("LIFO pops the last term when its margin cannot be met", {
  f <- lifo_fixture()
  tot <- rails_totals(f$ref, vars = f$vars, order = 2)

  ## Selection order here is agegroup:sex then sex:income. Make a sex:income
  ## margin unreachable -- larger than the population total -- so raking the
  ## full model must fail while the model without it still converges.
  target <- grep("^sex.*:income", names(tot), value = TRUE)[1]
  bad <- tot
  bad[target] <- f$n * 5

  fit <- suppressWarnings(rails(f$np, f$ref, f$vars, pop_totals = bad,
                                aggregated = TRUE, start = 1, scope = 2,
                                lifo = "descending", verbose = FALSE))

  expect_length(fit$terms_selected, 2L)          # both terms were selected
  expect_equal(fit$n_used, fit$n_base + 1L)      # but LIFO popped one
  expect_equal(setdiff(fit$terms_used, fit$terms_base), "agegroup:sex")
  expect_false(fit$converged)
  expect_false(anyNA(weights(fit)))
})

test_that("LIFO warns when it has to drop back", {
  f <- lifo_fixture()
  tot <- rails_totals(f$ref, vars = f$vars, order = 2)
  target <- grep("^sex.*:income", names(tot), value = TRUE)[1]
  bad <- tot
  bad[target] <- f$n * 5

  w <- capture_warnings(
    rails(f$np, f$ref, f$vars, pop_totals = bad, aggregated = TRUE,
          start = 1, scope = 2, lifo = "descending", verbose = FALSE))

  expect_true(any(grepl("LIFO dropped back to 1", w)))
  ## survey's own "Failed to converge" chatter is muffled; its text is folded
  ## into our warning instead of firing once per attempt.
  expect_false(any(grepl("^Failed to converge", w)))
  expect_true(any(grepl("Failed to converge", w)))
})

test_that("an unknown lifo direction is rejected", {
  f <- lifo_fixture()
  expect_error(
    rails(f$np, f$ref, f$vars, aggregated = TRUE, lifo = "sideways"),
    "'arg' should be one of"
  )
})

test_that("survey signals non-convergence as an error, not just a warning", {
  ## rails() treats a raking failure as "this model did not converge" by
  ## catching an error. That is only safe because survey::calibrate() stops
  ## rather than returning unconverged weights: grake() warns and sets
  ## attr(g, "failed"), and calibrate.survey.design2() turns that into
  ## stop("Calibration failed") whenever force = FALSE and no trim= was given.
  ##
  ## If a future survey release downgraded that to a warning -- or if someone
  ## added trim= or force=TRUE to the call in rails_engine() -- non-convergent
  ## weights would be silently accepted, and the muffling would hide the only
  ## remaining signal. This test fails first if that ever happens.
  d <- data.frame(x = factor(c("a", "a", "b", "b")), count = c(10, 10, 10, 10))
  des <- survey::svydesign(id = ~1, weights = ~count, data = d)
  impossible <- c("(Intercept)" = 40, "xb" = 400)   # a margin above the total

  expect_error(
    suppressWarnings(survey::calibrate(des, ~x, impossible, calfun = "raking",
                                       maxit = 50, epsilon = rep(1e-7, 2))),
    "Calibration failed"
  )

  ## And the warning that precedes it really is only a warning.
  expect_warning(
    try(survey::calibrate(des, ~x, impossible, calfun = "raking",
                          maxit = 50, epsilon = rep(1e-7, 2)), silent = TRUE),
    "Failed to converge"
  )
})
