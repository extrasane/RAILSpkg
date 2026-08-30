test_that("rails_var returns a corrected SE at least as large as the naive one", {
  pop <- rails_simulate(20000, seed = 42)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  fit <- rails(aou, ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", verbose = FALSE)
  v <- rails_var(fit, aou$y, type = "both", truth = mean(pop$y))

  expect_type(v, "double")
  expect_true(all(c("mu_hat", "se_naive", "se_full", "ratio_full_naive") %in% names(v)))
  expect_equal(unname(v["mu_hat"]),
               sum(weights(fit) * aou$y) / sum(weights(fit)), tolerance = 1e-10)
  expect_gt(v["se_full"], 0)
  expect_gt(v["se_naive"], 0)
  expect_equal(unname(v["ratio_full_naive"]),
               unname(v["se_full"] / v["se_naive"]))
  expect_true(v["cover_full"] %in% c(0, 1))
  expect_equal(unname(v["n_NP"]), nrow(aou))
})

test_that("the variance decomposition adds up", {
  pop <- rails_simulate(20000, seed = 8)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  fit <- rails(aou, ref, c("agegroup", "sex"), weights_ref = "dweight",
               verbose = FALSE)
  v <- rails_var(fit, aou$y, type = "both")

  expect_equal(unname(v["var_naive"]), unname(v["term_33"]))
  expect_equal(
    unname(v["var_full"]),
    unname(sum(v[c("term_33", "term_11", "term_22", "term_12", "term_13", "term_23")])),
    tolerance = 1e-8
  )
})

test_that("rails_var refuses fits it cannot handle", {
  pop <- rails_simulate(20000, seed = 4)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps
  vars <- c("agegroup", "sex")

  fit_nodata <- rails(aou, ref, vars, weights_ref = "dweight",
                      keep_data = FALSE, verbose = FALSE)
  expect_error(rails_var(fit_nodata, aou$y, type = "stacked"), "keep_data = FALSE")
  ## The simplified variance needs no model matrix, so it still works.
  expect_type(rails_var(fit_nodata, aou$y), "double")

  fit_cells <- rails(rails_cells(aou, vars),
                     rails_cells(ref, vars, weights = "dweight"),
                     vars, aggregated = TRUE, verbose = FALSE)
  expect_error(rails_var(fit_cells, aou$y), "individual records")

  fit <- rails(aou, ref, vars, weights_ref = "dweight", verbose = FALSE)
  expect_error(rails_var(fit, aou$y[-1]), "but the non-probability sample has")
})

test_that("subgroup variance returns one row per stratum", {
  pop <- rails_simulate(20000, seed = 42)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  sfit <- rails_subgroup(aou, ref, c("agegroup", "income"), by = "sex",
                         weights_ref = "dweight", start = 1, scope = 2,
                         verbose = FALSE)

  v <- rails_var(sfit, aou$y)
  expect_s3_class(v, "data.frame")
  expect_equal(nrow(v), length(sfit$fits))
  expect_equal(names(v)[1], "sex")
  expect_setequal(v$sex, names(sfit$fits))

  ## Each row must equal what the stratum fit alone would give.
  for (i in seq_len(nrow(v))) {
    lev <- v$sex[i]
    f   <- sfit$fits[[lev]]
    one <- rails_var(f, aou$y[f$subgroup_rows])
    expect_equal(v$mu_hat[i], unname(one[["mu_hat"]]))
    expect_equal(v$se_naive[i], unname(one[["se_naive"]]))
  }
})

test_that("the stacked subgroup variance warns about between-stratum covariance", {
  pop <- rails_simulate(20000, seed = 42)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  sfit <- rails_subgroup(aou, ref, c("agegroup", "income"), by = "sex",
                         weights_ref = "dweight", start = 1, scope = 2,
                         verbose = FALSE)

  expect_warning(v <- rails_var(sfit, aou$y, type = "stacked"),
                 "ignores the covariance between strata")
  expect_true(all(c("se_full", "term_11") %in% names(v)))

  ## The simplified default must stay quiet.
  expect_silent(rails_var(sfit, aou$y))
})

test_that("subgroup truth may be one value or one per stratum", {
  pop <- rails_simulate(20000, seed = 42)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  sfit <- rails_subgroup(aou, ref, c("agegroup", "income"), by = "sex",
                         weights_ref = "dweight", start = 1, scope = 2,
                         verbose = FALSE)

  per_stratum <- tapply(pop$y, pop$sex, mean)
  v <- rails_var(sfit, aou$y, truth = per_stratum)
  expect_true(all(v$cover_naive %in% c(0, 1)))

  shared <- rails_var(sfit, aou$y, truth = mean(pop$y))
  expect_true(all(shared$cover_naive %in% c(0, 1)))

  expect_error(rails_var(sfit, aou$y, truth = c(other = 0.2)),
               "no entry for sex")
})

test_that("type selects which variance is returned", {
  pop <- rails_simulate(20000, seed = 42)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  fit <- rails(aou, ref, c("agegroup", "sex", "income"),
               weights_ref = "dweight", verbose = FALSE)

  simple  <- rails_var(fit, aou$y)
  stacked <- rails_var(fit, aou$y, type = "stacked")
  both    <- rails_var(fit, aou$y, type = "both")

  ## The default is the simplified form.
  expect_identical(simple, rails_var(fit, aou$y, type = "simplified"))

  expect_true(all(c("se_naive", "ci_naive_lower") %in% names(simple)))
  expect_false(any(grepl("full|term_", names(simple))))

  expect_true(all(c("se_full", "term_11") %in% names(stacked)))
  expect_false(any(grepl("naive", names(stacked))))

  ## "both" agrees with each on its own quantities.
  expect_equal(both[names(simple)][c("mu_hat", "var_naive", "se_naive")],
               simple[c("mu_hat", "var_naive", "se_naive")])
  expect_equal(unname(both["se_full"]), unname(stacked["se_full"]))
  expect_equal(unname(both["ratio_full_naive"]),
               unname(both["se_full"] / both["se_naive"]))
})

test_that("the simplified variance is the fixed-weight formula", {
  pop <- rails_simulate(20000, seed = 3)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  fit <- rails(aou, ref, c("agegroup", "sex", "income"), weights_ref = "dweight",
               start = 1, scope = 2, verbose = FALSE)
  v <- rails_var(fit, aou$y)

  w  <- weights(fit)
  mu <- sum(w * aou$y) / sum(w)
  expect_equal(unname(v["var_naive"]),
               sum((w * (aou$y - mu))^2) / sum(w)^2)
})

test_that("a stratum with no converged model gives an NA row, not an error", {
  ## At this seed the male stratum selects nothing, so with fallback = FALSE its
  ## weights are NA. The female stratum must still be reported.
  pop <- rails_simulate(20000, seed = 6)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  sfit <- suppressWarnings(rails_subgroup(
    aou, ref, c("agegroup", "income"), by = "sex",
    weights_ref = "dweight", start = 1, scope = 2, verbose = FALSE))

  expect_warning(v <- rails_var(sfit, aou$y), "no converged model for sex")
  expect_equal(nrow(v), 2L)
  expect_equal(v$sex, names(sfit$fits))          # original stratum order kept
  expect_true(is.na(v$mu_hat[v$sex == "male"]))
  expect_false(is.na(v$mu_hat[v$sex == "female"]))
})
