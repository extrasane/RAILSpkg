test_that("rails_subgroup fits every level and pools the weights", {
  pop <- rails_simulate(20000, seed = 7)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  fit <- rails_subgroup(aou, ref, c("agegroup", "income"), by = "sex",
                        weights_ref = "dweight", verbose = FALSE)

  expect_s3_class(fit, "rails_subgroup_fit")
  expect_equal(length(fit$fits), length(unique(aou$sex)))
  expect_equal(length(weights(fit)), nrow(aou))
  expect_false(anyNA(weights(fit)))

  ## Each stratum is scaled to its own reference total.
  got  <- tapply(weights(fit), aou$sex, sum)
  want <- tapply(ref$dweight, ref$sex, sum)
  expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-4)
})

test_that("the stratifier cannot also be a covariate", {
  pop <- rails_simulate(500, seed = 9)
  expect_error(
    rails_subgroup(pop, pop, c("agegroup", "sex"), by = "sex"),
    "must not also appear"
  )
})

test_that("rails_design hands weights to survey", {
  pop <- rails_simulate(3000, seed = 12)
  aou <- pop[pop$aou == 1, ]
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  fit <- rails(aou, ref, c("agegroup", "sex"), weights_ref = "dweight",
               verbose = FALSE)
  des <- rails_design(fit, aou)

  expect_s3_class(des, "survey.design")
  est <- survey::svymean(~y, des)
  expect_equal(as.numeric(est),
               sum(weights(fit) * aou$y) / sum(weights(fit)), tolerance = 1e-8)

  expect_error(rails_design(fit, aou[-1, ]), "rows but the fit carries")
})

test_that("deprecated shims still work and warn", {
  pop <- rails_simulate(20000, seed = 42)
  vars <- c("agegroup", "sex", "income")
  ref <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  cells_np  <- as.data.frame(rails_cells(pop[pop$aou == 1, ], vars))
  cells_ref <- as.data.frame(rails_cells(ref, vars, weights = "dweight"))
  tot <- rails_totals(cells_ref, vars = vars, order = 3)

  expect_warning(
    out <- suppressMessages(
      fun.rails.threeway(cells_np, cells_ref, tot, names_univar = vars)
    ),
    "deprecated"
  )
  expect_true(all(c("d_rails", "d_cal1", "d_cal2", "d_nps1", "d_nps2",
                    "d_nps1_rake", "d_nps2_rake", "d_unweighted",
                    "selected_terms", "calibrated_terms") %in% names(out)))
  expect_equal(nrow(out), nrow(cells_np))

  expect_warning(fun.out(c(1, 2, 3)), "deprecated")
})

test_that("rails() and fun.rails.threeway() give identical weights", {
  pop  <- rails_simulate(20000, seed = 42)
  vars <- c("agegroup", "sex", "income")
  aou  <- pop[pop$aou == 1, ]
  ref  <- pop[pop$s == 1, ]
  ref$dweight <- 1 / ref$ps

  cells_np  <- as.data.frame(rails_cells(aou, vars))
  cells_ref <- as.data.frame(rails_cells(ref, vars, weights = "dweight"))
  tot <- rails_totals(cells_ref, vars = vars, order = 3)

  new <- rails(cells_np, cells_ref, vars, pop_totals = tot,
               aggregated = TRUE, verbose = FALSE)
  old <- suppressWarnings(suppressMessages(
    fun.rails.threeway(cells_np, cells_ref, tot, names_univar = vars)
  ))

  expect_equal(weights(new), old$d_rails)
  expect_equal(paste(new$terms_used, collapse = " + "),
               unique(old$calibrated_terms))
})
