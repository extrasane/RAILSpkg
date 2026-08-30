## Regression against the pre-package implementation.
##
## `reference-original.rds` holds inputs and the outputs produced by running the
## original sources -- RAILS_Edit/RAILS/R/fun_rails_threeway.R,
## fun_sub_rails_threeway.R, and Functions/supplementary_functions.R -- which are
## what produced the published results. These tests assert the package
## reproduces those numbers exactly, so a refactor that changes an answer fails
## here rather than silently in someone's analysis.
##
## Tolerance is 0 throughout: bit-for-bit, not merely close. Every call passes
## lifo = "ascending", the direction the application code used; the package
## default is now "descending", per the manuscripts.

ref <- readRDS(test_path("reference-original.rds"))

global_fit <- function() {
  suppressWarnings(rails(
    ref$cells_np, ref$cells_ref, ref$vars,
    pop_totals = ref$pop_totals, aggregated = TRUE,
    benchmarks = TRUE, lifo = "ascending", verbose = FALSE
  ))
}

test_that("rails() reproduces fun.rails.threeway() exactly", {
  fit <- global_fit()

  ## as.numeric() on the reference: the original leaked survey's names and `eta`
  ## attribute into its output, which the package now strips.
  expect_equal(weights(fit), as.numeric(ref$global$d_rails), tolerance = 0)
  expect_identical(paste(c(fit$terms_base, fit$terms_selected), collapse = " + "),
                   ref$global$selected_terms)
  expect_identical(paste(fit$terms_used, collapse = " + "),
                   ref$global$calibrated_terms)
})

test_that("the benchmark methods reproduce exactly too", {
  fit <- global_fit()
  for (nm in c("d_unweighted", "d_cal1", "d_cal2", "d_nps1", "d_nps2",
               "d_nps1_rake", "d_nps2_rake")) {
    expect_equal(fit$benchmarks[[nm]], as.numeric(ref$global[[nm]]),
                 tolerance = 0, info = nm)
  }
})

test_that("the deprecated shim returns the original wide data frame", {
  old <- suppressWarnings(suppressMessages(
    fun.rails.threeway(ref$cells_np, ref$cells_ref, ref$pop_totals,
                       names_univar = ref$vars)
  ))

  expect_equal(old$d_rails, as.numeric(ref$global$d_rails), tolerance = 0)
  expect_equal(old$d_cal2,  as.numeric(ref$global$d_cal2),  tolerance = 0)
  expect_identical(unique(old$selected_terms),   ref$global$selected_terms)
  expect_identical(unique(old$calibrated_terms), ref$global$calibrated_terms)
})

test_that("rails_subgroup() reproduces fun.sub.rails.threeway() exactly", {
  s  <- ref$subgroup
  by <- s$by

  fit <- suppressWarnings(rails_subgroup(
    s$cells_np, s$cells_ref, s$vars, by = by,
    aggregated = TRUE, lifo = "ascending", verbose = FALSE
  ))

  expect_named(fit$fits, as.character(sort(unique(s$cells[[by]]))))

  for (lev in names(fit$fits)) {
    got  <- fit$fits[[lev]]$weights
    rows <- s$cells[[by]] == lev
    key_ref <- do.call(paste, s$cells[rows, s$vars])
    key_got <- do.call(paste, fit$fits[[lev]]$cells[s$vars])
    want <- as.numeric(s$d_rails[rows][match(key_got, key_ref)])

    expect_equal(got, want, tolerance = 0, info = paste(by, "=", lev))
  }
})

test_that("a stratum the original could not fit is still not fitted", {
  ## sex = 0 selects no three-way term, so the published code returns NA
  ## there. fallback = FALSE must preserve that rather than quietly
  ## substituting the raked starting model.
  s   <- ref$subgroup
  fit <- suppressWarnings(rails_subgroup(
    s$cells_np, s$cells_ref, s$vars, by = s$by,
    aggregated = TRUE, lifo = "ascending", verbose = FALSE
  ))

  na_ref <- tapply(is.na(s$d_rails), s$cells[[s$by]], all)
  na_got <- vapply(fit$fits, function(f) all(is.na(f$weights)), logical(1))

  expect_equal(as.logical(na_got[names(na_ref)]), as.logical(na_ref))
  expect_true(any(na_ref))   # the degenerate case is genuinely exercised
  expect_false(all(na_ref))  # and so is the converged one
})

test_that("rails_var() reproduces fun.rails.var() exactly", {
  v   <- ref$variance
  fit <- suppressWarnings(rails(v$np, v$ref, ref$vars, weights_ref = v$d_ref,
                                lifo = "ascending", verbose = FALSE))
  expect_false(anyNA(weights(fit)))

  got <- rails_var(fit, v$np$y, type = "both", truth = v$truth)

  for (nm in c("mu_hat", "var_naive", "se_naive", "var_full", "se_full",
               "ratio_full_naive", "term_11", "term_22", "term_33",
               "term_12", "term_13", "term_23", "N_hat", "n_NP", "n_P",
               "ci_naive_lower", "ci_naive_upper",
               "ci_full_lower", "ci_full_upper")) {
    expect_equal(unname(got[[nm]]), unname(v$out[[nm]]), tolerance = 0, info = nm)
  }

  ## Coverage indicators are renamed but must agree.
  expect_equal(unname(got[["cover_naive"]]), unname(v$out[["cover_naive_temp"]]))
  expect_equal(unname(got[["cover_full"]]),  unname(v$out[["cover_full_temp"]]))
})

test_that("type = 'simplified' matches the original's naive variance", {
  v   <- ref$variance
  fit <- suppressWarnings(rails(v$np, v$ref, ref$vars, weights_ref = v$d_ref,
                                lifo = "ascending", verbose = FALSE))
  got <- rails_var(fit, v$np$y, truth = v$truth)

  expect_equal(unname(got[["var_naive"]]), unname(v$out[["var_naive"]]),
               tolerance = 0)
  expect_equal(unname(got[["ci_naive_lower"]]), unname(v$out[["ci_naive_lower"]]),
               tolerance = 0)
  expect_equal(unname(got[["cover_naive"]]), unname(v$out[["cover_naive_temp"]]))
})
