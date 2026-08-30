test_that("rails_cells collapses to unique covariate patterns", {
  d <- data.frame(
    a = factor(c("x", "x", "y", "y", "y")),
    b = factor(c(1, 1, 1, 2, 2))
  )
  cells <- rails_cells(d, c("a", "b"))

  expect_s3_class(cells, "rails_cells")
  expect_equal(nrow(cells), 3L)
  expect_equal(sum(cells$weight), nrow(d))
  expect_equal(sort(cells$weight), c(1, 2, 2))
})

test_that("row_cell maps every record back to its cell", {
  d <- data.frame(a = factor(c("x", "y", "x")), b = factor(c(1, 1, 1)))
  cells <- rails_cells(d, c("a", "b"))
  idx <- attr(cells, "row_cell")

  expect_equal(length(idx), nrow(d))
  expect_equal(as.character(cells$a[idx]), as.character(d$a))
})

test_that("weights are summed, not counted, when supplied", {
  d <- data.frame(a = factor(c("x", "x", "y")), w = c(2, 3, 10))
  cells <- rails_cells(d, "a", weights = "w")
  expect_equal(cells$weight[cells$a == "x"], 5)
  expect_equal(cells$weight[cells$a == "y"], 10)

  ## A vector and a column name must agree.
  expect_equal(cells$weight, rails_cells(d, "a", weights = d$w)$weight)
})

test_that("factor levels survive aggregation of a subset", {
  d <- data.frame(a = factor(c("x", "y", "z")))
  cells <- rails_cells(d[d$a != "z", , drop = FALSE], "a")
  expect_equal(levels(cells$a), c("x", "y", "z"))
})

test_that("missing covariates are dropped with a warning", {
  d <- data.frame(a = factor(c("x", NA, "y")))
  expect_warning(cells <- rails_cells(d, "a"), "dropped 1 row")
  expect_equal(nrow(cells), 2L)
  expect_true(is.na(attr(cells, "row_cell")[2]))
})

test_that("rails_totals names match a calibration model matrix", {
  pop <- rails_simulate(400, seed = 11)
  cells <- rails_cells(pop, c("agegroup", "sex"))
  tot <- rails_totals(cells, order = 2)

  mm <- stats::model.matrix(~ agegroup + sex + agegroup:sex, pop)
  expect_true(all(colnames(mm) %in% names(tot)))
  ## One-way margins must equal the raw counts.
  expect_equal(unname(tot["(Intercept)"]), nrow(pop))
})

test_that("informative errors on malformed input", {
  d <- data.frame(a = 1:3)
  expect_error(rails_cells(d, "nope"), "missing the variable")
  expect_error(rails_totals(d, vars = "a"), "no `weight` column")
})
