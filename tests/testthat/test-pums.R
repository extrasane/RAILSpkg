## No test here touches the network: pums_recode() and pums_reference() are
## pure, and pums_fetch() is only exercised for its input validation.

raw_records <- function() {
  data.frame(
    PWGTP   = c(100, 250,  90, 300, 120,  75),
    AGEP    = c( 30,  70,  16,  50,  22,  80),   # one minor, to be filtered
    SEX     = c(  1,   2,   2,   2,   1,   1),
    HINCP   = c(42000, 120000, 20000, 60000, 90000, -70000),
    SCHL    = c( 21,  16,  12,  19,   9,  22),
    TEN     = c(  1,   3,   2,   4,   3,   1),
    HISP    = c(  1,   5,   1,   1,   1,   1),
    RACWHT  = c(  1,   0,   0,   0,   1,   1),
    RACBLK  = c(  0,   0,   1,   0,   0,   0),
    RACASN  = c(  0,   0,   0,   1,   0,   0),
    REGION  = c(  1,   4,   3,   2,   3,   1)
  )
}

test_that("pums_recode bands every variable as documented", {
  out <- pums_recode(raw_records())

  expect_equal(nrow(out), 5L)                       # the 16-year-old is gone
  expect_true(all(out$age > 17))

  expect_equal(as.character(out$agegroup),
               c("25-44", "65-74", "45-64", "18-24", "75+"))
  expect_equal(as.character(out$sex),
               c("Male", "Female", "Female", "Male", "Male"))
  expect_equal(as.character(out$homeown),
               c("Own", "Rent", "Others", "Rent", "Own"))
  expect_equal(as.character(out$region),
               c("Northeast", "West", "Midwest", "South", "Northeast"))
})

test_that("HISP is Hispanic for every code except 1", {
  out <- pums_recode(raw_records())
  ## Record 2 has HISP = 5 and so is Hispanic regardless of its race codes.
  expect_equal(as.character(out$race_eth)[2], "Hispanic")
  expect_equal(as.character(out$race_eth)[1], "NH White")
  expect_equal(as.character(out$race_eth)[3], "NH Asian")
})

test_that("income bands are half-open and very negative income is missing", {
  out <- pums_recode(raw_records())
  expect_equal(as.character(out$income),
               c("35k-50k", ">100k", "50k-75k", "75k-100k", NA))

  ## A value exactly on a break belongs to the upper band.
  edge <- raw_records()[1, ]
  edge$HINCP <- 35000
  expect_equal(as.character(pums_recode(edge)$income), "35k-50k")
  edge$HINCP <- 34999
  expect_equal(as.character(pums_recode(edge)$income), "<35k")
})

test_that("education bands follow the SCHL cut points", {
  d <- raw_records()[1, ]
  band <- function(schl) {
    d$SCHL <- schl
    as.character(pums_recode(d)$edu)
  }
  expect_equal(band(11), "Less than highschool")
  expect_equal(band(12), "Some highschool")
  expect_equal(band(16), "Highschool graduate")
  expect_equal(band(17), "Highschool graduate")
  expect_equal(band(18), "Some college")
  expect_equal(band(21), "College graduate or advanced")
})

test_that("factor levels are fixed, so separate extracts line up", {
  a <- pums_recode(raw_records()[1, ])
  b <- pums_recode(raw_records()[2, ])
  for (v in c("agegroup", "sex", "race_eth", "income", "edu", "homeown",
              "region")) {
    expect_identical(levels(a[[v]]), levels(b[[v]]), info = v)
  }
})

test_that("adults_only = FALSE keeps everyone", {
  expect_equal(nrow(pums_recode(raw_records(), adults_only = FALSE)), 6L)
})

test_that("ST is mapped to a USPS code", {
  d <- raw_records()
  d$ST <- c(6, 36, 48, 12, 4, 99)
  out <- suppressWarnings(pums_recode(d))
  expect_equal(out$state, c("CA", "NY", "FL", "AZ", NA))
  expect_warning(pums_recode(d), "outside the 50 states")
})

test_that("a non-2022 vintage warns before it recodes", {
  expect_warning(pums_recode(raw_records(), year = 2019),
                 "codings here follow the 2022")
})

test_that("missing ACS columns are named, not silently tolerated", {
  d <- raw_records()
  d$SCHL <- NULL
  d$TEN  <- NULL
  expect_error(pums_recode(d), "SCHL, TEN")
})

test_that("pums_reference rescales weights to the pre-drop population", {
  d <- data.frame(
    weight   = c(100, 200, 300, 400),
    agegroup = factor(c("18-24", "25-44", "25-44", NA)),
    sex      = factor(c("Male", "Female", "Male", "Female"))
  )

  ref <- pums_reference(d, vars = c("agegroup", "sex"), order = 2)

  expect_s3_class(ref, "pums_reference")
  expect_equal(ref$n_dropped, 1L)
  expect_equal(ref$n_pop, 1000)              # includes the dropped record
  expect_equal(sum(ref$micro$weight), 1000)  # rescaled up from 600
  expect_equal(sum(ref$cells$weight), 1000)
  expect_equal(unname(ref$totals["(Intercept)"]), 1000)
})

test_that("rescale = FALSE leaves the retained weights alone", {
  d <- data.frame(
    weight   = c(100, 200, 300, 400),
    agegroup = factor(c("18-24", "25-44", "25-44", NA)),
    sex      = factor(c("Male", "Female", "Male", "Female"))
  )
  ref <- pums_reference(d, vars = c("agegroup", "sex"), order = 2,
                        rescale = FALSE)
  expect_equal(sum(ref$micro$weight), 600)
  expect_equal(ref$n_pop, 600)
})

test_that("pums_reference output feeds rails() directly", {
  pop <- rails_simulate(20000, seed = 42)
  vars <- c("agegroup", "sex", "income")

  aou  <- pop[pop$aou == 1, ]
  refd <- pop[pop$s == 1, ]
  refd$weight <- 1 / refd$ps

  ref <- pums_reference(refd, vars = vars, order = 3)

  fit <- suppressWarnings(rails(
    rails_cells(aou, vars), ref$cells, vars,
    pop_totals = ref$totals, aggregated = TRUE, verbose = FALSE
  ))

  expect_false(anyNA(weights(fit)))
  ## The weights reproduce the reference margins, which is the whole point.
  got  <- tapply(weights(fit) * rails_cells(aou, vars)$weight,
                 rails_cells(aou, vars)$sex, sum)
  want <- tapply(ref$cells$weight, ref$cells$sex, sum)
  expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-6)
})

test_that("print.pums_reference reports the drop and the rescale", {
  d <- data.frame(
    weight   = c(100, 200, 300, 400),
    agegroup = factor(c("18-24", "25-44", "25-44", NA)),
    sex      = factor(c("Male", "Female", "Male", "Female"))
  )
  ref <- pums_reference(d, vars = c("agegroup", "sex"), order = 2)
  expect_output(print(ref), "dropped for missing covariates")
  expect_output(print(ref), "population : 1,000")
})

test_that("pums_fetch validates before it reaches the network", {
  expect_error(pums_fetch(key = ""), "No Census API key")
  expect_error(pums_fetch(geography = "county"), "'arg' should be one of")
})

test_that("pums_variables documents what it fetches", {
  v <- pums_variables()
  expect_true(all(c("PWGTP", "AGEP", "SEX", "HINCP", "SCHL", "TEN", "HISP",
                    "REGION") %in% names(v)))
  expect_false("ST" %in% names(pums_variables("region")))
  expect_true("ST" %in% names(pums_variables("state")))
})
