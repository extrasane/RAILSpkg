## Capture the ORIGINAL implementation's output as a test fixture, so the
## package is pinned to it forever without the originals being present at test
## time. Re-run only if the reference implementation itself changes.

suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(survey); library(magrittr); library(rlang)
})

orig <- new.env(parent = globalenv())
for (f in list.files("C:/Work/REpo/RAILS_Edit/RAILS/R", pattern = "[.]R$",
                     full.names = TRUE)) sys.source(f, envir = orig)
varenv <- new.env(parent = globalenv())
sys.source("C:/Work/REpo/RAILS_Edit/Functions/supplementary_functions.R", envir = varenv)

suppressMessages(devtools::load_all("C:/Work/REpo/RAILS_package", quiet = TRUE))

vars <- c("agegroup", "sex", "income")
pop  <- rails_simulate(20000, seed = 42)
aou  <- pop[pop$aou == 1, ]
ref  <- pop[pop$s == 1, ]
d_ref <- 1 / ref$ps

cells_np  <- as.data.frame(rails_cells(aou, vars))
cells_ref <- as.data.frame(rails_cells(ref, vars, weights = d_ref))

mat_max <- Matrix::sparse.model.matrix(
  ~ agegroup + sex + income + agegroup:sex + agegroup:income + sex:income +
    agegroup:sex:income, data = cells_ref, keep.order = TRUE)
pop_totals <- setNames(as.numeric(Matrix::crossprod(mat_max, cells_ref$weight)),
                       colnames(mat_max))

## ---- Global estimator -----------------------------------------------------
g <- suppressWarnings(suppressMessages(
  orig$fun.rails.threeway(cells_np, cells_ref, pop_totals, names_univar = vars)))

## ---- Subgroup estimator ---------------------------------------------------
sub_vars <- c("agegroup", "income", "homeown")
cn_s <- as.data.frame(rails_cells(aou, c(sub_vars, "sex")))
cr_s <- as.data.frame(rails_cells(ref, c(sub_vars, "sex"), weights = d_ref))
s <- suppressWarnings(suppressMessages(
  orig$fun.sub.rails.threeway(cn_s, cr_s, subgroup_var = "sex",
                              names_univar = sub_vars)))

## ---- Variance -------------------------------------------------------------
fit <- suppressWarnings(rails(aou, ref, vars, weights_ref = d_ref, verbose = FALSE))
cat_temp <- fit$formula_used
X_NP <- model.matrix(cat_temp, aou[, vars, drop = FALSE])
X_P  <- model.matrix(cat_temp, ref[, vars, drop = FALSE])
v <- varenv$fun.rails.var(
  dt_aou = aou[, vars, drop = FALSE], dt_s = ref[, vars, drop = FALSE],
  y_aou = aou$y, cat_temp = cat_temp, w_NP = weights(fit), d_P = d_ref,
  pi_NP = plogis(as.numeric(X_NP %*% fit$theta)),
  pi_P  = plogis(as.numeric(X_P  %*% fit$theta)),
  T_pop = RAILS:::create_v3(fit$pop_totals, colnames(X_NP)),
  tempy = mean(pop$y))

fixture <- list(
  vars       = vars,
  cells_np   = cells_np,
  cells_ref  = cells_ref,
  pop_totals = pop_totals,
  global = list(
    d_rails      = g$d_rails,
    d_unweighted = g$d_unweighted, d_cal1 = g$d_cal1, d_cal2 = g$d_cal2,
    d_nps1 = g$d_nps1, d_nps2 = g$d_nps2,
    d_nps1_rake = g$d_nps1_rake, d_nps2_rake = g$d_nps2_rake,
    selected_terms   = unique(g$selected_terms),
    calibrated_terms = unique(g$calibrated_terms)
  ),
  subgroup = list(
    vars = sub_vars, by = "sex",
    cells_np = cn_s, cells_ref = cr_s,
    cells = s[, c(sub_vars, "sex")], d_rails = s$d_rails,
    selected_terms = tapply(s$selected_terms, s$sex, function(z) unique(z)),
    calibrated_terms = tapply(s$calibrated_terms, s$sex, function(z) unique(z))
  ),
  variance = list(
    np    = cbind(aou[, vars, drop = FALSE], y = aou$y),
    ref   = ref[, vars, drop = FALSE],
    d_ref = d_ref,
    truth = mean(pop$y),
    out   = v
  ),
  provenance = paste(
    "Produced by RAILS_Edit/RAILS/R/*.R (fun.rails.threeway,",
    "fun.sub.rails.threeway) and RAILS_Edit/Functions/supplementary_functions.R",
    "(fun.rails.var) on R", getRversion(), "-- the pre-package implementations."
  )
)

path <- "C:/Work/REpo/RAILS_package/tests/testthat/reference-original.rds"
saveRDS(fixture, path, compress = "xz")
cat("wrote", path, "-", round(file.size(path) / 1024, 1), "KB\n")
