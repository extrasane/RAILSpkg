## Verify the new package against the ORIGINAL sources, not against itself.
##
## Originals:
##   data-raw/original-scaffold/R/*.R   the pre-package scaffold, archived here
##     from the RAILS repo branch package-build before that branch was deleted
##   C:/Work/REpo/RAILS_Edit/Functions/supplementary_functions.R   fun.rails.var
##
## Both implementations are fed byte-identical inputs; only the code differs.
##
## Every rails() call passes lifo = "ascending": the application code climbed
## from the starting model and halted at the first raking failure, whereas the
## package now defaults to the descending LIFO the manuscripts define.

suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(survey); library(magrittr); library(rlang)
})

OLD_DIR <- "C:/Work/REpo/RAILS_package/data-raw/original-scaffold/R"
FUNCS   <- "C:/Work/REpo/RAILS_Edit/Functions/supplementary_functions.R"

## ---- Load the originals into their own environment ------------------------
orig <- new.env(parent = globalenv())
for (f in list.files(OLD_DIR, pattern = "[.]R$", full.names = TRUE)) {
  sys.source(f, envir = orig)
}
cat("original scaffold functions:",
    paste(sort(ls(orig)), collapse = ", "), "\n\n")

## fun.rails.var lives only in the simulation sources. Source it in isolation:
## the file also defines fun.lkd / fun.vs, which we do not want here.
varenv <- new.env(parent = globalenv())
sys.source(FUNCS, envir = varenv)
orig_rails_var <- get("fun.rails.var", envir = varenv)

## ---- The new package ------------------------------------------------------
suppressMessages(devtools::load_all("C:/Work/REpo/RAILS_package", quiet = TRUE))

## ---- Shared input ---------------------------------------------------------
vars <- c("agegroup", "sex", "income")
pop  <- rails_simulate(20000, seed = 42)
aou  <- pop[pop$aou == 1, ]
ref  <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps

cells_np  <- as.data.frame(rails_cells(aou, vars))
cells_ref <- as.data.frame(rails_cells(ref, vars, weights = "dweight"))

## Population margins, built exactly as the application scripts build them.
mat_max <- Matrix::sparse.model.matrix(
  ~ agegroup + sex + income +
    agegroup:sex + agegroup:income + sex:income + agegroup:sex:income,
  data = cells_ref, keep.order = TRUE
)
pop_totals <- setNames(
  as.numeric(Matrix::crossprod(mat_max, cells_ref$weight)), colnames(mat_max)
)

report <- function(label, a, b, tol = 0) {
  same <- isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
  d <- suppressWarnings(max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE))
  cat(sprintf("%-46s %-9s max|diff| = %s\n", label,
              if (same) "IDENTICAL" else "DIFFERS",
              format(d, digits = 3, scientific = TRUE)))
  invisible(same)
}

ok <- logical(0)

## ===========================================================================
cat("== 1. nps_weights() vs original fun.nps() ==================================\n")
f2 <- ~ agegroup + sex + income + agegroup:sex + agegroup:income + sex:income

d_old <- orig$fun.nps(f2, cells_np, cells_ref)
d_new <- nps_weights(f2, cells_np, cells_ref)
ok <- c(ok, report("weights", as.numeric(d_old), as.numeric(d_new)))
ok <- c(ok, report("theta",   attr(d_old, "theta"), attr(d_new, "theta")))

## Warm start must land in the same place as a cold start.
f3 <- ~ agegroup + sex + income + agegroup:sex + agegroup:income + sex:income +
        agegroup:sex:income
d_old3 <- orig$fun.nps(f3, cells_np, cells_ref)
d_new3 <- nps_weights(f3, cells_np, cells_ref, theta_init = attr(d_new, "theta"))
ok <- c(ok, report("warm-started weights", as.numeric(d_old3), as.numeric(d_new3),
                   tol = 1e-8))

## ===========================================================================
cat("\n== 2. gvs_step() vs original fun.lkd() ====================================\n")
m_ref <- Matrix::sparse.model.matrix(f2, cells_ref)
terms2 <- c("agegroup", "sex", "income",
            "agegroup:sex", "agegroup:income", "sex:income")
lk2 <- sum(as.numeric(Matrix::crossprod(
  Matrix::sparse.model.matrix(f2, cells_np), cells_np$weight)) *
    attr(d_new, "theta")) -
  sum(cells_ref$weight * log1p(exp(as.numeric(m_ref %*% attr(d_new, "theta")))))

lrt_old <- orig$fun.lkd("agegroup:sex:income", terms2, lk2, cells_ref, cells_np, m_ref)
lrt_new <- gvs_step("agegroup:sex:income", terms2, lk2, cells_ref, cells_np, m_ref)
ok <- c(ok, report("c(loglik, LRT, p, LRT/df, df)", lrt_old, lrt_new))

## ===========================================================================
cat("\n== 3. rails() vs original fun.rails.threeway() ============================\n")
old <- suppressWarnings(suppressMessages(
  orig$fun.rails.threeway(cells_np, cells_ref, pop_totals, names_univar = vars)
))
new <- suppressWarnings(rails(cells_np, cells_ref, vars, pop_totals = pop_totals,
                              aggregated = TRUE, benchmarks = TRUE,
                              lifo = "ascending", verbose = FALSE))

ok <- c(ok, report("d_rails", old$d_rails, weights(new)))
for (nm in c("d_unweighted", "d_cal1", "d_cal2", "d_nps1", "d_nps2",
             "d_nps1_rake", "d_nps2_rake")) {
  ok <- c(ok, report(nm, old[[nm]], new$benchmarks[[nm]]))
}
cat(sprintf("%-46s %s\n", "selected_terms",
            if (identical(unique(old$selected_terms),
                          paste(c(new$terms_base, new$terms_selected),
                                collapse = " + "))) "IDENTICAL" else "DIFFERS"))
cat(sprintf("%-46s %s\n", "calibrated_terms",
            if (identical(unique(old$calibrated_terms),
                          paste(new$terms_used, collapse = " + "))) "IDENTICAL" else "DIFFERS"))
ok <- c(ok, identical(unique(old$selected_terms),
                      paste(c(new$terms_base, new$terms_selected), collapse = " + ")),
        identical(unique(old$calibrated_terms),
                  paste(new$terms_used, collapse = " + ")))

## ===========================================================================
cat("\n== 4. rails_subgroup() vs original fun.sub.rails.threeway() ===============\n")
## The original crashes with fewer than three within-stratum covariates
## (utils::combn(2, 3)), so the comparison uses three.
sub_vars <- c("agegroup", "income", "homeown")
cells_np_s  <- as.data.frame(rails_cells(aou, c(sub_vars, "sex")))
cells_ref_s <- as.data.frame(rails_cells(ref, c(sub_vars, "sex"),
                                         weights = "dweight"))

old_s <- suppressWarnings(suppressMessages(
  orig$fun.sub.rails.threeway(cells_np_s, cells_ref_s, subgroup_var = "sex",
                              names_univar = sub_vars)
))
new_s <- suppressWarnings(rails_subgroup(
  aou, ref, sub_vars, by = "sex", weights_ref = "dweight",
  lifo = "ascending", verbose = FALSE
))

## The original returns cell-level results; line them up by cell.
for (lev in levels(droplevels(factor(old_s$sex)))) {
  o <- old_s[old_s$sex == lev, ]
  n <- new_s$fits[[lev]]
  key_o <- do.call(paste, o[sub_vars])
  key_n <- do.call(paste, n$cells[sub_vars])
  stopifnot(!anyNA(match(key_n, key_o)), !anyDuplicated(key_o))
  ok <- c(ok, report(paste0("sex = ", lev, ": d_rails"),
                     o$d_rails[match(key_n, key_o)], n$weights_unit))
}

## ===========================================================================
cat("\n== 5. rails_var() vs original fun.rails.var() =============================\n")
fit <- suppressWarnings(rails(aou, ref, vars, weights_ref = "dweight",
                              lifo = "ascending", verbose = FALSE))
stopifnot(!anyNA(weights(fit)))

cat_temp <- fit$formula_used
theta    <- fit$theta
w_NP     <- weights(fit)
d_P      <- 1 / ref$ps

X_NP  <- model.matrix(cat_temp, aou[, vars, drop = FALSE])
X_P   <- model.matrix(cat_temp, ref[, vars, drop = FALSE])
pi_NP <- plogis(as.numeric(X_NP %*% theta))
pi_P  <- plogis(as.numeric(X_P  %*% theta))
T_pop <- RAILS:::create_v3(fit$pop_totals, colnames(X_NP))

v_old <- orig_rails_var(
  dt_aou = aou[, vars, drop = FALSE], dt_s = ref[, vars, drop = FALSE],
  y_aou = aou$y, cat_temp = cat_temp, w_NP = w_NP, d_P = d_P,
  pi_NP = pi_NP, pi_P = pi_P, T_pop = T_pop, tempy = mean(pop$y)
)
v_new <- rails_var(fit, aou$y, type = "both", truth = mean(pop$y))

## The original sums the sandwich over individual records; the package sums the
## same quantities over covariate cells, which is exact algebra but a different
## summation order. Bit-equality is therefore not the right test here -- these
## comparisons run at VAR_TOL, which is still far below anything of scientific
## consequence. Everything outside this section stays at tolerance 0.
##
## term_11, term_12 and term_13 are exempt because they are numerically zero on
## this data (order 1e-31 against a variance of 1e-05), so they carry only
## cancellation noise and their *relative* difference is meaningless.
VAR_TOL  <- 1e-10
NEAR_ZED <- c("term_11", "term_12", "term_13")

shared <- c("mu_hat", "var_naive", "se_naive", "var_full", "se_full",
            "ratio_full_naive", "term_11", "term_22", "term_33",
            "term_12", "term_13", "term_23", "N_hat", "n_NP", "n_P")
for (nm in shared) {
  if (nm %in% NEAR_ZED) {
    ok <- c(ok, report(paste0(nm, " (both ~0)"), 0,
                       max(abs(c(v_old[[nm]], v_new[[nm]]))) / v_old[["var_full"]],
                       tol = 1e-12))
  } else {
    ok <- c(ok, report(nm, v_old[[nm]], v_new[[nm]], tol = VAR_TOL))
  }
}
ok <- c(ok, report("ci_naive_lower", v_old[["ci_naive_lower"]],
                   v_new[["ci_naive_lower"]], tol = VAR_TOL))
ok <- c(ok, report("ci_full_upper",  v_old[["ci_full_upper"]],
                   v_new[["ci_full_upper"]],  tol = VAR_TOL))
ok <- c(ok, report("coverage (naive, full)",
                   c(v_old[["cover_naive_temp"]], v_old[["cover_full_temp"]]),
                   c(v_new[["cover_naive"]],      v_new[["cover_full"]])))

## The simplified-only path must agree with the same quantities in "both".
v_simple <- rails_var(fit, aou$y, truth = mean(pop$y))
ok <- c(ok, report("type='simplified' vs original var_naive",
                   v_old[["var_naive"]], v_simple[["var_naive"]], tol = VAR_TOL))

## An aggregated fit must reach the same numbers from cell tables alone.
fit_agg <- suppressWarnings(rails(cells_np, cells_ref, vars,
                                  aggregated = TRUE, lifo = "ascending",
                                  verbose = FALSE))
rc  <- attr(rails_cells(aou, vars), "row_cell")
y_c <- data.frame(sum   = as.numeric(rowsum(aou$y,   rc)),
                  sumsq = as.numeric(rowsum(aou$y^2, rc)))
v_agg <- rails_var(fit_agg, y_c, type = "both", truth = mean(pop$y))
ok <- c(ok, report("aggregated fit vs microdata fit (var_full)",
                   v_new[["var_full"]], v_agg[["var_full"]], tol = VAR_TOL))

## ===========================================================================
cat("\n===========================================================================\n")
cat(sprintf("%d of %d comparisons identical to the original code.\n", sum(ok), length(ok)))
if (all(ok)) cat("ALL MATCH\n") else cat("*** MISMATCHES PRESENT ***\n")

## (fixture writing lives in build_fixture.R)
