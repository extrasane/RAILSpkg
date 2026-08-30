# Subgroup RAILS

``` r

library(RAILS)
```

## Why not just subset

Suppose the question is not national prevalence but prevalence *by
region*, or *by sex*. The obvious move is to fit
[`rails()`](https://extrasane.github.io/RAILSpkg/reference/rails.md)
once on the whole sample and then sum the weights within each subgroup.

That is not the same estimator, and it is weaker in two ways. A single
global fit selects one set of interactions for the whole population, so
a term that matters only in one stratum has to be strong enough to
survive a test run on everyone. And a global fit rakes to national
margins, which constrains the weights to reproduce national totals — not
each subgroup’s own totals.

[`rails_subgroup()`](https://extrasane.github.io/RAILSpkg/reference/rails_subgroup.md)
instead runs the whole procedure independently within each level: its
own propensity model, its own selection path, its own stepwise walk, and
its own population margins built from that stratum’s reference records.

## Fitting

``` r

pop <- rails_simulate(20000, seed = 7)

aou <- pop[pop$aou == 1, ]
ref <- pop[pop$s == 1, ]
ref$dweight <- 1 / ref$ps
```

`by` names the stratifier. It must be a column of both samples, and it
must not also appear in `vars` — within a stratum it is constant, so it
cannot be a covariate there.

### Choosing a stratifier

Which variable you stratify by is worth a moment’s thought, because it
decides what is left for the selection step to find.

The reason is that conditioning removes structure. Whatever you hold
fixed can no longer interact with anything: an interaction between the
stratifier and another variable does not disappear, it simply reappears
as a main effect inside each stratum, where there is nothing left to
detect.

That matters here because this synthetic population is built with two
interactions driving who volunteers — age with sex, and sex with income
— and nothing else. Hold sex fixed and both of them flatten: age×sex
becomes a plain age effect, sex×income becomes a plain income effect,
and the only candidate the forward test could still consider inside a
stratum, age×income, was never in the population to begin with. Every
stratum would correctly select nothing, and the vignette would
demonstrate nothing.

Hold *age* fixed instead and sex and income both still vary, with the
sex×income term intact in each stratum. So that is what the fit below
does: stratify by `agegroup`, and let the selection step look for
interactions among `sex` and `income`.

None of this is a limitation of the method — it is a property of how
this particular example was simulated. But the same question is worth
asking of real data: if the structure you care about runs through the
variable you are stratifying by, you will not find it inside the strata.

``` r

fit <- rails_subgroup(
  np   = aou,
  ref  = ref,
  vars = c("sex", "income"),
  by   = "agegroup",
  weights_ref = "dweight",
  start = 1,          # as in the getting-started vignette
  scope = 2,
  fallback = TRUE,    # explained below
  verbose = FALSE
)

fit
#> Subgroup RAILS fit
#>   stratifier : agegroup (3 of 3 level(s) fitted)
#>   covariates : sex, income
#> 
#>   agegroup = young: 1 selected, 1 raked, 8 cells
#>   agegroup = middle: 0 selected, 0 raked, 8 cells
#>   agegroup = older: 1 selected, 1 raked, 8 cells
#> 
#>   weights: n = 3826, sum = 19,944.8
```

Note in the printout that the strata need not select the same terms.
That is the point of the method: the covariate structure driving
volunteer participation can genuinely differ between subpopulations, and
forcing one model on all of them gives up that information. Here two
strata pick up `sex:income` and the third does not — the term is real,
but not every stratum has the sample size to detect it.

``` r

lapply(fit$fits, function(f) f$terms_selected)
#> $young
#> [1] "sex:income"
#> 
#> $middle
#> character(0)
#> 
#> $older
#> [1] "sex:income"
```

### What happens when a stratum selects nothing

That third stratum is why `fallback = TRUE` appears in the call above,
and it is worth understanding rather than copying.

RAILS builds its weights in two stages: selection proposes terms, then
those terms are raked in one at a time. If selection proposes nothing —
or if the very first raking step fails to converge — there is no
higher-order model to hand back, and something has to be returned
anyway.

By default RAILS returns `NA`, which is what the published code does.
That is a deliberate choice: it makes a non-result impossible to miss,
so nobody reports starting-model weights as though a selection step had
confirmed them.

`fallback = TRUE` returns the raked starting model in that case instead.
The weights are perfectly usable — they are simply the starting model’s
weights, with nothing higher-order added. The two settings on the same
data:

``` r

strict  <- suppressWarnings(rails_subgroup(
  aou, ref, c("sex", "income"), by = "agegroup", weights_ref = "dweight",
  start = 1, scope = 2, verbose = FALSE))            # fallback = FALSE

c(strict  = sum(is.na(weights(strict))),
  relaxed = sum(is.na(weights(fit))))
#>  strict relaxed 
#>    1669       0
```

Either way `converged` records which happened, so the two are always
distinguishable after the fact:

``` r

vapply(fit$fits, function(f) f$converged, logical(1))
#>  young middle  older 
#>   TRUE   TRUE   TRUE
```

In a real analysis, `FALSE` is the safer default and `TRUE` is the
pragmatic one for a pipeline that has to produce a weight for every
record. This vignette uses `TRUE` so that all three strata have weights
to show.

## Using the weights

The returned object carries one pooled weight vector aligned to the rows
of the non-probability sample, so it is used exactly like a global
fit’s:

``` r

w <- weights(fit)

## Each stratum is scaled to its own reference total.
round(rbind(
  reference = tapply(ref$dweight, ref$agegroup, sum),
  rails     = tapply(w, aou$agegroup, sum)
))
#>           young middle older
#> reference  5887   9067  4991
#> rails      5887   9067  4991
```

Subgroup estimates, against the truth we happen to know:

``` r

est <- tapply(seq_along(w), aou$agegroup,
              function(i) sum(w[i] * aou$y[i]) / sum(w[i]))

round(rbind(
  truth      = tapply(pop$y, pop$agegroup, mean),
  unweighted = tapply(aou$y, aou$agegroup, mean),
  rails      = est
), 4)
#>             young middle  older
#> truth      0.1473 0.1984 0.2498
#> unweighted 0.1755 0.2271 0.2961
#> rails      0.1367 0.1929 0.2622
```

The individual `rails_fit` objects are available if a stratum needs
inspecting on its own:

``` r

summary(fit$fits[[1]])
#> RAILS fit summary
#> 
#> Covariates: sex, income
#> Model     : 2 term(s) in the starting model, 1 candidate(s) in scope, alpha 0.05
#> 
#> Selected terms (1):
#>   1. sex:income
#> 
#> Raked to 3 term(s); full walk completed.
#> 
#> Weight diagnostics:
#>                sum                var non_positive_ratio          below_one 
#>          5886.8416            27.2550             0.0000             0.0000 
#>                min                max 
#>             3.5003            23.6906 
#> 
#> Design effect (Kish): 1.444
#> Effective sample size: 520.718
```

## Standard errors

[`rails_var()`](https://extrasane.github.io/RAILSpkg/reference/rails_var.md)
takes a subgroup fit and returns one row per stratum:

``` r

rails_var(fit, aou$y, truth = tapply(pop$y, pop$agegroup, mean))
#>   agegroup    mu_hat    var_naive   se_naive ci_naive_lower ci_naive_upper
#> 1    young 0.1366742 0.0001751015 0.01323259      0.1107388      0.1626097
#> 2   middle 0.1928909 0.0001101003 0.01049287      0.1723252      0.2134565
#> 3    older 0.2622384 0.0001853732 0.01361518      0.2355531      0.2889236
#>   cover_naive    N_hat n_NP
#> 1           1 5886.842  752
#> 2           1 9066.539 1669
#> 3           1 4991.417 1405
```

`truth` may be a single value or, as here, one per stratum named by
level.

The default is the **simplified** variance, and for subgroup RAILS that
default is the one the paper reports. Residual bias persists at the
subgroup level, so a tighter interval would overstate what the estimator
actually delivers.

The **stacked** form is available and does apply within a stratum — each
stratum has its own propensity model, its own selected terms and its own
margins, so the estimating equations really are separate. What the
per-stratum sandwich does not carry is the covariance *between* strata,
and it says so:

``` r

rails_var(fit, aou$y, type = "stacked")[, c("agegroup", "mu_hat", "se_full")]
#> Warning: rails_var(): the stacked variance is computed within each stratum, so
#> it ignores the covariance between strata. Per-stratum intervals are usable;
#> anything pooling or differencing across strata is not.
#>   agegroup    mu_hat    se_full
#> 1    young 0.1366742 0.01283694
#> 2   middle 0.1928909 0.01027202
#> 3    older 0.2622384 0.01338498
```

So those intervals are usable one stratum at a time. Anything that pools
strata, or takes a difference between two of them, needs the
between-stratum covariance that is not computed here — bootstrap the
whole procedure if you need that.

## Empty strata

A level present in the volunteer sample but absent from the reference
sample cannot be fitted: there are no margins to rake to. Those levels
are skipped with a warning and their records receive `NA` weights,
rather than being silently dropped or pooled into a neighbouring
stratum. Check for them before aggregating:

``` r

sum(is.na(w))
#> [1] 0
```
