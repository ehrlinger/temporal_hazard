# Summarize a hazard model

Returns a compact summary of a `hazard` object, including model
metadata, fit diagnostics, and coefficient-level statistics when
available.

## Usage

``` r
# S3 method for class 'hazard'
summary(object, ...)
```

## Arguments

- object:

  A `hazard` object.

- ...:

  Unused; for S3 compatibility.

## Value

An object of class `summary.hazard`.

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for model fitting,
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
for predictions.

[`vignette("fitting-hazard-models")`](https://ehrlinger.github.io/temporal_hazard/articles/fitting-hazard-models.md)
for fitting workflows,
[`vignette("inference-diagnostics")`](https://ehrlinger.github.io/temporal_hazard/articles/inference-diagnostics.md)
for bootstrap CIs and diagnostics.

## Examples

``` r
# ── Single-phase Weibull summary ────────────────────────────────────
fit <- hazard(time = rexp(30, 0.5), status = rep(1L, 30),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
summary(fit)
#> hazard model summary
#>   observations: 30 
#>   predictors:   0 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -51.3802 
#>   evaluations: fn=14, gr=6
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.4751371 0.08497102 5.591755 2.247861e-08
#> nu 1.0748592 0.15491372 6.938438 3.964594e-12

# \donttest{
# ── Multiphase model summary ────────────────────────────────────────
set.seed(42)
n   <- 200
dat <- data.frame(
  time   = rexp(n, rate = 0.25) + 0.01,
  status = rbinom(n, size = 1, prob = 0.65)
)
fit_mp <- hazard(
  survival::Surv(time, status) ~ 1,
  data   = dat,
  dist   = "multiphase",
  phases = list(
    early = hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0),
    late  = hzr_phase("cdf",      t_half = 5,   nu = 1, m = 0)
  ),
  fit     = TRUE,
  control = list(n_starts = 3, maxit = 500)
)
summary(fit_mp)
#> Multiphase hazard model (2 phases)
#>   observations: 200 
#>   predictors:   0 
#>   dist:         multiphase 
#>   phase 1:      early - cdf (early risk)
#>   phase 2:      late - cdf (late risk)
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -389.632 
#>   evaluations: fn=234, gr=40
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                   estimate std_error z_stat p_value
#>   log_mu     -4.432043e+00        NA     NA      NA
#>   log_t_half -7.237202e-01        NA     NA      NA
#>   nu          1.913460e-03        NA     NA      NA
#>   m           8.278897e-06        NA     NA      NA
#> 
#>   Phase: late (cdf)
#>                estimate std_error z_stat p_value
#>   log_mu      2.4487765        NA     NA      NA
#>   log_t_half  4.2496498        NA     NA      NA
#>   nu          1.6407281        NA     NA      NA
#>   m          -0.7935705        NA     NA      NA
# }
```
