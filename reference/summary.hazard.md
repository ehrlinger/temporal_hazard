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
#>   log-lik:      -48.3228 
#>   evaluations: fn=15, gr=6
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.5201972 0.09124047 5.701386 1.188371e-08
#> nu 1.0961938 0.16189266 6.771115 1.277937e-11

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
#> Warning: NaNs produced
#> Multiphase hazard model (2 phases)
#>   observations: 200 
#>   predictors:   0 
#>   dist:         multiphase 
#>   phase 1:      early - cdf (early risk)
#>   phase 2:      late - cdf (late risk)
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -388.918 
#>   evaluations: fn=100, gr=32
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                  estimate   std_error      z_stat      p_value
#>   log_mu     -3.979438679 0.569883164   -6.982903 2.891427e-12
#>   log_t_half -0.943134742 0.007166541 -131.602509 0.000000e+00
#>   nu          0.003751501         NaN          NA           NA
#>   m           4.152983242         NaN          NA           NA
#> 
#>   Phase: late (cdf)
#>                estimate std_error     z_stat      p_value
#>   log_mu      2.7627013 1.1630799  2.3753324 1.753315e-02
#>   log_t_half  5.0259890 2.8740152  1.7487691 8.033094e-02
#>   nu          2.6355015 3.4025002  0.7745779 4.385891e-01
#>   m          -0.7663702 0.1108225 -6.9152967 4.668851e-12
# }
```
