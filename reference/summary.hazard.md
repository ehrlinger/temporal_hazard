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
#>   log-lik:      -53.9889 
#>   evaluations: fn=18, gr=5
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.4490094 0.08664337 5.182271 2.192003e-07
#> nu 1.0022854 0.13850336 7.236542 4.602695e-13

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
#>   phase 2:      late - cdf (early risk)
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -392.666 
#>   evaluations: fn=83, gr=35
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                estimate std_error      z_stat   p_value
#>   log_mu     -0.8392255 10.016157 -0.08378718 0.9332256
#>   log_t_half  1.0149845 14.579562  0.06961694 0.9444986
#>   nu          1.6897930  8.365820  0.20198772 0.8399263
#>   m          -0.1753914  1.087213 -0.16132203 0.8718398
#> 
#>   Phase: late (cdf)
#>                estimate std_error     z_stat   p_value
#>   log_mu     2.86567863 2.2399311 1.27936015 0.2007703
#>   log_t_half 5.41924041 3.6513364 1.48417999 0.1377612
#>   nu         2.27025536 1.5995548 1.41930456 0.1558102
#>   m          0.02566328 0.6731866 0.03812208 0.9695903
# }
```
