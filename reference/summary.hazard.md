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
#>   log-lik:      -48.2997 
#>   evaluations: fn=14, gr=6
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.4903919 0.07632445 6.425096 1.317866e-10
#> nu 1.2365251 0.17655229 7.003733 2.492310e-12

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
#>   log-lik:      -392.68 
#>   evaluations: fn=77, gr=35
#> 
#> Coefficients (internal scale):
#> 
#>   Phase: early (cdf)
#>                estimate std_error    z_stat      p_value
#>   log_mu      2.7838717 0.5125649  5.431257 5.595838e-08
#>   log_t_half  5.3884274       NaN        NA           NA
#>   nu          2.9533613       NaN        NA           NA
#>   m          -0.4742678 0.2193793 -2.161862 3.062879e-02
#> 
#>   Phase: late (cdf)
#>                 estimate std_error      z_stat   p_value
#>   log_mu     -1.34802281 2.2951022 -0.58734762 0.5569703
#>   log_t_half  0.49996554 2.9695810  0.16836232 0.8662983
#>   nu          1.40824517 1.4278706  0.98625543 0.3240078
#>   m          -0.00592305 0.2012198 -0.02943573 0.9765171
# }
```
