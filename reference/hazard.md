# Build and optionally fit a hazard model

Creates a `hazard` object and optionally fits it via maximum likelihood.
This mirrors the argument-oriented workflow of the legacy HAZARD C/SAS
implementation: supply starting values in `theta` and the function will
optimize to produce fitted estimates.

## Usage

``` r
hazard(
  formula = NULL,
  data = NULL,
  time = NULL,
  status = NULL,
  time_lower = NULL,
  time_upper = NULL,
  x = NULL,
  time_windows = NULL,
  theta = NULL,
  dist = "weibull",
  fit = FALSE,
  control = list(),
  ...
)
```

## Arguments

- formula:

  Optional formula of the form `Surv(time, status) ~ predictors`. When
  provided, overrides direct time/status/x arguments and extracts from
  data. Example:
  `hazard(Surv(time, status) ~ x1 + x2, data = df, dist = "weibull", fit = TRUE)`.

- data:

  Optional data frame containing variables referenced in formula.

- time:

  Numeric follow-up time vector.

- status:

  Numeric or logical event indicator vector.

- time_lower:

  Optional numeric lower bound vector for censoring intervals. Used when
  `status == 2` (interval-censored); defaults to `time` if NULL.

- time_upper:

  Optional numeric upper bound vector for censoring intervals. Used when
  `status %in% c(-1, 2)`; defaults to `time` if NULL.

- x:

  Optional design matrix (or data frame coercible to matrix).

- time_windows:

  Optional numeric vector of strictly positive cut points for piecewise
  time-varying coefficients. When provided, each predictor column in `x`
  is expanded into one column per time window so each window gets its
  own coefficient.

- theta:

  Optional numeric coefficient vector (starting values for
  optimization).

- dist:

  Character baseline distribution label (default "weibull").

- fit:

  Logical; if TRUE and theta is provided, fit the model via ML (default
  TRUE).

- control:

  Named list of control options (see Details).

- ...:

  Additional named arguments retained for parity with legacy calling
  conventions.

## Details

Control parameters:

- `maxit`: Maximum iterations (default 1000)

- `reltol`: Relative parameter change tolerance (default 1e-5)

- `abstol`: Absolute gradient norm tolerance (default 1e-6)

- `method`: Optimization method: "bfgs" or "nm" (default "bfgs")

- `condition`: Condition number control (default 14)

- `nocov`, `nocor`: Suppress covariance/correlation output (legacy;
  no-op in M2)

Censoring status coding:

- 1: Exact event at time

- 0: Right-censored at time

- -1: Left-censored with upper bound at time_upper \\or time\\

- 2: Interval-censored in the interval \\time_lower, time_upper\\

Time-varying coefficients:

- If `time_windows` is supplied, predictors are expanded to piecewise
  window interactions so each window has its own coefficient vector.

- This is implemented as design-matrix expansion, so the existing
  likelihood engines remain unchanged.

## Examples

``` r
# ── Univariable Weibull ──────────────────────────────────────────────
set.seed(1)
time   <- rexp(50, rate = 0.3)
status <- sample(0:1, 50, replace = TRUE, prob = c(0.3, 0.7))
fit <- hazard(time = time, status = status,
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
summary(fit)
#> hazard model summary
#>   observations: 50 
#>   predictors:   0 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -89.5173 
#>   evaluations: fn=23, gr=7
#> 
#> Coefficients:
#>     estimate  std_error   z_stat      p_value
#> mu 0.2190092 0.02728228 8.027524 9.945983e-16
#> nu 1.3389540 0.15829050 8.458840 2.700510e-17

# ── Formula interface with covariates ────────────────────────────────
set.seed(1001)
n   <- 180
dat <- data.frame(
  time   = rexp(n, rate = 0.35) + 0.05,
  status = rbinom(n, size = 1, prob = 0.6),
  age    = rnorm(n, mean = 62, sd = 11),
  nyha   = sample(1:4, n, replace = TRUE),
  shock  = rbinom(n, size = 1, prob = 0.18)
)

fit2 <- hazard(
  survival::Surv(time, status) ~ age + nyha + shock,
  data    = dat,
  theta   = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist    = "weibull",
  fit     = TRUE,
  control = list(maxit = 300)
)
summary(fit2)
#> hazard model summary
#>   observations: 180 
#>   predictors:   3 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   converged:    TRUE 
#>   log-lik:      -293.553 
#>   evaluations: fn=36, gr=8
#> 
#> Coefficients:
#>          estimate   std_error      z_stat      p_value
#> mu    0.121860518 0.062260594  1.95726558 5.031625e-02
#> nu    1.143730632 0.084302528 13.56697905 6.285984e-42
#> beta1 0.001716551 0.008807986  0.19488579 8.454824e-01
#> beta2 0.156208724 0.090597558  1.72420457 8.467092e-02
#> beta3 0.017352702 0.362918437  0.04781433 9.618642e-01

# \donttest{
# ── Parametric survival with Kaplan-Meier overlay (requires hvtiPlotR)
if (requireNamespace("hvtiPlotR", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(hvtiPlotR)

  # Parametric curve on a fine grid at median covariate profile
  t_grid   <- seq(0.05, max(dat$time), length.out = 80)
  curve_df <- data.frame(
    time = t_grid, age = median(dat$age), nyha = 2, shock = 0
  )
  curve_df$survival <- predict(fit2, newdata = curve_df, type = "survival")

  # Kaplan-Meier empirical overlay
  km    <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
  km_df <- data.frame(time = km$time, estimate = km$surv * 100)

  hz_obj <- hv_hazard(
    curve_data       = transform(curve_df, survival = survival * 100),
    x_col            = "time",
    estimate_col     = "survival",
    empirical        = km_df,
    emp_x_col        = "time",
    emp_estimate_col = "estimate"
  )

  plot(hz_obj) +
    ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(
      x     = "Years after surgery",
      y     = "Freedom from death (%)",
      title = "Parametric survival vs Kaplan-Meier"
    ) +
    hv_theme_manuscript()
}

# }
```
