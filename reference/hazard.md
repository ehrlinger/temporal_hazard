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
# Fit a Weibull hazard model
set.seed(1)
time <- rexp(50, rate = 0.3)
status <- sample(0:1, 50, replace = TRUE, prob = c(0.3, 0.7))
fit <- hazard(time = time, status = status,
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
print(fit)
#> hazard object
#>   observations: 50 
#>   predictors:   0 
#>   dist:         weibull 
#>   engine:       native-r-m2 
#>   log-lik:      9353.55 
#>   converged:    FALSE 

# Formula interface
df <- data.frame(time = time, status = status, x1 = rnorm(50))
fit2 <- hazard(survival::Surv(time, status) ~ x1, data = df,
               theta = c(0.3, 1.0, 0), dist = "weibull", fit = TRUE)
```
