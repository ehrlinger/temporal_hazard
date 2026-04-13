# Predict from a hazard model object

Produces prediction outputs from a `hazard` object. Supports multiple
prediction types including linear predictor, hazard, survival
probability, and cumulative hazard.

## Usage

``` r
# S3 method for class 'hazard'
predict(
  object,
  newdata = NULL,
  type = c("hazard", "linear_predictor", "survival", "cumulative_hazard"),
  decompose = FALSE,
  ...
)
```

## Arguments

- object:

  A `hazard` object.

- newdata:

  Optional matrix or data frame of predictors. For types requiring time
  (e.g., "survival", "cumulative_hazard"), newdata should include a
  `time` column, or time will be taken from the fitted object's data.

- type:

  Prediction type:

  - `"linear_predictor"`: Linear predictor η = x·β (not available for
    multiphase)

  - `"hazard"`: Hazard scale exp(η) (not available for multiphase)

  - `"survival"`: Survival probability S(t\|x) = exp(-H(t\|x))

  - `"cumulative_hazard"`: Cumulative hazard H(t\|x) at event times

- decompose:

  Logical; if `TRUE` and the model is multiphase, return a data frame
  with per-phase cumulative hazard contributions alongside the total.
  Ignored for single-distribution models. Default `FALSE`.

- ...:

  Unused; included for S3 compatibility.

## Value

Numeric vector of predictions.

## Details

For Weibull models with survival or cumulative_hazard predictions:

- Cumulative hazard: H(t\|x) = (μ·t)^ν · exp(η)

- Survival: S(t\|x) = exp(-H(t\|x))

Time values must be positive and finite. If newdata contains a `time`
column, it will be used; otherwise, the time vector from the fitted
object is used. For models fit with `time_windows`, predictions for
`type = "linear_predictor"` or `"hazard"` also require time values (via
`newdata$time` or fitted-time fallback) so window-specific coefficients
can be selected.

## See also

[`hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/hazard.md)
for model fitting,
[`summary.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/summary.hazard.md)
for model summaries,
[`hzr_phase()`](https://ehrlinger.github.io/temporal_hazard/reference/hzr_phase.md)
for multiphase temporal shapes.

[`vignette("prediction-visualization")`](https://ehrlinger.github.io/temporal_hazard/articles/prediction-visualization.md)
for detailed prediction workflows including decomposed hazard plots and
patient-specific curves.

## Examples

``` r
# ── Basic predictions ────────────────────────────────────────────────
set.seed(1)
fit <- hazard(time = rexp(50, 0.3), status = rep(1L, 50),
              theta = c(0.3, 1.0), dist = "weibull", fit = TRUE)
predict(fit, type = "survival")
#>  [1] 0.508448836 0.315117497 0.909632142 0.913801355 0.704022509 0.034429522
#>  [7] 0.297904379 0.635876619 0.407732732 0.908685617 0.245861859 0.504733067
#> [13] 0.295096862 0.003729697 0.364948062 0.373066533 0.134498594 0.565324696
#> [19] 0.772692686 0.605268524 0.070992893 0.572924045 0.803140735 0.619329758
#> [25] 0.937239852 0.968075472 0.611315229 0.007471798 0.318195648 0.389681823
#> [31] 0.232969073 0.981596565 0.781843506 0.267479469 0.868327070 0.378412625
#> [37] 0.797695039 0.524953194 0.510431884 0.845615583 0.354514703 0.376046945
#> [43] 0.276615587 0.289749395 0.626388405 0.798022001 0.276332062 0.390676459
#> [49] 0.652266897 0.113525051
predict(fit, newdata = data.frame(time = c(1, 2, 5)),
        type = "cumulative_hazard")
#> [1] 0.2244716 0.5138491 1.5356737

# ── Patient-specific survival curves ─────────────────────────────────
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
  data  = dat,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist  = "weibull", fit = TRUE
)

new_patients <- data.frame(
  time = c(0.5, 1.5, 3.0),
  age  = c(50, 65, 75),
  nyha = c(1, 3, 4),
  shock = c(0, 0, 1)
)
# Compute predictions from the clean covariate frame before adding columns
surv   <- predict(fit2, newdata = new_patients, type = "survival")
cumhaz <- predict(fit2, newdata = new_patients, type = "cumulative_hazard")
new_patients$survival          <- surv
new_patients$cumulative_hazard <- cumhaz
new_patients
#>   time age nyha shock  survival cumulative_hazard
#> 1  0.5  50    1     0 0.9494097        0.05191482
#> 2  1.5  65    3     0 0.7743185        0.25577202
#> 3  3.0  75    4     1 0.5046535        0.68388317

# \donttest{
# ── Grouped survival curves ───────────────────────────────────────
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  t_grid <- seq(0.05, max(dat$time), length.out = 80)
  profiles <- data.frame(
    label = c("Low risk (age 50, NYHA I)",
              "High risk (age 75, NYHA IV)"),
    age   = c(50, 75),
    nyha  = c(1, 4),
    shock = c(0, 1)
  )

  curve_list <- lapply(seq_len(nrow(profiles)), function(i) {
    nd <- data.frame(
      time  = t_grid,
      age   = profiles$age[i],
      nyha  = profiles$nyha[i],
      shock = profiles$shock[i]
    )
    nd$survival <- predict(fit2, newdata = nd, type = "survival") * 100
    nd$profile  <- profiles$label[i]
    nd
  })
  curve_df <- do.call(rbind, curve_list)

  ggplot(curve_df, aes(time, survival, colour = profile)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Months after surgery",
         y = "Freedom from death (%)",
         title = "Predicted survival by risk profile",
         colour = NULL) +
    theme_minimal()
}

# }

# \donttest{
# ── Multiphase predictions with decomposition ────────────────────
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

t_grid <- seq(0.01, max(dat$time) * 0.9, length.out = 100)
nd     <- data.frame(time = t_grid)

# Overall survival
predict(fit_mp, newdata = nd, type = "survival")
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.99998987   0.96791010   0.92109875   0.88058438   0.84325372   0.80741007 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.77242464   0.73813734   0.70457964   0.67185110   0.64006400   0.60931921 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.57969679   0.55125376   0.52402520   0.49802681   0.47325795   0.44970467 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.42734250   0.40613894   0.38605555   0.36704972   0.34907610   0.33208776 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.31603717   0.30087687   0.28656009   0.27304112   0.26027567   0.24822107 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.23683642   0.22608267   0.21592271   0.20632131   0.19724516   0.18866281 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.18054462   0.17286267   0.16559074   0.15870420   0.15217991   0.14599621 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.14013276   0.13457053   0.12929166   0.12427946   0.11951828   0.11499348 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.11069132   0.10659897   0.10270441   0.09899635   0.09546424   0.09209820 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.08888895   0.08582780   0.08290660   0.08011770   0.07745394   0.07490857 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.07247527   0.07014811   0.06792150   0.06579019   0.06374926   0.06179405 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.05992021   0.05812361   0.05640039   0.05474689   0.05315968   0.05163550 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.05017130   0.04876418   0.04741143   0.04611045   0.04485883   0.04365425 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.04249453   0.04137763   0.04030158   0.03926455   0.03826478   0.03730062 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.03637049   0.03547291   0.03460647   0.03376982   0.03296168   0.03218086 
#> early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.03142620   0.03069660   0.02999103   0.02930851   0.02864808   0.02800886 
#> early.log_mu early.log_mu early.log_mu early.log_mu 
#>   0.02739000   0.02679068   0.02621012   0.02564760 

# Per-phase decomposed cumulative hazard
decomp <- predict(fit_mp, newdata = nd,
                  type = "cumulative_hazard", decompose = TRUE)
head(decomp)
#>        time        total        early         late
#> 1 0.0100000 1.013249e-05 3.595706e-14 1.013249e-05
#> 2 0.3177112 3.261607e-02 2.637804e-02 6.238033e-03
#> 3 0.6254224 8.218803e-02 6.104339e-02 2.114464e-02
#> 4 0.9331336 1.271695e-01 8.444137e-02 4.272814e-02
#> 5 1.2408448 1.704874e-01 1.007901e-01 6.969728e-02
#> 6 1.5485560 2.139236e-01 1.128482e-01 1.010754e-01
# }
```
