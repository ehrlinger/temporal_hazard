# Example: Hazard Prediction (hp.\*)

``` r
library(TemporalHazard)
library(survival)
library(hvtiPlotR)
```

These examples correspond to the `hp.*` SAS files: prediction from a
fitted hazard model including survival curves, patient-level risk, and
group comparisons. These map to
[`predict.hazard()`](https://ehrlinger.github.io/temporal_hazard/reference/predict.hazard.md)
and `hvtiPlotR` plotting output.

------------------------------------------------------------------------

## hp.death.AVC — Predicted survival after AVC repair

**SAS source:** `examples/hp.death.AVC.sas`

**Purpose:** Compute and plot predicted survival and hazard curves from
the fitted AVC model.

### SAS (original)

``` sas
%HAZPRED(
PROC HAZARD DATA=OUTEST IN=AVCS OUT=HPOUT;
     TIME INT_DEAD;
     PREDICT SURVIVAL HAZARD;
);
```

### R translation (runnable M4 example)

``` r
set.seed(3201)
n <- 240
avcs <- data.frame(
  INT_DEAD = rexp(n, rate = 0.55) + 0.05,
  DEAD = rbinom(n, size = 1, prob = 0.52),
  AGE = rnorm(n, mean = 6, sd = 3),
  STATUS = sample(1:4, n, replace = TRUE),
  MAL = rbinom(n, 1, 0.2)
)

fit_hm <- hazard(
  Surv(INT_DEAD, DEAD) ~ AGE + STATUS + MAL,
  data = avcs,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist = "weibull",
  fit = TRUE,
  control = list(maxit = 300)
)

time_grid <- seq(0.1, 4.0, length.out = 40)
base_profile <- data.frame(
  time = time_grid,
  AGE = median(avcs$AGE),
  STATUS = 2,
  MAL = 0
)

pred_grid <- transform(
  base_profile,
  survival = predict(fit_hm, newdata = base_profile, type = "survival"),
  cumulative_hazard = predict(fit_hm, newdata = base_profile, type = "cumulative_hazard")
)

head(pred_grid)
```

------------------------------------------------------------------------

## hp.death.AVC.hm1 / hp.death.AVC.hm2 — Predictions from two models

**SAS sources:** `examples/hp.death.AVC.hm1.sas`,
`examples/hp.death.AVC.hm2.sas`

**Purpose:** Compare predicted hazard curves from two fitted
multivariable models (e.g., alternative covariate specifications).

### R translation (runnable M4 example)

``` r
fit_hm1 <- hazard(
  Surv(INT_DEAD, DEAD) ~ AGE + STATUS,
  data = avcs,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0),
  dist = "weibull",
  fit = TRUE
)

fit_hm2 <- hazard(
  Surv(INT_DEAD, DEAD) ~ AGE + STATUS + MAL,
  data = avcs,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist = "weibull",
  fit = TRUE
)

compare_grid <- data.frame(
  time = rep(seq(0.1, 4.0, length.out = 25), times = 2),
  AGE = median(avcs$AGE),
  STATUS = 3,
  MAL = rep(c(0, 1), each = 25),
  model = rep(c("Model1", "Model2"), each = 25)
)

compare_grid$survival <- c(
  predict(fit_hm1, newdata = subset(compare_grid, model == "Model1")[, c("time", "AGE", "STATUS")], type = "survival"),
  predict(fit_hm2, newdata = subset(compare_grid, model == "Model2")[, c("time", "AGE", "STATUS", "MAL")], type = "survival")
)

head(compare_grid)
```

------------------------------------------------------------------------

## hp.death.COMPARISON — Model comparison

**SAS source:** `examples/hp.death.COMPARISON.sas`

**Purpose:** Overlay predicted survival/hazard from multiple models or
covariate scenarios on a single plot.

### R translation (runnable M4 example)

``` r
scenario_grid <- rbind(
  data.frame(profile = "low-risk", time = seq(0.1, 4.0, length.out = 20), AGE = quantile(avcs$AGE, 0.25), STATUS = 1, MAL = 0),
  data.frame(profile = "high-risk", time = seq(0.1, 4.0, length.out = 20), AGE = quantile(avcs$AGE, 0.75), STATUS = 4, MAL = 1)
)

scenario_grid$survival <- predict(
  fit_hm,
  newdata = scenario_grid[, c("time", "AGE", "STATUS", "MAL")],
  type = "survival"
)
scenario_grid$hazard_multiplier <- predict(
  fit_hm,
  newdata = scenario_grid[, c("time", "AGE", "STATUS", "MAL")],
  type = "hazard"
)

head(scenario_grid)
```

------------------------------------------------------------------------

## hp.dthip.PAIVS.time — Prediction with time-varying risk

**SAS source:** `examples/hp.dthip.PAIVS.time.sas`

**Purpose:** Predict patient-specific hazard over time for the PAIVS
(Pulmonary Atresia with Intact Ventricular Septum) dataset.

### R translation

``` r
set.seed(42)

# Synthetic PAIVS-like example for time-varying risk demonstration.
n <- 120
paivs <- data.frame(
     time = rexp(n, rate = 0.25) + 0.05,
     dead = rbinom(n, size = 1, prob = 0.65),
     vent_score = rnorm(n),
     shunt = rbinom(n, size = 1, prob = 0.4)
)

# Define piecewise windows for β(t): early, mid, late follow-up.
tv_windows <- c(1.0, 3.0)

# Weibull with time-varying coefficients via windowed design expansion.
# theta layout with two predictors and 3 windows:
# [mu, nu, vent_w1, shunt_w1, vent_w2, shunt_w2, vent_w3, shunt_w3]
fit_paivs <- hazard(
     time = paivs$time,
     status = paivs$dead,
     x = paivs[, c("vent_score", "shunt")],
     time_windows = tv_windows,
     theta = c(0.4, 1.1, 0, 0, 0, 0, 0, 0),
     dist = "weibull",
     fit = TRUE
)

# Time-specific risk summaries.
pred_df <- paivs
pred_df$lp <- predict(
     fit_paivs,
     newdata = pred_df[, c("time", "vent_score", "shunt")],
     type = "linear_predictor"
)
pred_df$haz_mult <- predict(
     fit_paivs,
     newdata = pred_df[, c("time", "vent_score", "shunt")],
     type = "hazard"
)

head(pred_df[, c("time", "lp", "haz_mult")])
```
