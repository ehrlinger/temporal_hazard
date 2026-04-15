## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  eval     = has_pkg
)


## -----------------------------------------------------------------------------
#| label: setup
library(TemporalHazard)
library(survival)
library(ggplot2)


## -----------------------------------------------------------------------------
#| label: kul-data
data(cabgkul)
str(cabgkul)


## -----------------------------------------------------------------------------
#| label: kul-weibull
fit_kul <- hazard(
  Surv(int_dead, dead) ~ 1,
  data  = cabgkul,
  dist  = "weibull",
  theta = c(mu = 0.10, nu = 1.0),
  fit   = TRUE
)

fit_kul


## -----------------------------------------------------------------------------
#| label: fig-kul-km
#| fig-cap: "Weibull parametric survival vs. Kaplan-Meier (CABG, KU Leuven)"
#| fig-width: 7
#| fig-height: 4.5
t_grid <- seq(0.01, max(cabgkul$int_dead) * 0.9, length.out = 200)
nd     <- data.frame(time = t_grid)
surv   <- predict(fit_kul, newdata = nd, type = "survival") * 100

km     <- survfit(Surv(int_dead, dead) ~ 1, data = cabgkul)
km_df  <- data.frame(time = km$time, survival = km$surv * 100)

ggplot() +
  geom_step(data = km_df, aes(time, survival, colour = "Kaplan-Meier"),
            linewidth = 0.5) +
  geom_line(data = data.frame(time = t_grid, survival = surv),
            aes(time, survival, colour = "Weibull"), linewidth = 1) +
  scale_colour_manual(values = c("Weibull" = "#0072B2",
                                 "Kaplan-Meier" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Months after CABG", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: avc-data
data(avc)
avc <- na.omit(avc)
str(avc)


## -----------------------------------------------------------------------------
#| label: avc-mv
fit_avc <- hazard(
  Surv(int_dead, dead) ~ age + status + mal + com_iv + inc_surg + orifice,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.20, nu = 1.0, rep(0, 6)),
  fit   = TRUE,
  control = list(maxit = 500)
)

fit_avc


## -----------------------------------------------------------------------------
#| label: avc-mp
fit_mp <- hazard(
  Surv(int_dead, dead) ~ 1,
  data   = avc,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant")
  ),
  fit     = TRUE,
  control = list(n_starts = 5, maxit = 1000)
)

summary(fit_mp)


## -----------------------------------------------------------------------------
#| label: fig-avc-compare
#| fig-cap: "Single-phase Weibull vs. multiphase model against Kaplan-Meier (AVC)"
#| fig-width: 7
#| fig-height: 4.5
t_grid <- seq(0.01, max(avc$int_dead) * 0.95, length.out = 200)
nd     <- data.frame(time = t_grid)

km_avc <- survfit(Surv(int_dead, dead) ~ 1, data = avc)
km_df  <- data.frame(time = km_avc$time, survival = km_avc$surv * 100)

fit_wb <- hazard(
  Surv(int_dead, dead) ~ 1, data = avc, dist = "weibull",
  theta = c(mu = 0.20, nu = 1.0), fit = TRUE
)

surv_wb <- predict(fit_wb, newdata = nd, type = "survival") * 100
surv_mp <- predict(fit_mp, newdata = nd, type = "survival") * 100

plot_df <- rbind(
  data.frame(time = t_grid, survival = surv_wb, Model = "Single Weibull"),
  data.frame(time = t_grid, survival = surv_mp, Model = "Multiphase (2-phase)")
)

ggplot() +
  geom_step(data = km_df, aes(time, survival), colour = "grey50",
            linewidth = 0.5) +
  geom_line(data = plot_df, aes(time, survival, colour = Model),
            linewidth = 1) +
  scale_colour_manual(values = c("Single Weibull" = "#E69F00",
                                 "Multiphase (2-phase)" = "#0072B2")) +
  scale_y_continuous(limits = c(0, 100)) +
  annotate("text", x = max(t_grid) * 0.6, y = 95, label = "KM (grey)",
           size = 3, colour = "grey50") +
  labs(x = "Months after AVC repair", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: valves-death
data(valves)
valves <- na.omit(valves)

fit_death <- hazard(
  Surv(int_dead, dead) ~ age_cop + nyha + mechvalv,
  data  = valves,
  dist  = "weibull",
  theta = c(mu = 0.10, nu = 1.0, rep(0, 3)),
  fit   = TRUE,
  control = list(maxit = 500)
)

fit_death


## -----------------------------------------------------------------------------
#| label: valves-pve
fit_pve <- hazard(
  Surv(int_pve, pve) ~ age_cop + nve + mechvalv,
  data  = valves,
  dist  = "weibull",
  theta = c(mu = 0.02, nu = 1.0, rep(0, 3)),
  fit   = TRUE,
  control = list(maxit = 500)
)

fit_pve

