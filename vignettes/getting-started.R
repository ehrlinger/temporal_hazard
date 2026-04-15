## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = has_pkg
)


## -----------------------------------------------------------------------------
library(TemporalHazard)
library(survival)

set.seed(1001)
n <- 180
dat <- data.frame(
  time = rexp(n, rate = 0.35) + 0.05,
  status = rbinom(n, size = 1, prob = 0.6),
  age = rnorm(n, mean = 62, sd = 11),
  nyha = sample(1:4, n, replace = TRUE),
  shock = rbinom(n, size = 1, prob = 0.18)
)

fit <- TemporalHazard::hazard(
  Surv(time, status) ~ age + nyha + shock,
  data = dat,
  theta = c(mu = 0.25, nu = 1.10, beta1 = 0, beta2 = 0, beta3 = 0),
  dist = "weibull",
  fit = TRUE,
  control = list(maxit = 300)
)

summary(fit)


## -----------------------------------------------------------------------------
new_patients <- data.frame(
  time = c(0.5, 1.5, 3.0),
  age = c(50, 65, 75),
  nyha = c(1, 3, 4),
  shock = c(0, 0, 1)
)

pred_input <- new_patients

new_patients$linear_predictor <- predict(fit, newdata = pred_input, type = "linear_predictor")
new_patients$hazard_multiplier <- predict(fit, newdata = pred_input, type = "hazard")
new_patients$survival <- predict(fit, newdata = pred_input, type = "survival")
new_patients$cumulative_hazard <- predict(fit, newdata = pred_input, type = "cumulative_hazard")

new_patients


## -----------------------------------------------------------------------------
#| label: fig-survival
#| fig-cap: "Parametric Weibull survival curve with Kaplan-Meier overlay"
#| fig-width: 7
#| fig-height: 5
library(ggplot2)

# Parametric curve on a fine grid
t_grid <- seq(0.05, max(dat$time), length.out = 80)
curve_df <- data.frame(
  time  = t_grid,
  age   = median(dat$age),
  nyha  = 2,
  shock = 0
)
curve_df$survival <- predict(fit, newdata = curve_df, type = "survival") * 100

# Kaplan-Meier empirical overlay
km <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
km_df <- data.frame(time = km$time, survival = km$surv * 100)

ggplot() +
  geom_step(data = km_df, aes(time, survival, colour = "Kaplan-Meier")) +
  geom_line(data = curve_df, aes(time, survival,
                                  colour = "Parametric (Weibull)")) +
  scale_colour_manual(
    values = c("Parametric (Weibull)" = "#0072B2",
               "Kaplan-Meier"         = "#D55E00")
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Months after surgery", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: multiphase-fit
# CABGKUL is the benchmark dataset for 3-phase decomposition (n = 5,880)
data(cabgkul)

fit_mp <- hazard(
  Surv(int_dead, dead) ~ 1,
  data   = cabgkul,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                          fixed = "shapes"),
    constant = hzr_phase("constant"),
    late     = hzr_phase("g3",  tau = 1, gamma = 3, alpha = 1, eta = 1,
                          fixed = "shapes")
  ),
  fit     = TRUE,
  control = list(n_starts = 5, maxit = 1000)
)

summary(fit_mp)


## -----------------------------------------------------------------------------
#| label: fig-hazard-phases
#| fig-cap: "Additive phase decomposition: total hazard (solid) = early + constant + late (dashed)"
#| fig-width: 7
#| fig-height: 4.5
t_grid <- seq(0.01, max(cabgkul$int_dead) * 0.95, length.out = 200)
nd     <- data.frame(time = t_grid)

# decompose = TRUE returns per-phase cumulative hazard columns
decomp <- predict(fit_mp, newdata = nd, type = "cumulative_hazard",
                  decompose = TRUE)

# Numerical differentiation: h(t) ≈ ΔH(t) / Δt
num_hazard <- function(cumhaz, time) {
  dt <- diff(time)
  dH <- diff(cumhaz)
  c(dH[1] / dt[1], dH / dt)
}

h_long <- rbind(
  data.frame(time = t_grid, hazard = num_hazard(decomp$early, t_grid),
             Phase = "Early"),
  data.frame(time = t_grid, hazard = num_hazard(decomp$constant, t_grid),
             Phase = "Constant"),
  data.frame(time = t_grid, hazard = num_hazard(decomp$late, t_grid),
             Phase = "Late"),
  data.frame(time = t_grid, hazard = num_hazard(decomp$total, t_grid),
             Phase = "Total")
)
h_long$Phase <- factor(h_long$Phase,
                       levels = c("Total", "Early", "Constant", "Late"))

ggplot(h_long, aes(time, hazard, colour = Phase, linetype = Phase)) +
  geom_line(aes(linewidth = Phase)) +
  scale_colour_manual(values = c(Total = "#222222", Early = "#E69F00",
                                 Constant = "#56B4E9", Late = "#CC79A7")) +
  scale_linetype_manual(values = c(Total = "solid", Early = "dashed",
                                   Constant = "dashed", Late = "dashed")) +
  scale_linewidth_manual(values = c(Total = 1.3, Early = 0.7,
                                    Constant = 0.7, Late = 0.7)) +
  labs(x = "Months after CABG", y = "Hazard rate",
       colour = "Phase", linetype = "Phase", linewidth = "Phase") +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| label: fig-mp-survival
#| fig-cap: "Multiphase parametric survival vs. Kaplan-Meier"
#| fig-width: 7
#| fig-height: 4.5
surv_df <- data.frame(
  time     = t_grid,
  survival = predict(fit_mp, newdata = nd, type = "survival") * 100
)

km    <- survfit(Surv(int_dead, dead) ~ 1, data = cabgkul)
km_df <- data.frame(time = km$time, survival = km$surv * 100)

ggplot() +
  geom_step(data = km_df, aes(time, survival, colour = "Kaplan-Meier"),
            linewidth = 0.6) +
  geom_line(data = surv_df, aes(time, survival, colour = "Multiphase model"),
            linewidth = 1.1) +
  scale_colour_manual(
    values = c("Multiphase model" = "#0072B2", "Kaplan-Meier" = "#D55E00")
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Months after CABG", y = "Freedom from death (%)",
       colour = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
#| eval: false
# hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 1)   # Early risk (bounded)
# hzr_phase("constant")                                   # Flat background rate
# hzr_phase("cdf",      t_half = 10,  nu = 1, m = 0)    # Late risk (CDF-based)
# hzr_phase("g3",       tau = 1, gamma = 3, alpha = 1,   # Late risk (G3 power law)
#                        eta = 1)


## -----------------------------------------------------------------------------
TemporalHazard::hzr_log1pexp(c(-2, 0, 2))

