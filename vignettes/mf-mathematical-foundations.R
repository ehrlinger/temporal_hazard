## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = has_pkg,
  fig.width = 7,
  fig.height = 4.5
)


## -----------------------------------------------------------------------------
#| label: setup
library(TemporalHazard)


## -----------------------------------------------------------------------------
#| label: basic-decompos
t_grid <- seq(0.01, 10, length.out = 200)

# Standard sigmoidal (Case 1: m > 0, nu > 0)
d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)

# Verify the half-life property: G(t_half) = 0.5
d_half <- hzr_decompos(3, t_half = 3, nu = 2, m = 1)
cat("G(t_half) =", round(d_half$G, 6), "\n")

# Verify the identity: h(t) = g(t) / (1 - G(t))
h_check <- d$g / (1 - d$G)
cat("Max |h - g/(1-G)| =", max(abs(d$h - h_check)), "\n")


## -----------------------------------------------------------------------------
#| label: six-cases-table
#| echo: false
cases <- data.frame(
  Case = c("1", "1L", "2", "2L", "3", "3L"),
  m     = c("> 0", "= 0", "< 0", "< 0", "> 0", "= 0"),
  nu    = c("> 0", "> 0", "> 0", "= 0", "< 0", "< 0"),
  Behavior = c(
    "Standard sigmoidal",
    "Weibull-like (exponential limit as m -> 0)",
    "Heavy-tailed",
    "Exponential decay",
    "Bounded cumulative",
    "Bounded exponential"
  )
)
knitr::kable(cases, caption = "Six valid parameter sign combinations")


## -----------------------------------------------------------------------------
#| label: fig-six-cases
#| fig-cap: "G(t), g(t), and h(t) for all six valid decomposition cases"
#| fig-height: 8
if (has_ggplot) {
  library(ggplot2)

  params <- list(
    "Case 1: m=1, nu=2"     = list(nu = 2,   m = 1),
    "Case 1L: m=0, nu=2"    = list(nu = 2,   m = 0),
    "Case 2: m=-0.5, nu=2"  = list(nu = 2,   m = -0.5),
    "Case 2L: m=-0.5, nu=0" = list(nu = 0,   m = -0.5),
    "Case 3: m=1, nu=-0.5"  = list(nu = -0.5, m = 1),
    "Case 3L: m=0, nu=-0.5" = list(nu = -0.5, m = 0)
  )

  t_grid <- seq(0.01, 10, length.out = 200)
  rows <- list()

  for (nm in names(params)) {
    p <- params[[nm]]
    d <- hzr_decompos(t_grid, t_half = 3, nu = p$nu, m = p$m)
    rows <- c(rows, list(data.frame(
      time = rep(t_grid, 3),
      value = c(d$G, d$g, d$h),
      quantity = rep(c("G(t) — CDF", "g(t) — density", "h(t) — hazard"),
                     each = length(t_grid)),
      case = nm
    )))
  }

  df <- do.call(rbind, rows)

  # Cap hazard display at reasonable values for readability
  df$value[df$quantity == "h(t) — hazard" & df$value > 5] <- NA

  ggplot(df, aes(x = time, y = value)) +
    geom_line(linewidth = 0.6) +
    facet_grid(quantity ~ case, scales = "free_y") +
    labs(x = "Time", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(size = 7))
}


## -----------------------------------------------------------------------------
#| label: phase-construction
# Classic three-phase pattern
early <- hzr_phase("cdf",      t_half = 0.5, nu = 2, m = 0)
const <- hzr_phase("constant")
late  <- hzr_phase("hazard",   t_half = 5,   nu = 1, m = 0)

print(early)
print(const)
print(late)


## -----------------------------------------------------------------------------
#| label: fig-phase-shapes
#| fig-cap: "Phi(t) and phi(t) for each phase type"
#| fig-height: 5
if (has_ggplot) {
  t_grid <- seq(0.01, 12, length.out = 200)

  phi_early <- hzr_phase_cumhaz(t_grid, t_half = 0.5, nu = 2, m = 0,
                                 type = "cdf")
  phi_late  <- hzr_phase_cumhaz(t_grid, t_half = 5, nu = 1, m = 0,
                                 type = "hazard")
  phi_const <- hzr_phase_cumhaz(t_grid, type = "constant")

  dphi_early <- hzr_phase_hazard(t_grid, t_half = 0.5, nu = 2, m = 0,
                                  type = "cdf")
  dphi_late  <- hzr_phase_hazard(t_grid, t_half = 5, nu = 1, m = 0,
                                  type = "hazard")
  dphi_const <- hzr_phase_hazard(t_grid, type = "constant")

  df <- data.frame(
    time = rep(t_grid, 6),
    value = c(phi_early, phi_late, phi_const,
              dphi_early, dphi_late, dphi_const),
    phase = rep(rep(c("cdf (early)", "hazard (late)", "constant"),
                    each = length(t_grid)), 2),
    quantity = rep(c("Phi(t) — cumulative", "phi(t) — instantaneous"),
                   each = 3 * length(t_grid))
  )

  ggplot(df, aes(x = time, y = value, colour = phase)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
    scale_colour_manual(values = c(
      "cdf (early)" = "#0072B2",
      "hazard (late)" = "#D55E00",
      "constant" = "#009E73"
    )) +
    labs(x = "Time", y = NULL, colour = "Phase type") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}


## -----------------------------------------------------------------------------
#| label: fig-additive-demo
#| fig-cap: "Three-phase additive cumulative hazard: total = early + constant + late"
#| fig-height: 4
if (has_ggplot) {
  t_grid <- seq(0.01, 15, length.out = 300)

  mu_early <- 0.3
  mu_const <- 0.05
  mu_late  <- 0.2

  H_early <- mu_early * hzr_phase_cumhaz(t_grid, t_half = 0.5, nu = 2,
                                           m = 0, type = "cdf")
  H_const <- mu_const * hzr_phase_cumhaz(t_grid, type = "constant")
  H_late  <- mu_late  * hzr_phase_cumhaz(t_grid, t_half = 5, nu = 1,
                                           m = 0, type = "hazard")
  H_total <- H_early + H_const + H_late

  df <- data.frame(
    time = rep(t_grid, 4),
    cumhaz = c(H_early, H_const, H_late, H_total),
    component = rep(c("Early (cdf)", "Constant", "Late (hazard)", "Total"),
                    each = length(t_grid))
  )
  df$component <- factor(df$component,
    levels = c("Total", "Early (cdf)", "Constant", "Late (hazard)"))

  ggplot(df, aes(x = time, y = cumhaz, colour = component,
                 linewidth = component)) +
    geom_line() +
    scale_colour_manual(values = c(
      "Total" = "black", "Early (cdf)" = "#0072B2",
      "Constant" = "#009E73", "Late (hazard)" = "#D55E00"
    )) +
    scale_linewidth_manual(values = c(
      "Total" = 1.2, "Early (cdf)" = 0.6,
      "Constant" = 0.6, "Late (hazard)" = 0.6
    )) +
    labs(x = "Time", y = "Cumulative hazard H(t)",
         colour = NULL, linewidth = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}


## -----------------------------------------------------------------------------
#| label: log1mexp-demo
# Naive log(1 - exp(-x)) fails near zero
x_small <- 1e-15
cat("Naive:  ", log(1 - exp(-x_small)), "\n")
cat("Stable: ", hzr_log1mexp(x_small), "\n")


## -----------------------------------------------------------------------------
#| label: global-covariates
#| eval: false
# # Both phases use age and nyha from the global formula
# hazard(
#   survival::Surv(time, status) ~ age + nyha,
#   data   = dat,
#   dist   = "multiphase",
#   phases = list(
#     early = hzr_phase("cdf",    t_half = 0.5, nu = 2, m = 0),
#     late  = hzr_phase("hazard", t_half = 5,   nu = 1, m = 0)
#   ),
#   fit = TRUE
# )


## -----------------------------------------------------------------------------
#| label: phase-specific-covariates
#| eval: false
# # Early risk depends on surgical factors; late risk on chronic conditions
# hazard(
#   survival::Surv(time, status) ~ age + nyha + shock,
#   data   = dat,
#   dist   = "multiphase",
#   phases = list(
#     early = hzr_phase("cdf",    t_half = 0.5, nu = 2, m = 0,
#                       formula = ~ age + shock),
#     late  = hzr_phase("hazard", t_half = 5,   nu = 1, m = 0,
#                       formula = ~ age + nyha)
#   ),
#   fit = TRUE
# )


## -----------------------------------------------------------------------------
#| label: time-varying-demo
#| eval: false
# # Two windows: [0, 2] and (2, infinity)
# hazard(
#   survival::Surv(time, status) ~ age + nyha,
#   data         = dat,
#   theta        = c(mu = 0.25, nu = 1.1,
#                    age_w1 = 0, nyha_w1 = 0,
#                    age_w2 = 0, nyha_w2 = 0),
#   dist         = "weibull",
#   time_windows = 2,
#   fit          = TRUE
# )

