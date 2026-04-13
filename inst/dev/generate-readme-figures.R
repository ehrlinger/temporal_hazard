# generate-readme-figures.R
# Run this script from the package root to regenerate man/figures/ PNGs.
# Requires: TemporalHazard (installed), ggplot2, survival

library(TemporalHazard)
library(ggplot2)
library(survival)

# -- Load the CABGKUL dataset (lazy-loaded with the package) -----------------
data(cabgkul)

# -- Fit a 3-phase multiphase model on CABG death ----------------------------
fit <- hazard(
  Surv(int_dead, dead) ~ 1,
  data   = cabgkul,
  dist   = "multiphase",
  phases = list(
    early    = hzr_phase("cdf",      t_half = 0.5, nu = 1, m = 1),
    constant = hzr_phase("constant"),
    late     = hzr_phase("cdf",      t_half = 10,  nu = 1, m = 1)
  ),
  fit     = TRUE,
  control = list(n_starts = 10, maxit = 2000, reltol = 1e-10)
)

# -- Time grid for prediction ------------------------------------------------
t_grid <- seq(0.01, max(cabgkul$int_dead) * 0.95, length.out = 200)
nd     <- data.frame(time = t_grid)

# -- Decomposed cumulative hazard → numerical hazard rate --------------------
decomp_H <- predict(fit, newdata = nd, type = "cumulative_hazard",
                    decompose = TRUE)

# Numerical differentiation: h(t) ≈ ΔH(t) / Δt
num_hazard <- function(cumhaz, time) {
  dt <- diff(time)
  dH <- diff(cumhaz)
  c(dH[1] / dt[1], dH / dt)
}

# -- Figure 1: Additive phase hazards ----------------------------------------
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

h_long <- rbind(
  data.frame(time = t_grid, hazard = num_hazard(decomp_H$early, t_grid),
             Phase = "Early"),
  data.frame(time = t_grid, hazard = num_hazard(decomp_H$constant, t_grid),
             Phase = "Constant"),
  data.frame(time = t_grid, hazard = num_hazard(decomp_H$late, t_grid),
             Phase = "Late"),
  data.frame(time = t_grid, hazard = num_hazard(decomp_H$total, t_grid),
             Phase = "Total")
)
h_long$Phase <- factor(h_long$Phase,
                       levels = c("Total", "Early", "Constant", "Late"))

phase_colours <- c(
  "Total"    = "#222222",
  "Early"    = "#E69F00",
  "Constant" = "#56B4E9",
  "Late"     = "#CC79A7"
)

p1 <- ggplot(h_long, aes(time, hazard, colour = Phase, linetype = Phase,
                          linewidth = Phase)) +
  geom_line() +
  scale_colour_manual(values = phase_colours) +
  scale_linetype_manual(values = c(Total = "solid", Early = "dashed",
                                   Constant = "dashed", Late = "dashed")) +
  scale_linewidth_manual(values = c(Total = 1.3, Early = 0.7,
                                    Constant = 0.7, Late = 0.7)) +
  labs(x = "Months after CABG",
       y = "Hazard rate",
       title = "Additive phase decomposition of hazard (CABGKUL, n = 5,880)",
       colour = "Phase", linetype = "Phase", linewidth = "Phase") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("man/figures/readme-hazard-phases.png", p1,
       width = 7, height = 4.5, dpi = 150, bg = "white")

# -- Figure 2: Overall survival with KM overlay ------------------------------
surv_total <- predict(fit, newdata = nd, type = "survival")
surv_df    <- data.frame(time = t_grid, survival = surv_total * 100)

km      <- survfit(Surv(int_dead, dead) ~ 1, data = cabgkul)
km_df   <- data.frame(time = km$time, survival = km$surv * 100)

p2 <- ggplot() +
  geom_step(data = km_df, aes(time, survival, colour = "Kaplan-Meier"),
            linewidth = 0.6) +
  geom_line(data = surv_df, aes(time, survival, colour = "Multiphase model"),
            linewidth = 1.1) +
  scale_colour_manual(
    values = c("Multiphase model" = "#0072B2", "Kaplan-Meier" = "#D55E00")
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Months after CABG",
       y = "Freedom from death (%)",
       title = "Multiphase parametric survival vs. Kaplan-Meier (CABGKUL)",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("man/figures/readme-survival.png", p2,
       width = 7, height = 4.5, dpi = 150, bg = "white")

message("Figures written to man/figures/")
