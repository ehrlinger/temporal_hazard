## -----------------------------------------------------------------------------
#| include: false
has_pkg <- requireNamespace("TemporalHazard", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  eval = has_pkg
)


## -----------------------------------------------------------------------------
#| label: setup
library(TemporalHazard)
if (has_ggplot2) library(ggplot2)

data(avc)
avc <- na.omit(avc)
str(avc)




## -----------------------------------------------------------------------------
#| label: weibull-fit
fit_weib <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
  data  = avc,
  dist  = "weibull",
  theta = c(mu = 0.05, nu = 0.5),
  fit   = TRUE
)
summary(fit_weib)




## -----------------------------------------------------------------------------
#| label: multiphase-fit
fit_mp <- hazard(
  survival::Surv(int_dead, dead) ~ 1,
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
#| label: screening
covariates <- c("age", "status", "mal", "com_iv", "orifice",
                "inc_surg", "opmos")

screen_results <- data.frame(
  variable = covariates,
  coef     = numeric(length(covariates)),
  p_value  = numeric(length(covariates))
)

for (i in seq_along(covariates)) {
  v <- covariates[i]
  fml <- as.formula(paste("dead ~", v))
  lg <- glm(fml, data = avc, family = binomial)
  s <- summary(lg)$coefficients
  if (nrow(s) >= 2) {
    screen_results$coef[i] <- s[2, 1]
    screen_results$p_value[i] <- s[2, 4]
  }
}

screen_results <- screen_results[order(screen_results$p_value), ]
screen_results$significant <- ifelse(screen_results$p_value < 0.05,
                                      "*", "")
screen_results




## -----------------------------------------------------------------------------
#| label: multivariable-fit
fit_mv <- hazard(
  survival::Surv(int_dead, dead) ~ age + status + mal + com_iv,
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
summary(fit_mv)

