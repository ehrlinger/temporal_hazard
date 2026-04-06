#' Generate and save golden fixtures for parity testing
#'
#' Creates synthetic datasets and fits Weibull models using the R implementation,
#' storing results as golden fixtures for validation.
#'
#' @keywords internal

.hzr_create_synthetic_golden_fixtures <- function(output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- system.file("fixtures", package = "TemporalHazard")
    if (!nzchar(output_dir)) {
      output_dir <- "inst/fixtures"
    }
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Set seed for reproducibility
  set.seed(42)
  
  # ===== Fixture 1: Univariable Weibull (shape estimation only) =====
  # Simulate from Weibull with known parameters
  n <- 100
  mu_true <- 0.5
  nu_true <- 1.2
  
  # Generate event times from Weibull
  time <- (-log(runif(n)) / mu_true) ^ (1 / nu_true)
  
  # Add random censoring (30% censoring rate)
  cens_time <- runif(n, 0, quantile(time, 0.8))
  status <- as.integer(time <= cens_time)
  time <- pmin(time, cens_time)
  
  # Fit univariable model
  fit_hz_uni <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 0.3, nu = 1.0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Extract fixture: parameter estimates, log-likelihood, convergence info
  fixture_hz_uni <- list(
    description = "Univariable Weibull hazard: shape estimation only",
    data = list(time = time, status = status, n = n, events = sum(status)),
    seed = 42,
    true_params = list(mu = mu_true, nu = nu_true),
    fit = if (!is.null(fit_hz_uni$fit)) {
      list(
        theta = fit_hz_uni$fit$theta,
        logl = fit_hz_uni$fit$objective,
        converged = fit_hz_uni$fit$converged == 0,
        vcov = fit_hz_uni$fit$vcov
      )
    } else NULL,
    timestamp = Sys.time()
  )
  
  # Save fixture
  fixture_file_1 <- file.path(output_dir, "hz_univariate.rds")
  saveRDS(fixture_hz_uni, fixture_file_1)
  message("✓ Generated fixture: hz_univariate.rds")
  
  # ===== Fixture 2: Multivariable Weibull with 2 covariates =====
  # Design matrix with 2 covariates
  x <- cbind(
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  
  # True parameters
  mu_true_mv <- 0.6
  nu_true_mv <- 1.1
  beta_true <- c(-0.3, 0.5)  # Covariate effects
  
  # Generate from Weibull with covariates
  eta <- x %*% beta_true
  time_mv <- (-log(runif(n)) / (mu_true_mv * exp(eta))) ^ (1 / nu_true_mv)
  
  # Add censoring (25%)
  cens_time_mv <- runif(n, 0, quantile(time_mv, 0.75))
  status_mv <- as.integer(time_mv <= cens_time_mv)
  time_mv <- pmin(time_mv, cens_time_mv)
  
  # Fit multivariable model
  fit_hm_multi <- hazard(
    time = time_mv,
    status = status_mv,
    x = x,
    theta = c(mu = 0.4, nu = 1.0, beta_X1 = -0.2, beta_X2 = 0.3),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Extract fixture
  fixture_hm_multi <- list(
    description = "Multivariable Weibull hazard: shape + 2 covariates",
    data = list(time = time_mv, status = status_mv, x = x, n = n, events = sum(status_mv)),
    seed = 42,
    true_params = list(mu = mu_true_mv, nu = nu_true_mv, beta = beta_true),
    fit = if (!is.null(fit_hm_multi$fit)) {
      list(
        theta = fit_hm_multi$fit$theta,
        logl = fit_hm_multi$fit$objective,
        converged = fit_hm_multi$fit$converged == 0,
        vcov = fit_hm_multi$fit$vcov
      )
    } else NULL,
    timestamp = Sys.time()
  )
  
  # Save fixture
  fixture_file_2 <- file.path(output_dir, "hm_multivariate.rds")
  saveRDS(fixture_hm_multi, fixture_file_2)
  message("✓ Generated fixture: hm_multivariate.rds")
  
  # ===== Fixture 3: High-covariate small sample (edge case) =====
  n_small <- 20
  x_small <- cbind(
    X1 = rnorm(n_small),
    X2 = rnorm(n_small),
    X3 = rnorm(n_small)
  )
  
  time_small <- rexp(n_small, rate = 0.5)
  status_small <- rbinom(n_small, 1, 0.6)
  
  fit_edge <- hazard(
    time = time_small,
    status = status_small,
    x = x_small,
    theta = c(mu = 0.5, nu = 1.2, beta_X1 = 0, beta_X2 = 0, beta_X3 = 0),
    dist = "weibull",
    fit = TRUE,
    control = list(maxit = 300, reltol = 1e-5)
  )
  
  fixture_edge <- list(
    description = "Edge case: small sample (n=20) with 3 covariates",
    data = list(time = time_small, status = status_small, x = x_small, n = n_small, events = sum(status_small)),
    seed = 42,
    fit = if (!is.null(fit_edge$fit)) {
      list(
        theta = fit_edge$fit$theta,
        logl = fit_edge$fit$objective,
        converged = fit_edge$fit$converged == 0,
        vcov = fit_edge$fit$vcov
      )
    } else NULL,
    timestamp = Sys.time()
  )
  
  fixture_file_3 <- file.path(output_dir, "hm_edge_case.rds")
  saveRDS(fixture_edge, fixture_file_3)
  message("✓ Generated fixture: hm_edge_case.rds")
  
  invisible(list(
    fixture_hz_uni = fixture_hz_uni,
    fixture_hm_multi = fixture_hm_multi,
    fixture_edge = fixture_edge
  ))
}
