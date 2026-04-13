#' @importFrom stats runif quantile rnorm rexp rbinom rlnorm approx
#' @keywords internal
NULL

# golden_fixtures.R - Synthetic golden fixture generation for parity testing
#
# PURPOSE
# -------
# Golden fixtures are pre-fitted model results saved as .rds files in
# inst/fixtures/.  They serve as reference outputs: each test run refits
# the model on the same data and compares estimates to the stored values.
# This catches regressions when the likelihood, gradient, or optimizer changes.
#
# REGENERATING FIXTURES
# ---------------------
# Fixtures must be regenerated whenever the model parameterisation changes
# in a non-backward-compatible way (e.g. a change to mu/nu encoding).
# Run each generator interactively after devtools::load_all():
#
#   .hzr_create_synthetic_golden_fixtures()   # Weibull variants
#   .hzr_create_loglogistic_golden_fixture()  # Log-logistic univariate
#   .hzr_create_lognormal_golden_fixture()    # Log-normal univariate
#
# All generators use set.seed(42) internally for reproducibility.
# Commit the updated .rds files alongside any code changes.
#
# ADDING A NEW FIXTURE
# --------------------
# 1. Write a .hzr_create_<dist>_golden_fixture() function following the
#    existing pattern (set seed, simulate, fit, saveRDS).
# 2. Call it once to write the .rds file.
# 3. Add a corresponding parity test in tests/testthat/test-parity-core.R.

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
  message("Generated fixture: hz_univariate.rds")
  
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
  message("Generated fixture: hm_multivariate.rds")
  
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
  message("Generated fixture: hm_edge_case.rds")
  
  invisible(list(
    fixture_hz_uni = fixture_hz_uni,
    fixture_hm_multi = fixture_hm_multi,
    fixture_edge = fixture_edge
  ))
}

#' Create and save log-logistic golden fixture
#'
#' Generates a synthetic dataset and fits a log-logistic model,
#' storing results as a golden fixture for validation.
#'
#' @keywords internal

.hzr_create_loglogistic_golden_fixture <- function(output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- system.file("fixtures", package = "TemporalHazard")
    if (!nzchar(output_dir)) {
      output_dir <- "inst/fixtures"
    }
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Set seed for reproducibility
  set.seed(42)
  
  # Generate univariable log-logistic sample
  n <- 80
  alpha_true <- 1.0
  beta_true <- 1.3
  
  # Generate event times from log-logistic: if U ~ Uniform(0,1),
  # then T = (U / (1-U))^(1/beta) / alpha ~ LogLogistic(alpha, beta)
  u <- runif(n)
  time <- (u / (1 - u)) ^ (1 / beta_true) / alpha_true
  
  # Add random censoring (35% censoring rate)
  cens_time <- runif(n, 0, quantile(time, 0.75))
  status <- as.integer(time <= cens_time)
  time <- pmin(time, cens_time)
  
  # Fit log-logistic model
  fit_hz_ll <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(log_alpha = 0.0, log_beta = 0.2),
    dist = "loglogistic",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  # Extract fixture
  fixture_hz_ll <- list(
    description = "Univariable log-logistic hazard: shape estimation",
    data = list(time = time, status = status, n = n, events = sum(status)),
    seed = 42,
    true_params = list(alpha = alpha_true, beta = beta_true),
    fit = if (!is.null(fit_hz_ll$fit)) {
      list(
        theta = fit_hz_ll$fit$theta,
        logl = fit_hz_ll$fit$objective,
        converged = fit_hz_ll$fit$converged == 0,
        vcov = fit_hz_ll$fit$vcov
      )
    } else NULL,
    timestamp = Sys.time()
  )
  
  # Save fixture
  fixture_file <- file.path(output_dir, "hz_loglogistic.rds")
  saveRDS(fixture_hz_ll, fixture_file)
  message("Generated fixture: hz_loglogistic.rds")
  
  invisible(fixture_hz_ll)
}

#' Create and save log-normal golden fixture
#'
#' @keywords internal

.hzr_create_lognormal_golden_fixture <- function(output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- system.file("fixtures", package = "TemporalHazard")
    if (!nzchar(output_dir)) {
      output_dir <- "inst/fixtures"
    }
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  set.seed(42)
  
  n <- 80
  mu_true <- 1.0
  sigma_true <- 0.8
  
  # Generate log-normal times: log(T) ~ Normal(mu, sigma)
  time <- rlnorm(n, meanlog = mu_true, sdlog = sigma_true)
  
  # Add random censoring (~30% censoring rate)
  cens_time <- runif(n, 0, quantile(time, 0.80))
  status <- as.integer(time <= cens_time)
  time <- pmin(time, cens_time)
  
  # Fit log-normal model (theta: mu, log(sigma))
  fit_hz_ln <- hazard(
    time = time,
    status = status,
    x = NULL,
    theta = c(mu = 1.0, log_sigma = 0.0),
    dist = "lognormal",
    fit = TRUE,
    control = list(maxit = 500, reltol = 1e-6)
  )
  
  fixture_hz_ln <- list(
    description = "Univariable log-normal hazard: location + scale estimation",
    data = list(time = time, status = status, n = n, events = sum(status)),
    seed = 42,
    true_params = list(mu = mu_true, sigma = sigma_true),
    fit = if (!is.null(fit_hz_ln$fit)) {
      list(
        theta = fit_hz_ln$fit$theta,
        logl = fit_hz_ln$fit$objective,
        converged = fit_hz_ln$fit$converged == 0,
        vcov = fit_hz_ln$fit$vcov
      )
    } else NULL,
    timestamp = Sys.time()
  )
  
  fixture_file <- file.path(output_dir, "hz_lognormal.rds")
  saveRDS(fixture_hz_ln, fixture_file)
  message("Generated fixture: hz_lognormal.rds")
  
  invisible(fixture_hz_ln)
}


# ============================================================================
# Multiphase golden fixtures
# ============================================================================

#' Create golden fixture for multiphase parity testing (synthetic)
#'
#' Generates a synthetic 3-phase dataset with known parameters, fits the
#' multiphase model, and saves the result.  Used for round-trip regression
#' testing of the optimizer.
#'
#' @keywords internal

.hzr_create_multiphase_golden_fixture <- function(output_dir = NULL) {

  if (is.null(output_dir)) {
    output_dir <- system.file("fixtures", package = "TemporalHazard")
    if (!nzchar(output_dir)) {
      output_dir <- "inst/fixtures"
    }
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  set.seed(42)

  # --- Simulate 3-phase survival data -----------------------------------------
  # True model: early CDF phase + constant background + late hazard phase
  n <- 500
  t_grid <- seq(0.01, 20, length.out = n)

  # True parameters (user scale)
  mu_early <- 0.05
  mu_const <- 0.002
  mu_late  <- 1e-6

  # Shape parameters for early and late phases
  t_half_early <- 0.5; nu_early <- 1; m_early <- 1
  t_half_late  <- 3;   nu_late  <- 1; m_late  <- 1

  # Generate cumulative hazard at each time for Poisson-thinning
  # H(t) = mu_e * G(t) + mu_c * t + mu_l * (-log(1-G_l(t)))
  d_early <- hzr_decompos(t_grid, t_half = t_half_early, nu = nu_early, m = m_early)
  d_late  <- hzr_decompos(t_grid, t_half = t_half_late,  nu = nu_late,  m = m_late)

  H_true <- mu_early * d_early$G + mu_const * t_grid +
    mu_late * (-log(pmax(1 - d_late$G, .Machine$double.xmin)))
  S_true <- exp(-H_true)

  # Inverse-CDF sampling: U ~ Uniform(0,1), find t such that S(t) = U
  u <- runif(n)
  time <- approx(S_true, t_grid, xout = u, rule = 2)$y

  # Add censoring (~30%)
  cens_time <- runif(n, 0, quantile(time, 0.80))
  status <- as.integer(time <= cens_time)
  time <- pmin(time, cens_time)

  # Remove any zero or negative times

  keep <- time > 0
  time <- time[keep]; status <- status[keep]

  # --- Fit 3-phase model with R -----------------------------------------------
  fit <- hazard(
    time = time, status = status,
    dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf",      t_half = t_half_early, nu = nu_early, m = m_early),
      constant = hzr_phase("constant"),
      late     = hzr_phase("hazard",   t_half = t_half_late,  nu = nu_late,  m = m_late)
    ),
    fit = TRUE,
    control = list(n_starts = 5, maxit = 1000, reltol = 1e-8)
  )

  fixture <- list(
    description = "Synthetic 3-phase multiphase model (early CDF + constant + late hazard)",
    data = list(time = time, status = status, n = length(time), events = sum(status)),
    seed = 42,
    true_params = list(
      mu_early = mu_early, t_half_early = t_half_early, nu_early = nu_early, m_early = m_early,
      mu_const = mu_const,
      mu_late  = mu_late,  t_half_late  = t_half_late,  nu_late  = nu_late,  m_late  = m_late
    ),
    phases = list(
      early    = list(type = "cdf",      t_half = t_half_early, nu = nu_early, m = m_early),
      constant = list(type = "constant"),
      late     = list(type = "hazard",   t_half = t_half_late,  nu = nu_late,  m = m_late)
    ),
    fit = if (!is.null(fit$fit)) {
      list(
        theta = fit$fit$theta,
        logl  = fit$fit$objective,
        converged = fit$fit$converged,
        vcov  = fit$fit$vcov
      )
    } else NULL,
    timestamp = Sys.time()
  )

  fixture_file <- file.path(output_dir, "mp_synthetic_3phase.rds")
  saveRDS(fixture, fixture_file)
  message("Generated fixture: mp_synthetic_3phase.rds")

  invisible(fixture)
}


#' Create golden fixture for C binary parity (KUL CABG dataset)
#'
#' Loads the KUL cardiac surgery dataset from inst/extdata/cabgkul.csv and
#' stores it alongside the C HAZARD binary reference output for parity testing.
#' This fixture does NOT fit a model — it stores the data and C reference
#' values so that tests can evaluate the R log-likelihood at C-estimated
#' parameters and/or fit the model and compare.
#'
#' @keywords internal

.hzr_create_c_reference_kul_fixture <- function(output_dir = NULL) {

  if (is.null(output_dir)) {
    output_dir <- system.file("fixtures", package = "TemporalHazard")
    if (!nzchar(output_dir)) {
      output_dir <- "inst/fixtures"
    }
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load KUL dataset (lazy-loaded with the package)
  data("cabgkul", package = "TemporalHazard", envir = environment())
  dat <- cabgkul
  stopifnot(nrow(dat) == 5880, sum(dat$dead) == 545)

  # C binary reference output from hz.deadp.KUL.lst
  # All shapes FIXED; only 3 log(mu) parameters estimated.
  c_reference <- list(
    description = paste(
      "C HAZARD binary output for KUL Primary Isolated CABG,",
      "3-phase (Early/Constant/Late), shapes fixed, intercept-only."
    ),
    n = 5880L,
    events = 545L,
    log_likelihood = -3740.52,
    iterations = 9L,

    # Estimated parameters (internal/log scale)
    theta = c(E0 = -3.77955, C0 = -7.2258, L0 = -16.6578),

    # Standard errors
    std_err = c(E0 = 0.09381214, C0 = 0.09312647, L0 = 0.1157703),

    # User-scale mu estimates
    mu = c(MUE = 0.02283304, MUC = 0.0007275697, MUL = 5.829177e-08),

    # Fixed shape parameters
    early_shapes = list(delta = 0, t_half = 0.2, nu = 1, m = 1),
    late_shapes  = list(tau = 1, gamma = 3, alpha = 1, eta = 1),

    # Asymptotic variance-covariance matrix
    vcov = matrix(
      c(0.008800718, -0.00116984,  0.0007409538,
       -0.00116984,   0.008672539, -0.00564017,
        0.0007409538, -0.00564017,  0.01340277),
      nrow = 3, byrow = TRUE,
      dimnames = list(c("E0", "C0", "L0"), c("E0", "C0", "L0"))
    )
  )

  fixture <- list(
    description = "KUL CABG dataset with C HAZARD binary 3-phase reference output",
    data = list(
      time   = dat$int_dead,
      status = dat$dead,
      n      = nrow(dat),
      events = sum(dat$dead)
    ),
    c_reference = c_reference,
    timestamp = Sys.time()
  )

  fixture_file <- file.path(output_dir, "mp_c_reference_kul.rds")
  saveRDS(fixture, fixture_file)
  message("Generated fixture: mp_c_reference_kul.rds")

  invisible(fixture)
}
