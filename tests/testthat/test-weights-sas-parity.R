# SAS parity test for FRACTIONAL (non-integer) observation weights.
# Roadmap 7a / FIXTURE-GAP-LIST B5: confirm SAS PROC HAZARD's WEIGHT
# statement matches hazard(weights = ...) for non-integer weights.
#
# R-side correctness is already proven without SAS by the additive-split
# and linear-scaling invariants in test-weights.R.  This test adds the
# external SAS confirmation, and like the other parity tests it skips
# gracefully when the capture fixture is absent, so the package installs
# and passes CI with no SAS output.
#
# To enable: run the SAS template in
#   inst/extdata/weights-fixtures/payload/avc-weighted.sas  (via run.sh)
# then convert the listing with
#   source("inst/extdata/weights-fixtures/parse-weighted-lst.R")
# which writes inst/fixtures/weighted-avc-fractional.rds.

.weighted_parity_tolerance <- list(
  logLik      = 0.5,    # absolute
  coefficient = 0.02    # relative; covariate betas for shared phases
)

# Re-fit the exact SAS specification in R: two-phase early-CDF + constant,
# AGE in both phases, fractional IPW weights, conservation of events on.
# Uses the identical committed avc-weighted.csv that the SAS run imports.
.weighted_parity_run <- function() {
  csv <- system.file("extdata", "weights-fixtures", "payload",
                     "avc-weighted.csv", package = "TemporalHazard")
  if (!nzchar(csv) || !file.exists(csv)) {
    testthat::skip("avc-weighted.csv payload not found")
  }
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  d <- d[stats::complete.cases(d[, c("int_dead", "dead", "age", "ipw")]), ]

  set.seed(20260611L)
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant", formula = ~ age)),
    weights = d$ipw, fit = TRUE,
    control = list(n_starts = 3L, maxit = 500L, conserve = TRUE)))
}


test_that("AVC fractional-weight fit: R covariate estimates match SAS", {
  skip_on_cran()
  fix <- .hzr_load_weighted_fixture("avc-fractional")
  if (is.null(fix)) {
    skip("Weighted parity fixture weighted-avc-fractional.rds not found")
  }

  r_fit  <- .weighted_parity_run()
  r_coef <- coef(r_fit)

  # SAS covariate rows (exclude the intercepts E0 / C0); map each to R's
  # "phase.variable" coefficient name (e.g. "early.age").
  sas <- fix$final$coef
  sas_cov <- sas[!sas$variable %in% c("e0", "c0"), ]

  for (i in seq_len(nrow(sas_cov))) {
    nm <- paste0(sas_cov$phase[i], ".", sas_cov$variable[i])
    expect_true(nm %in% names(r_coef),
                label = paste0("R model has coefficient ", nm))
    expect_equal(unname(r_coef[[nm]]), sas_cov$estimate[i],
                 tolerance = .weighted_parity_tolerance$coefficient,
                 label = paste0(nm, " (R ", round(r_coef[[nm]], 5),
                                " vs SAS ", sas_cov$estimate[i], ")"))
  }
})

test_that("AVC fractional-weight fit: R log-likelihood matches SAS", {
  skip_on_cran()
  fix <- .hzr_load_weighted_fixture("avc-fractional")
  if (is.null(fix)) {
    skip("Weighted parity fixture weighted-avc-fractional.rds not found")
  }

  r_fit <- .weighted_parity_run()
  expect_equal(r_fit$fit$objective, fix$final$logLik,
               tolerance = .weighted_parity_tolerance$logLik,
               label = paste0("R logLik (", round(r_fit$fit$objective, 3),
                              ") vs SAS logLik (", fix$final$logLik, ")"))
})


# Sanity: the loader/validator work on a synthetic fixture (no SAS output
# required).  Mirrors the stepwise fixture-validator tests.
test_that("weighted fixture validator rejects malformed input", {
  bad <- list(meta = list(dataset = "avc"), final = list(logLik = -100))
  expect_error(
    TemporalHazard:::.hzr_validate_weighted_fixture(bad),
    "Invalid weighted parity fixture"
  )
})

test_that("weighted fixture validator accepts a well-formed skeleton", {
  ok <- list(
    meta = list(
      dataset = "avc", dist = "multiphase", weight_var = "ipw",
      weight_kind = "fractional", sas_version = "9.4",
      captured_on = "2026-06-11", proc_hazard = "PROC HAZARD ... WEIGHT IPW"
    ),
    final = list(
      coef = data.frame(
        variable = c("e0", "age", "c0", "age"),
        phase    = c("early", "early", "constant", "constant"),
        estimate = c(-0.97, -0.008, -7.88, 0.0018),
        stringsAsFactors = FALSE
      ),
      logLik = -175.5
    )
  )
  expect_silent(TemporalHazard:::.hzr_validate_weighted_fixture(ok))
})
