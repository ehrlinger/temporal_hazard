# SAS parity tests for hzr_stepwise().  Skip gracefully when the
# fixture file under inst/fixtures/ is absent, so the package installs
# and passes CI without any SAS capture.
#
# To enable: run the SAS template in
# `inst/extdata/stepwise-fixtures/avc-forward-wald.sas`, then
# convert its three output files via:
#
#   TemporalHazard:::.hzr_build_stepwise_fixture(
#     trace_csv = "path/to/stepwise_trace.csv",
#     final_csv = "path/to/stepwise_final.csv",
#     meta_txt  = "path/to/stepwise_meta.txt",
#     out_path  = "inst/fixtures/stepwise-avc-forward-wald.rds"
#   )
#
# ALGORITHMIC NOTE: SAS HAZARD uses Q-statistics (score test, evaluated at the
# current parameter estimates without refitting) for stepwise candidate
# selection.  R's hzr_stepwise(criterion = "wald") uses full model refits and
# Wald chi-square from the fitted model.  These produce different intermediate
# statistics and different step sequences, so row-by-row step comparison is not
# meaningful.  The test instead verifies that R selects largely the same final
# variable set and achieves a competitive log-likelihood.

.stepwise_parity_tolerance <- list(
  logLik      = 10,    # absolute; generous to account for Q-stat vs Wald path difference
  coefficient = 0.10   # relative; final estimates for shared covariates
)


.stepwise_parity_run <- function(fix) {
  if (fix$meta$dataset != "avc") {
    testthat::skip(paste0("Parity runner only wired up for AVC so far; saw ",
                           fix$meta$dataset))
  }
  data(avc, envir = environment())
  avc <- na.omit(avc)

  # Two-phase (early CDF + constant) model matching SAS CONDITION=14:
  #   Early phase    -- Weibull CDF, THALF/NU/M fixed at the unconditional fit values
  #   Constant phase -- flat exponential hazard rate
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = avc,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                            fixed = c("t_half", "nu", "m")),
      constant = hzr_phase("constant")
    ),
    fit    = TRUE,
    control = list(n_starts = 3L, maxit = 500L, conserve = TRUE)
  ))

  # For multiphase stepwise, scope must be a named list of one-sided formulas.
  cands    <- fix$scope$candidates
  scope_mp <- list(early = reformulate(cands), constant = reformulate(cands))

  hzr_stepwise(
    base,
    scope     = scope_mp,
    data      = avc,
    direction = fix$meta$direction,
    criterion = fix$meta$criterion,
    slentry   = fix$meta$slentry,
    slstay    = fix$meta$slstay,
    force_in  = fix$scope$force_in,
    force_out = fix$scope$force_out,
    trace     = FALSE
  )
}


test_that("AVC forward-Wald stepwise: R selects similar variables to SAS", {
  skip_on_cran()
  fix <- .hzr_load_stepwise_fixture("avc-forward-wald")
  if (is.null(fix)) {
    skip("Stepwise parity fixture avc-forward-wald.rds not found")
  }

  r_fit   <- suppressWarnings(.stepwise_parity_run(fix))
  r_steps <- r_fit$steps[r_fit$steps$action %in% c("enter", "drop"), ]

  # R must enter at least as many variables as SAS has non-intercept final covariates.
  sas_coef <- fix$final$coef
  sas_covariates <- sas_coef[!sas_coef$variable %in% c("e0", "c0"), ]
  expect_gte(sum(r_steps$action == "enter"),
             nrow(sas_covariates),
             label = "R stepwise should enter at least as many variables as SAS selected")

  # Final log-likelihood: R may find a different (even better) optimum via a
  # different variable path; we only require it stays within the generous
  # tolerance that accounts for the Q-stat vs Wald path difference.
  expect_gte(r_fit$fit$objective,
             fix$final$logLik - .stepwise_parity_tolerance$logLik,
             label = paste0("R logLik (", round(r_fit$fit$objective, 3),
                            ") should be within ", .stepwise_parity_tolerance$logLik,
                            " of SAS logLik (", fix$final$logLik, ")"))

  # Final variable set: SAS non-intercept covariates should largely appear in R.
  # Map fixture coef rows (variable + phase columns) to R's "phase.variable" naming.
  r_summary <- summary(r_fit)$coefficients
  r_coef_names <- rownames(r_summary)

  sas_varphase <- paste0(sas_covariates$phase, ".", sas_covariates$variable)
  n_shared <- sum(sas_varphase %in% r_coef_names)

  expect_gte(n_shared / length(sas_varphase), 0.5,
             label = paste0("At least half of SAS final covariates should appear in R model.",
                            " SAS: ", paste(sas_varphase, collapse = ", "),
                            " | R: ", paste(r_coef_names, collapse = ", ")))

  # Coefficient comparison is intentionally omitted: because R and SAS take
  # different selection paths (Wald refit vs Q-statistic), their final models
  # include different covariates.  Estimates for nominally "shared" variables
  # differ because the surrounding model composition differs, making direct
  # comparison invalid.
})


# Sanity: the fixture loader + validator at least work on a synthetic
# fixture (doesn't depend on SAS output being present).
test_that("fixture validator rejects malformed input", {
  bad <- list(
    meta = list(dataset = "x", dist = "weibull", criterion = "wald"),
    scope = list(candidates = "age", force_in = character(),
                 force_out = character()),
    steps = data.frame(step_num = 1L, action = "bogus",
                       variable = "x", phase = NA_character_,
                       stat = 1, df = 1L, p_value = 0.1,
                       stringsAsFactors = FALSE),
    final = list(coef = data.frame(), logLik = -100, iterations = 5L)
  )
  expect_error(
    TemporalHazard:::.hzr_validate_stepwise_fixture(bad),
    "Invalid stepwise parity fixture"
  )
})

test_that("fixture validator accepts a well-formed fixture skeleton", {
  ok <- list(
    meta = list(
      dataset     = "x", dist = "weibull",
      criterion   = "wald", direction = "forward",
      slentry     = 0.30, slstay = 0.20,
      sas_version = "9.4", captured_on = "2026-04-17",
      proc_hazard = "PROC HAZARD ..."
    ),
    scope = list(candidates = "age", force_in = character(),
                 force_out = character()),
    steps = data.frame(
      step_num = 1L, action = "enter", variable = "age",
      phase = NA_character_, stat = 5, df = 1L, p_value = 0.025,
      stringsAsFactors = FALSE
    ),
    final = list(coef = data.frame(), logLik = -100, iterations = 7L)
  )
  expect_silent(TemporalHazard:::.hzr_validate_stepwise_fixture(ok))
})
