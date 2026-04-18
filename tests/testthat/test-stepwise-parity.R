# SAS parity tests for hzr_stepwise().  Skip gracefully when the
# fixture file under inst/fixtures/ is absent, so the package installs
# and passes CI without any SAS capture.
#
# To enable: run the SAS template in
# `inst/extdata/stepwise-fixtures/cabgkul-forward-wald.sas`, then
# convert its three output files via:
#
#   TemporalHazard:::.hzr_build_stepwise_fixture(
#     trace_csv = "path/to/stepwise_trace.csv",
#     final_csv = "path/to/stepwise_final.csv",
#     meta_txt  = "path/to/stepwise_meta.txt",
#     out_path  = "inst/fixtures/stepwise-cabgkul-forward-wald.rds"
#   )

.stepwise_parity_tolerance <- list(
  stat        = 1e-3,   # Wald chi-square, relative
  p_value     = 1e-3,   # absolute
  logLik      = 1e-3,   # absolute
  coefficient = 1e-2    # relative
)


.stepwise_parity_run <- function(fix) {
  if (fix$meta$dataset != "cabgkul") {
    testthat::skip(paste0("Parity runner only wired up for CABGKUL so far; saw ",
                           fix$meta$dataset))
  }
  data(cabgkul, envir = environment())

  # Base model: intercept-only, same distribution as SAS, fit = TRUE
  base <- hazard(
    Surv(int_dead, dead) ~ 1,
    data  = cabgkul,
    theta = c(0.5, 1.0),
    dist  = fix$meta$dist,
    fit   = TRUE
  )

  hzr_stepwise(
    base,
    scope     = fix$scope$candidates,
    data      = cabgkul,
    direction = fix$meta$direction,
    criterion = fix$meta$criterion,
    slentry   = fix$meta$slentry,
    slstay    = fix$meta$slstay,
    force_in  = fix$scope$force_in,
    force_out = fix$scope$force_out,
    trace     = FALSE
  )
}


test_that("CABGKUL forward-Wald stepwise matches SAS output", {
  fix <- .hzr_load_stepwise_fixture("cabgkul-forward-wald")
  if (is.null(fix)) {
    skip("Stepwise parity fixture cabgkul-forward-wald.rds not found")
  }

  r_fit <- .stepwise_parity_run(fix)
  sas   <- fix$steps
  r_steps <- r_fit$steps[r_fit$steps$action %in% c("enter", "drop"), ]

  # 1. Same number of accepted steps
  expect_equal(nrow(r_steps), nrow(sas))

  # 2. Same (action, variable, phase) sequence, step by step
  expect_identical(r_steps$action,   sas$action)
  expect_identical(r_steps$variable, sas$variable)
  if ("phase" %in% names(sas)) {
    expect_identical(r_steps$phase, sas$phase)
  }

  # 3. Wald chi-square agreement.  The R driver stores z or chi^2 in
  # $stat depending on df; for df = 1 we square to compare.
  r_chi <- ifelse(r_steps$df == 1L, r_steps$stat ^ 2, r_steps$stat)
  expect_equal(r_chi, sas$stat,
               tolerance = .stepwise_parity_tolerance$stat)

  # 4. p-value agreement
  expect_equal(r_steps$p_value, sas$p_value,
               tolerance = .stepwise_parity_tolerance$p_value,
               scale = 1)

  # 5. Final log-likelihood agreement
  expect_equal(r_fit$fit$objective, fix$final$logLik,
               tolerance = .stepwise_parity_tolerance$logLik,
               scale = 1)

  # 6. Final coefficient table agreement (joined on variable name)
  r_summary <- summary(r_fit)$coefficients
  r_coef <- data.frame(
    variable  = rownames(r_summary),
    estimate  = r_summary$estimate,
    std_error = r_summary$std_error,
    stringsAsFactors = FALSE
  )
  sas_coef <- fix$final$coef[, c("variable", "estimate", "std_error")]

  common <- intersect(r_coef$variable, sas_coef$variable)
  expect_true(length(common) >= 1L,
              info = "no shared variables between R and SAS coef tables")

  for (v in common) {
    r_row   <- r_coef[r_coef$variable  == v, ]
    sas_row <- sas_coef[sas_coef$variable == v, ]
    expect_equal(r_row$estimate, sas_row$estimate,
                 tolerance = .stepwise_parity_tolerance$coefficient,
                 info = paste("estimate for", v))
    expect_equal(r_row$std_error, sas_row$std_error,
                 tolerance = .stepwise_parity_tolerance$coefficient,
                 info = paste("std_error for", v))
  }
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
