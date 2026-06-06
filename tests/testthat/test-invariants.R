# Tier-1 correctness invariants -- see inst/dev/CORRECTNESS-STRATEGY.md and
# helper-invariants.R. Single-engine, dataset-driven property tests: no SAS
# reference, no captured fixtures. Run over a model matrix spanning the four
# single distributions plus multiphase (conserve on/off, covariates,
# left-truncation, weights).

test_that("predicted survival and cumulative hazard are distributionally sane", {
  skip_on_cran()
  skip_if_not_installed("survival")
  for (m in .inv_models()) {
    .inv_assert_distribution(m)
  }
})

test_that("multiphase log-likelihood is self-consistent with returned coefficients", {
  skip_on_cran()
  skip_if_not_installed("survival")
  for (m in .inv_models()) {
    if (isTRUE(m$multiphase)) .inv_assert_self_consistent(m)
  }
})

test_that("Conservation of Events identity holds for conserve = TRUE multiphase fits", {
  skip_on_cran()
  skip_if_not_installed("survival")
  for (m in .inv_models()) {
    if (isTRUE(m$multiphase) && isTRUE(m$conserve)) .inv_assert_coe(m)
  }
})
