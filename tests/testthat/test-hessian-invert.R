test_that("well-conditioned PD Hessian inverts cleanly with diagnostics", {
  H <- matrix(c(4, 1, 1, 3), 2, 2)
  res <- expect_silent(.hzr_safe_solve(H))
  expect_equal(res$vcov, solve(H), tolerance = 1e-10)
  expect_true(res$pd)
  expect_true(is.finite(res$rcond) && res$rcond > 1e-3)
})

test_that("asymmetric input is symmetrized before inversion", {
  H_sym  <- matrix(c(4, 1, 1, 3), 2, 2)
  H_asym <- matrix(c(4, 0, 2, 3), 2, 2)  # (H+t(H))/2 == H_sym
  expect_equal(.hzr_safe_solve(H_asym)$vcov, solve(H_sym), tolerance = 1e-10)
})

test_that("ill-conditioned Hessian warns by name but still returns numbers", {
  H <- matrix(c(1, 1, 1, 1 + 1e-12), 2, 2)
  expect_warning(res <- .hzr_safe_solve(H), "ill-conditioned")
  expect_true(is.matrix(res$vcov))
  expect_true(is.finite(res$rcond))
})

test_that("non-PD Hessian falls back to solve and reports pd = FALSE", {
  # Indefinite (eigenvalues +/-): chol() fails, solve() succeeds.
  H <- matrix(c(1, 2, 2, 1), 2, 2)
  res <- suppressWarnings(.hzr_safe_solve(H))
  expect_false(res$pd)
  expect_true(is.matrix(res$vcov))
})

test_that("non-positive variance diagonals are NA'd with a warning", {
  # Indefinite matrix whose inverse has a negative diagonal entry.
  H <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_warning(res <- .hzr_safe_solve(H), "[Nn]on-positive variance")
  expect_true(any(is.na(diag(res$vcov))))
})

test_that("non-finite Hessian returns NA vcov with a warning", {
  H <- matrix(c(1, NA, NA, 1), 2, 2)
  expect_warning(res <- .hzr_safe_solve(H), "non-finite")
  expect_true(all(is.na(res$vcov)))
  expect_true(is.na(res$pd))
})

test_that("scalar / non-matrix input is handled as non-finite", {
  expect_warning(res <- .hzr_safe_solve(NA), "non-finite")
  expect_true(all(is.na(res$vcov)))
})

test_that("explicit tol drives the ill-conditioned threshold", {
  H <- matrix(c(4, 1, 1, 3), 2, 2)        # well-conditioned (rcond ~ 0.4)
  expect_silent(.hzr_safe_solve(H, tol = 1e-8))      # below rcond: no warning
  expect_warning(.hzr_safe_solve(H, tol = 0.5), "ill-conditioned")  # above rcond: warns
})

test_that(".hzr_optim_generic returns rcond and pd diagnostics", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  n <- 200L
  time <- rexp(n, rate = 0.5)
  status <- rep(1L, n)
  res <- .hzr_optim_generic(
    logl_fn = .hzr_logl_exponential,
    gradient_fn = .hzr_gradient_exponential,
    time = time, status = status,
    x = NULL, theta_start = c(log_rate = 0),
    control = list(maxit = 200, reltol = 1e-8, abstol = 1e-8),
    use_bounds = FALSE
  )
  expect_true(is.finite(res$rcond))
  expect_true(isTRUE(res$pd))
  expect_true(is.matrix(res$vcov) && all(is.finite(diag(res$vcov))))
})
