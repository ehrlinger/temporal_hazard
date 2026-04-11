# test-decomposition.R — Tests for hzr_decompos() and phase helpers
#
# Validates all 6 parameter cases, mathematical identities, and parity
# with the mixhazard package's decompos() function.

t_grid <- seq(0.1, 10, by = 0.1)

# ============================================================================
# Input validation
# ============================================================================

test_that("hzr_decompos rejects invalid inputs", {
  expect_error(hzr_decompos("a", t_half = 1, nu = 1, m = 0),
               "time must be numeric")
  expect_error(hzr_decompos(1:5, t_half = -1, nu = 1, m = 0),
               "t_half must be a positive scalar")
  expect_error(hzr_decompos(1:5, t_half = 0, nu = 1, m = 0),
               "t_half must be a positive scalar")
  expect_error(hzr_decompos(1:5, t_half = 1, nu = NA, m = 0),
               "nu must be a numeric scalar")
  expect_error(hzr_decompos(1:5, t_half = 1, nu = 1, m = Inf),
               "m must be a numeric scalar")
})

test_that("hzr_decompos rejects m < 0 AND nu < 0", {
  expect_error(hzr_decompos(1:5, t_half = 1, nu = -1, m = -1),
               "both m < 0 and nu < 0")
})


# ============================================================================
# Case 1: m > 0, nu > 0 (standard sigmoidal)
# ============================================================================

test_that("Case 1 (m>0, nu>0): G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))  # monotone non-decreasing
})

test_that("Case 1: G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 2, m = 1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 1: g >= 0 and h >= 0", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 1)
  expect_true(all(d$g >= -1e-15))
  expect_true(all(d$h >= -1e-15))
})


# ============================================================================
# Case 1L: m = 0, nu > 0 (Weibull-like)
# ============================================================================

test_that("Case 1L (m=0, nu>0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 2, m = 0)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 1L: G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Case 2: m < 0, nu > 0 (heavy-tailed)
# ============================================================================

test_that("Case 2 (m<0, nu>0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 2, m = -1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 2: G is monotone increasing", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = -1)
  expect_true(all(diff(d$G) >= -1e-12))
  expect_true(all(d$g >= -1e-15))
})


# ============================================================================
# Case 2L: m < 0, nu = 0 (exponential decay)
# ============================================================================

test_that("Case 2L (m<0, nu=0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = 0, m = -1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 2L: G is monotone increasing", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 0, m = -1)
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Case 3: m > 0, nu < 0 (bounded cumulative)
# ============================================================================

test_that("Case 3 (m>0, nu<0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = -2, m = 1)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 3: G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = -2, m = 1)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Case 3L: m = 0, nu < 0 (bounded exponential)
# ============================================================================

test_that("Case 3L (m=0, nu<0): G(t_half) = 0.5", {
  d <- hzr_decompos(3, t_half = 3, nu = -2, m = 0)
  expect_equal(d$G, 0.5, tolerance = 1e-10)
})

test_that("Case 3L: G is monotone increasing in [0,1]", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = -2, m = 0)
  expect_true(all(d$G >= 0 & d$G <= 1))
  expect_true(all(diff(d$G) >= -1e-12))
})


# ============================================================================
# Mathematical identities
# ============================================================================

test_that("h(t) = g(t) / (1 - G(t)) identity holds across all cases", {
  cases <- list(
    list(t_half = 3, nu =  2, m =  1),   # Case 1
    list(t_half = 3, nu =  2, m =  0),   # Case 1L
    list(t_half = 3, nu =  2, m = -1),   # Case 2
    list(t_half = 3, nu =  0, m = -1),   # Case 2L
    list(t_half = 3, nu = -2, m =  1),   # Case 3
    list(t_half = 3, nu = -2, m =  0)    # Case 3L
  )

  for (pars in cases) {
    d <- hzr_decompos(t_grid, t_half = pars$t_half,
                      nu = pars$nu, m = pars$m)
    # Only check where G(t) < 1 - eps (hazard blows up as G -> 1)
    valid <- d$G < 0.999
    if (sum(valid) > 0) {
      h_manual <- d$g[valid] / (1 - d$G[valid])
      expect_equal(d$h[valid], h_manual, tolerance = 1e-10,
                   label = paste0("h = g/(1-G) for m=", pars$m,
                                  ", nu=", pars$nu))
    }
  }
})

test_that("G(t) integrates consistently with g(t) via numerical check", {
  # For a smooth case, check that diff(G) ~ g * dt
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  dt <- diff(t_grid)
  dG <- diff(d$G)
  g_mid <- (d$g[-1] + d$g[-length(d$g)]) / 2  # midpoint rule
  ratio <- dG / (g_mid * dt)
  # Should be approximately 1 everywhere

  expect_true(all(abs(ratio - 1) < 0.05),
              label = "dG/dt approximately equals g(t)")
})


# ============================================================================
# hzr_phase_cumhaz
# ============================================================================

test_that("hzr_phase_cumhaz('cdf') returns G(t)", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "cdf")
  expect_equal(phi, d$G)
})

test_that("hzr_phase_cumhaz('hazard') returns -log(1 - G(t))", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "hazard")
  expected <- -log(1 - d$G)
  expect_equal(phi, expected, tolerance = 1e-12)
})

test_that("hzr_phase_cumhaz('constant') returns t", {
  phi <- hzr_phase_cumhaz(t_grid, type = "constant")
  expect_equal(phi, t_grid)
})

test_that("hzr_phase_cumhaz('hazard') is monotone increasing", {
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "hazard")
  expect_true(all(diff(phi) >= -1e-12))
})

test_that("hzr_phase_cumhaz('cdf') is bounded [0, 1]", {
  phi <- hzr_phase_cumhaz(t_grid, t_half = 3, nu = 2, m = 0, type = "cdf")
  expect_true(all(phi >= 0 & phi <= 1))
})


# ============================================================================
# hzr_phase_hazard
# ============================================================================

test_that("hzr_phase_hazard('cdf') returns g(t)", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_hazard(t_grid, t_half = 3, nu = 2, m = 0, type = "cdf")
  expect_equal(phi, d$g)
})

test_that("hzr_phase_hazard('hazard') returns h(t)", {
  d <- hzr_decompos(t_grid, t_half = 3, nu = 2, m = 0)
  phi <- hzr_phase_hazard(t_grid, t_half = 3, nu = 2, m = 0, type = "hazard")
  expect_equal(phi, d$h)
})

test_that("hzr_phase_hazard('constant') returns vector of 1s", {
  phi <- hzr_phase_hazard(t_grid, type = "constant")
  expect_equal(phi, rep(1, length(t_grid)))
})


# ============================================================================
# Parity with mixhazard::decompos()
# ============================================================================

test_that("hzr_decompos matches mixhazard::decompos() across all cases", {
  skip_if_not_installed("mixhazard")

  cases <- list(
    list(t_half = 3, nu =  2, m =  1),   # Case 1
    list(t_half = 3, nu =  2, m =  0),   # Case 1L
    list(t_half = 3, nu =  2, m = -1),   # Case 2
    list(t_half = 3, nu =  0, m = -1),   # Case 2L
    list(t_half = 3, nu = -2, m =  1),   # Case 3
    list(t_half = 3, nu = -2, m =  0)    # Case 3L
  )

  for (pars in cases) {
    ours   <- hzr_decompos(t_grid, t_half = pars$t_half,
                           nu = pars$nu, m = pars$m)
    theirs <- mixhazard::decompos(t_grid, thalf = pars$t_half,
                                  nu = pars$nu, m = pars$m,
                                  complet = 1)
    expect_equal(ours$G, theirs$capgt, tolerance = 1e-12,
                 label = paste0("G parity for m=", pars$m,
                                ", nu=", pars$nu))
    expect_equal(ours$g, theirs$gt, tolerance = 1e-12,
                 label = paste0("g parity for m=", pars$m,
                                ", nu=", pars$nu))
    expect_equal(ours$h, theirs$ht, tolerance = 1e-12,
                 label = paste0("h parity for m=", pars$m,
                                ", nu=", pars$nu))
  }
})


# ============================================================================
# Numerical stability edge cases
# ============================================================================

test_that("hzr_decompos handles very small times without error", {
  d <- hzr_decompos(1e-10, t_half = 3, nu = 2, m = 0)
  expect_true(is.finite(d$G))
  expect_true(is.finite(d$g))
  expect_true(is.finite(d$h))
})

test_that("hzr_decompos handles very large times without error", {
  d <- hzr_decompos(1e6, t_half = 3, nu = 2, m = 0)
  expect_true(is.finite(d$G))
  # g and h may be 0 or Inf at extreme times, but should not be NaN
  expect_false(is.nan(d$G))
})

test_that("hzr_decompos handles time = 0 gracefully (clamped)", {
  d <- hzr_decompos(0, t_half = 3, nu = 2, m = 0)
  expect_true(is.finite(d$G))
  expect_true(d$G >= 0)
})
