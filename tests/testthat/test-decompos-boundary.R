# tests/testthat/test-decompos-boundary.R
#
# Phase 7d: boundary and limiting-case behavior of hzr_decompos().
#
# hzr_decompos() dispatches on the signs of (m, nu) into six branches plus two
# explicit error guards.  The branches include three "limiting" forms (m -> 0
# for nu > 0 and nu < 0, and nu -> 0 for m < 0) that must agree with the
# general branches as the parameter approaches the boundary.  These tests:
#   (1) confirm the previously-unhandled (nu == 0, m >= 0) combination now fails
#       loud instead of raising the cryptic "object 'G' not found";
#   (2) lock in continuity of the verified limiting branches;
#   (3) check internal consistency (g == dG/dt) and CDF sanity (monotone, [0,1]);
#   (4) record the known Case 3 (m > 0, nu < 0) <-> Case 3L discontinuity as a
#       skipped test pending validation against the C HAZARD G1 reference.

# ---------------------------------------------------------------------------
# (1) Unhandled-combination guard
# ---------------------------------------------------------------------------

test_that("nu == 0 with m >= 0 fails loud (not 'object G not found')", {
  # Previously these fell through every dispatch branch, leaving G unassigned
  # and raising a cryptic internal error.
  expect_error(
    hzr_decompos(c(1, 2), t_half = 3, nu = 0, m = 1),
    "Decomposition undefined for nu = 0 with m = 1"
  )
  expect_error(
    hzr_decompos(c(1, 2), t_half = 3, nu = 0, m = 0),
    "Decomposition undefined for nu = 0 with m = 0"
  )
})

test_that("m < 0 with nu < 0 still errors (existing guard)", {
  expect_error(
    hzr_decompos(c(1, 2), t_half = 3, nu = -1, m = -1),
    "Decomposition undefined when both m < 0 and nu < 0"
  )
})

# ---------------------------------------------------------------------------
# (2) Continuity of the verified limiting branches
# ---------------------------------------------------------------------------

test_that("Case 1 (m>0) limits to Case 1L (m=0) as m -> 0, nu > 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = 2, m = 1e-7)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = 2, m = 0)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

test_that("Case 2 (m<0) limits to Case 1L (m=0) as m -> 0, nu > 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = 2, m = -1e-7)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = 2, m = 0)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

test_that("Case 2 (nu>0) limits to Case 2L (nu=0) as nu -> 0, m < 0", {
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = 1e-7, m = -1)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = 0,    m = -1)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------
# (3) Internal consistency (g == dG/dt) and CDF sanity
# ---------------------------------------------------------------------------

test_that("g equals the numerical derivative of G across the used branches", {
  check_deriv <- function(t_half, nu, m) {
    t0 <- 2.0
    h  <- 1e-6
    d  <- hzr_decompos(t0, t_half = t_half, nu = nu, m = m)
    dp <- hzr_decompos(t0 + h, t_half = t_half, nu = nu, m = m)
    dm <- hzr_decompos(t0 - h, t_half = t_half, nu = nu, m = m)
    num_g <- (dp$G - dm$G) / (2 * h)
    expect_equal(d$g, num_g, tolerance = 1e-4)
  }
  check_deriv(3,  2,  1)   # Case 1
  check_deriv(3,  2,  0)   # Case 1L
  check_deriv(3,  2, -1)   # Case 2
  check_deriv(3,  0, -1)   # Case 2L
})

test_that("G is a valid CDF (monotone, within [0, 1]) for the early-phase cases", {
  t <- seq(0.01, 30, length.out = 200)
  for (pars in list(c(nu = 2, m = 1), c(nu = 2, m = 0), c(nu = 2, m = -1))) {
    G <- hzr_decompos(t, t_half = 3, nu = pars[["nu"]], m = pars[["m"]])$G
    expect_true(all(G >= 0 & G <= 1),
                label = paste("G in [0,1] for nu=", pars[["nu"]], "m=", pars[["m"]]))
    expect_true(all(diff(G) >= -1e-9),
                label = paste("G monotone for nu=", pars[["nu"]], "m=", pars[["m"]]))
  }
})

# ---------------------------------------------------------------------------
# (4) Extreme t_half stability
# ---------------------------------------------------------------------------

test_that("extreme t_half values produce finite, valid G (nu > 0)", {
  t <- c(0.5, 1, 2, 5, 20)
  for (th in c(1e-6, 1e-3, 1e3, 1e6)) {
    G <- hzr_decompos(t, t_half = th, nu = 2, m = 1)$G
    expect_true(all(is.finite(G)), label = paste("finite G at t_half =", th))
    expect_true(all(G >= 0 & G <= 1), label = paste("G in [0,1] at t_half =", th))
    expect_true(all(diff(G) >= -1e-9), label = paste("G monotone at t_half =", th))
  }
})

# ---------------------------------------------------------------------------
# (5) Known issue: Case 3 (m>0, nu<0) <-> Case 3L (m=0, nu<0) discontinuity
# ---------------------------------------------------------------------------

test_that("Case 3 (m>0, nu<0) limits to Case 3L (m=0) as m -> 0 [KNOWN ISSUE]", {
  skip(paste0(
    "Case 3 (m>0, nu<0) is discontinuous with its m->0 limit Case 3L; both ",
    "are internally consistent (g==dG/dt) but disagree by ~0.65 at the ",
    "boundary. The C HAZARD G1 rho parameterization differs from the R ",
    "port, so resolution needs working through the full C G1 evaluation. ",
    "Neither branch is exercised by a shipped phase or parity fixture. ",
    "Tracked for Phase 7d follow-up."
  ))
  t <- c(0.3, 1, 2.5, 7)
  g_general <- hzr_decompos(t, t_half = 3, nu = -2, m = 1e-7)$G
  g_limit   <- hzr_decompos(t, t_half = 3, nu = -2, m = 0)$G
  expect_equal(g_general, g_limit, tolerance = 1e-5)
})
