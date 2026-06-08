# test-multiphase-hessian.R
# Analytic Hessian tests for multiphase models.
# Validates shape second-derivative helpers and .hzr_hessian_multiphase()
# against numDeriv references.

# ---------------------------------------------------------------------------
# Task 1: Shape second-derivative helpers
# ---------------------------------------------------------------------------

test_that(".hzr_phase_second_derivatives diagonals match numDeriv", {
  skip_if_not_installed("numDeriv")
  t <- c(0.2, 0.5, 1.0, 2.0, 4.0)
  t_half <- 1.0; nu <- 1.5; m <- 1.0

  d2 <- .hzr_phase_second_derivatives(t, t_half = t_half, nu = nu, m = m,
                                       type = "cdf")

  # Reference: numDeriv Hessian of Phi_i w.r.t. (log_t_half, nu, m)
  par0 <- c(log(t_half), nu, m)
  for (i in seq_along(t)) {
    phi_fn_i <- function(par) {
      hzr_phase_cumhaz(t[i], t_half = exp(par[1]), nu = par[2], m = par[3],
                       type = "cdf")
    }
    H_nd <- numDeriv::hessian(phi_fn_i, par0)
    expect_equal(d2$d2Phi_dlog_thalf2[i],       H_nd[1, 1], tolerance = 1e-3,
                 label = paste0("d2Phi/d(log_thalf)^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dnu2[i],              H_nd[2, 2], tolerance = 1e-3,
                 label = paste0("d2Phi/dnu^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dm2[i],               H_nd[3, 3], tolerance = 1e-3,
                 label = paste0("d2Phi/dm^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dlog_thalf_dnu[i],    H_nd[1, 2], tolerance = 1e-3,
                 label = paste0("d2Phi/d(log_thalf)dnu at t=", t[i]))
    expect_equal(d2$d2Phi_dlog_thalf_dm[i],     H_nd[1, 3], tolerance = 1e-3,
                 label = paste0("d2Phi/d(log_thalf)dm at t=", t[i]))
    expect_equal(d2$d2Phi_dnu_dm[i],            H_nd[2, 3], tolerance = 1e-3,
                 label = paste0("d2Phi/dnudm at t=", t[i]))
  }
})

test_that(".hzr_phase_second_derivatives phi second derivs match numDeriv", {
  skip_if_not_installed("numDeriv")
  t <- c(0.3, 0.8, 1.5, 3.0)
  t_half <- 0.5; nu <- 2.0; m <- 1.0
  d2 <- .hzr_phase_second_derivatives(t, t_half = t_half, nu = nu, m = m,
                                       type = "cdf")
  par0 <- c(log(t_half), nu, m)
  for (i in seq_along(t)) {
    haz_fn_i <- function(par) {
      hzr_phase_hazard(t[i], t_half = exp(par[1]), nu = par[2], m = par[3],
                       type = "cdf")
    }
    H_nd <- numDeriv::hessian(haz_fn_i, par0)
    expect_equal(d2$d2phi_dlog_thalf2[i], H_nd[1, 1], tolerance = 1e-3,
                 label = paste0("d2phi/d(log_thalf)^2 at t=", t[i]))
    expect_equal(d2$d2phi_dnu2[i],        H_nd[2, 2], tolerance = 1e-3,
                 label = paste0("d2phi/dnu^2 at t=", t[i]))
  }
})

test_that(".hzr_g3_phase_second_derivatives diagonals match numDeriv", {
  skip_if_not_installed("numDeriv")
  t <- c(0.5, 1.0, 3.0, 8.0)
  tau <- 2.0; gamma <- 1.5; alpha <- 0.8; eta <- 0.6

  d2 <- .hzr_g3_phase_second_derivatives(t, tau = tau, gamma = gamma,
                                          alpha = alpha, eta = eta)
  par0 <- c(log(tau), gamma, alpha, eta)
  for (i in seq_along(t)) {
    g3_fn_i <- function(par) {
      hzr_decompos_g3(t[i], tau = exp(par[1]), gamma = par[2],
                       alpha = par[3], eta = par[4])$G3
    }
    H_nd <- numDeriv::hessian(g3_fn_i, par0)
    expect_equal(d2$d2Phi_dlog_tau2[i], H_nd[1, 1], tolerance = 1e-2,
                 label = paste0("d2G3/d(log_tau)^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dgamma2[i],   H_nd[2, 2], tolerance = 1e-2,
                 label = paste0("d2G3/dgamma^2 at t=", t[i]))
    expect_equal(d2$d2Phi_dalpha2[i],   H_nd[3, 3], tolerance = 1e-2,
                 label = paste0("d2G3/dalpha^2 at t=", t[i]))
    expect_equal(d2$d2Phi_deta2[i],     H_nd[4, 4], tolerance = 1e-2,
                 label = paste0("d2G3/deta^2 at t=", t[i]))
  }
})
