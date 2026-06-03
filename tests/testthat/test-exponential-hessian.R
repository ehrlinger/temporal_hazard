test_that(".hzr_hessian_exponential matches numDeriv (no covariates)", {
  skip_if_not_installed("numDeriv")
  set.seed(2)
  n <- 300
  time <- rexp(n, rate = 0.7); status <- rbinom(n, 1, 0.7)
  theta <- c(log_rate = log(0.7))
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = NULL)
  H_an <- .hzr_hessian_exponential(theta, time, status, x = NULL)
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_exponential matches numDeriv (covariates + weights)", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  n <- 400
  x <- cbind(a = rnorm(n), b = rbinom(n, 1, 0.4))
  eta <- as.numeric(x %*% c(0.3, -0.5))
  time <- rexp(n, rate = exp(log(0.5) + eta)); status <- rbinom(n, 1, 0.8)
  w <- runif(n, 0.5, 2)
  theta <- c(log_rate = log(0.5), 0.3, -0.5)
  obj <- function(par) -.hzr_logl_exponential(par, time, status, x = x, weights = w)
  H_an <- .hzr_hessian_exponential(theta, time, status, x = x, weights = w)
  H_nd <- numDeriv::hessian(obj, theta)
  expect_equal(unname(H_an), unname(H_nd), tolerance = 1e-4)
})

test_that(".hzr_hessian_exponential returns NULL for left/interval censoring", {
  set.seed(4)
  n <- 50
  time <- rexp(n, 0.5); status <- rep(1L, n); status[1:5] <- 2L
  tl <- time; tu <- time + 0.5
  expect_null(.hzr_hessian_exponential(c(log_rate = 0), time, status,
                                       time_lower = tl, time_upper = tu, x = NULL))
  status2 <- rep(1L, n); status2[1:5] <- -1L
  expect_null(.hzr_hessian_exponential(c(log_rate = 0), time, status2,
                                       time_upper = tu, x = NULL))
})
