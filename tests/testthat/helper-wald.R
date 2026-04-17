# Shared fitting helper for Wald / candidate-score / stepwise tests.
# testthat auto-sources any helper-*.R file before running test-*.R.

.fit_weibull <- function(n = 200L, betas = 0.3, seed = 42L) {
  set.seed(seed)
  X <- matrix(rnorm(n * length(betas)), ncol = length(betas))
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  lp <- as.numeric(X %*% betas)
  time <- rexp(n, rate = exp(lp))
  status <- rep(1L, n)
  hazard(
    time   = time,
    status = status,
    x      = X,
    theta  = c(0.5, 1.0, rep(0, ncol(X))),
    dist   = "weibull",
    fit    = TRUE
  )
}
