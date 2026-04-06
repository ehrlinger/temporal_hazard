test_that("hzr_log1pexp is numerically stable", {
  x <- c(-1000, -10, 0, 10, 1000)
  y <- hzr_log1pexp(x)

  expect_true(all(is.finite(y)))
  expect_lt(abs(y[3] - log(2)), 1e-14)
  expect_lt(abs(y[4] - log1p(exp(10))), 1e-12)
  expect_lt(abs(y[5] - 1000), 1e-10)
})

test_that("hzr_log1mexp behaves for positive inputs", {
  x <- c(1e-8, 0.1, 1, 10)
  y <- hzr_log1mexp(x)

  expect_true(all(is.finite(y)))
  expect_equal(y[3], log(1 - exp(-1)), tolerance = 1e-12)
  expect_equal(y[4], log(1 - exp(-10)), tolerance = 1e-12)
})

test_that("hzr_log1mexp returns NA for non-positive inputs", {
  y <- hzr_log1mexp(c(-1, 0, 1))
  expect_true(is.na(y[1]))
  expect_true(is.na(y[2]))
  expect_false(is.na(y[3]))
})

test_that("hzr_clamp_prob bounds values", {
  p <- c(-1, 0, 0.25, 1, 2)
  y <- hzr_clamp_prob(p, eps = 1e-6)

  expect_equal(y[1], 1e-6)
  expect_equal(y[2], 1e-6)
  expect_equal(y[3], 0.25)
  expect_equal(y[4], 1 - 1e-6)
  expect_equal(y[5], 1 - 1e-6)
})

test_that("hzr_clamp_prob validates eps", {
  expect_error(hzr_clamp_prob(0.1, eps = 0))
  expect_error(hzr_clamp_prob(0.1, eps = 0.5))
  expect_error(hzr_clamp_prob(0.1, eps = NA_real_))
})
