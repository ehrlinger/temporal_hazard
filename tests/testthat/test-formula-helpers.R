# tests/testthat/test-formula-helpers.R
# Unit tests for .hzr_parse_formula() covering the four Surv() types it
# dispatches on (right, left, interval, counting-process) and its error
# paths.  These are the only tests that exercise R/formula-helpers.R
# directly -- integration tests elsewhere exercise it implicitly via
# hazard(), so a bug in the parser would only surface as odd downstream
# behaviour without these.

parse <- function(...) TemporalHazard:::.hzr_parse_formula(...)

# ---------------------------------------------------------------------------
# Right-censored: Surv(time, status)
# ---------------------------------------------------------------------------

test_that("right-censored Surv parses into time + status only", {
  df <- data.frame(
    t = c(1.5, 2.0, 3.25),
    d = c(1, 0, 1),
    x = c(0.1, 0.2, 0.3)
  )
  out <- parse(survival::Surv(t, d) ~ x, data = df)

  expect_equal(out$time,   df$t)
  expect_equal(out$status, df$d)
  expect_null(out$time_lower)
  expect_null(out$time_upper)
  expect_equal(as.numeric(out$x), df$x)
  expect_equal(colnames(out$x), "x")
})

# ---------------------------------------------------------------------------
# Left-censored: Surv(time, status, type = "left")
# ---------------------------------------------------------------------------

test_that("left-censored Surv sets time_upper but not time_lower", {
  df <- data.frame(
    t = c(1, 2, 3),
    d = c(1, 0, 1)
  )
  out <- parse(
    survival::Surv(t, d, type = "left") ~ 1,
    data = df
  )

  expect_equal(out$time,       df$t)
  expect_equal(out$status,     df$d)
  expect_null(out$time_lower)
  expect_equal(out$time_upper, df$t)
  expect_null(out$x)
})

# ---------------------------------------------------------------------------
# Interval-censored: Surv(time1, time2, type = "interval2")
# ---------------------------------------------------------------------------

test_that("interval-censored Surv extracts time_lower + time_upper + status", {
  # All rows interval-censored so Surv doesn't overwrite any bound.
  df <- data.frame(
    lo = c(1.0, 2.0, 3.0),
    hi = c(1.5, 2.8, 3.5)
  )
  out <- parse(
    survival::Surv(lo, hi, type = "interval2") ~ 1,
    data = df
  )

  # interval status codes coming out of Surv are 0/1/2/3; the parser
  # keeps them as-is, just routes lower/upper columns.
  expect_length(out$time, 3)
  expect_equal(out$time_lower, df$lo)
  expect_equal(out$time_upper, df$hi)
  expect_true(all(out$status %in% c(0, 1, 2, 3)))
})

# ---------------------------------------------------------------------------
# Counting-process (start-stop): Surv(start, stop, event)
# ---------------------------------------------------------------------------

test_that("counting-process Surv maps start->time_lower and stop->time", {
  df <- data.frame(
    start = c(0.0, 1.5, 0.0),
    stop  = c(1.5, 2.5, 3.0),
    event = c(0,   1,   1),
    age   = c(50,  50,  62)
  )
  out <- parse(survival::Surv(start, stop, event) ~ age, data = df)

  expect_equal(out$time_lower, df$start)
  expect_equal(out$time,       df$stop)
  expect_equal(out$status,     df$event)
  expect_null(out$time_upper)
  expect_equal(as.numeric(out$x), df$age)
})

# ---------------------------------------------------------------------------
# Predictor handling
# ---------------------------------------------------------------------------

test_that("intercept column is dropped from the design matrix", {
  df <- data.frame(t = c(1, 2, 3), d = c(1, 0, 1), a = c(10, 20, 30))
  out <- parse(survival::Surv(t, d) ~ a, data = df)
  expect_false("(Intercept)" %in% colnames(out$x))
  expect_equal(colnames(out$x), "a")
})

test_that("intercept-only RHS returns NULL design matrix", {
  df <- data.frame(t = c(1, 2, 3), d = c(1, 0, 1))
  out <- parse(survival::Surv(t, d) ~ 1, data = df)
  expect_null(out$x)
})

test_that("factor predictors expand to dummy columns", {
  df <- data.frame(
    t = c(1, 2, 3, 4),
    d = c(1, 0, 1, 0),
    g = factor(c("a", "b", "a", "c"))
  )
  out <- parse(survival::Surv(t, d) ~ g, data = df)
  # Factor with 3 levels -> 2 dummies (baseline absorbed)
  expect_equal(ncol(out$x), 2L)
  expect_true(all(grepl("^g", colnames(out$x))))
})

test_that("multi-term RHS preserves all predictors in order", {
  df <- data.frame(
    t = 1:4, d = c(1, 0, 1, 0),
    a = c(0.1, 0.2, 0.3, 0.4),
    b = c(1, 2, 3, 4),
    c = c(10, 20, 30, 40)
  )
  out <- parse(survival::Surv(t, d) ~ a + b + c, data = df)
  expect_equal(colnames(out$x), c("a", "b", "c"))
})

# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------

test_that("non-data.frame data is rejected", {
  expect_error(
    parse(survival::Surv(t, d) ~ 1, data = list(t = 1:3, d = c(1, 0, 1))),
    "data frame"
  )
})

test_that("non-formula input is rejected", {
  df <- data.frame(t = 1:3, d = c(1, 0, 1))
  expect_error(parse("Surv(t, d) ~ 1", data = df), "formula")
})

test_that("LHS without Surv() is rejected", {
  df <- data.frame(t = 1:3, d = c(1, 0, 1))
  expect_error(parse(t ~ d, data = df), "Surv\\(\\) call")
})

test_that("unknown predictor in RHS raises an informative error", {
  df <- data.frame(t = 1:3, d = c(1, 0, 1))
  expect_error(
    parse(survival::Surv(t, d) ~ nonexistent, data = df),
    "Failed to parse formula RHS"
  )
})
