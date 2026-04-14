# diagnostics.R — Model validation and goodness-of-fit utilities
#
# Implements SAS HAZARD macro equivalents for calibration and validation.

#' @importFrom stats predict
NULL

#' Decile-of-risk calibration
#'
#' Partition observations into groups (default 10) by predicted risk and
#' compare observed vs. expected event counts in each group.  Good
#' calibration means the two track each other across the risk spectrum.
#'
#' This implements the workflow of the SAS `deciles.hazard.sas` macro.
#' Patients are ranked by predicted cumulative hazard at a specified time
#' point, grouped into quantile bins, and each bin is tested with a
#' chi-square goodness-of-fit statistic.
#'
#' @param object A fitted `hazard` object (with `fit = TRUE`).
#' @param time Numeric scalar: the time point at which to evaluate
#'   predicted survival / cumulative hazard.  For example, `time = 12`
#'   evaluates 12-month predictions.
#' @param groups Integer: number of risk groups (default 10 for deciles).
#' @param status Optional numeric vector of event indicators (1 = event,
#'   0 = censored).
#'   If `NULL` (default), extracted from the fitted object's stored data.
#' @param event_time Optional numeric vector of observed event/censoring
#'   times.  If `NULL` (default), extracted from the fitted object.
#'
#' @return A data frame with one row per risk group and columns:
#' \describe{
#'   \item{group}{Integer group label (1 = lowest risk).}
#'   \item{n}{Number of observations in the group.}
#'   \item{events}{Observed event count.}
#'   \item{expected}{Expected event count (sum of cumulative hazard
#'     values, which equals the predicted number of events under the
#'     fitted model).}
#'   \item{observed_rate}{Observed event rate (events / n).}
#'   \item{expected_rate}{Expected event rate (expected / n).}
#'   \item{chi_sq}{Chi-square contribution: (events - expected)^2 /
#'     expected.}
#'   \item{p_value}{Two-sided p-value from the chi-square test for
#'     this group.}
#'   \item{mean_survival}{Mean predicted survival probability in the
#'     group.}
#'   \item{mean_cumhaz}{Mean predicted cumulative hazard in the group.}
#' }
#'
#' An attribute `"overall"` is attached with the overall chi-square
#' statistic, degrees of freedom, and p-value.
#'
#' @examples
#' \donttest{
#' data(avc)
#' avc <- na.omit(avc)
#' fit <- hazard(
#'   survival::Surv(int_dead, dead) ~ age + mal,
#'   data  = avc,
#'   dist  = "weibull",
#'   theta = c(mu = 0.01, nu = 0.5, beta_age = 0, beta_mal = 0),
#'   fit   = TRUE
#' )
#' cal <- hzr_deciles(fit, time = 120)
#' print(cal)
#' }
#'
#' @seealso [predict.hazard()] for the prediction types used internally.
#' @export
hzr_deciles <- function(object, time, groups = 10L,
                        status = NULL, event_time = NULL) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a fitted hazard object.", call. = FALSE)
  }

  if (is.null(object$fit$theta) ||
      (is.logical(object$fit$converged) && is.na(object$fit$converged))) {
    stop("'object' has no fitted parameters. Refit with fit = TRUE.",
         call. = FALSE)
  }
  if (!is.numeric(time) || length(time) != 1L || time <= 0) {
    stop("'time' must be a single positive number.", call. = FALSE)
  }
  groups <- as.integer(groups)
  if (groups < 2L) {
    stop("'groups' must be at least 2.", call. = FALSE)
  }


  # --- Extract observed data ------------------------------------------------
  if (is.null(status)) {
    status <- object$data$status
  }
  if (is.null(event_time)) {
    event_time <- object$data$time
  }
  n <- length(status)
  if (length(event_time) != n) {
    stop("'status' and 'event_time' must have the same length.",
         call. = FALSE)
  }

  # --- Predicted cumulative hazard and survival at the target time ----------
  # Build newdata with the target time and any covariates from the original fit
  if (!is.null(object$data$x) && ncol(object$data$x) > 0) {
    nd <- as.data.frame(object$data$x)
    nd$time <- time
  } else {
    nd <- data.frame(time = rep(time, n))
  }

  cumhaz <- predict(object, newdata = nd, type = "cumulative_hazard")
  survival <- exp(-cumhaz)

  # --- Rank into groups by predicted risk (cumulative hazard) ---------------
  # Higher cumhaz = higher risk; group 1 = lowest risk.
  # Use ntile-style assignment that handles ties gracefully (no duplicate

  # breaks).  When all predictions are identical (intercept-only model),
  # observations are distributed as evenly as possible across groups.
  ranks <- rank(cumhaz, ties.method = "first")
  group <- as.integer(cut(ranks,
                          breaks = seq(0, n, length.out = groups + 1L),
                          include.lowest = TRUE,
                          labels = seq_len(groups)))

  # --- Aggregate by group ---------------------------------------------------
  result <- data.frame(
    group = seq_len(groups),
    n = integer(groups),
    events = numeric(groups),
    expected = numeric(groups),
    observed_rate = numeric(groups),
    expected_rate = numeric(groups),
    chi_sq = numeric(groups),
    p_value = numeric(groups),
    mean_survival = numeric(groups),
    mean_cumhaz = numeric(groups)
  )

  for (g in seq_len(groups)) {
    idx <- which(group == g)
    ng <- length(idx)
    obs_events <- sum(status[idx] == 1)
    exp_events <- sum(cumhaz[idx])

    result$n[g] <- ng
    result$events[g] <- obs_events
    result$expected[g] <- exp_events
    result$observed_rate[g] <- if (ng > 0) obs_events / ng else NA_real_
    result$expected_rate[g] <- if (ng > 0) exp_events / ng else NA_real_
    result$mean_survival[g] <- mean(survival[idx])
    result$mean_cumhaz[g] <- mean(cumhaz[idx])

    # Per-group chi-square: (O - E)^2 / E
    if (exp_events > 0) {
      result$chi_sq[g] <- (obs_events - exp_events)^2 / exp_events
      # Two-sided p-value from chi-square with 1 df
      result$p_value[g] <- stats::pchisq(result$chi_sq[g], df = 1,
                                          lower.tail = FALSE)
    } else {
      result$chi_sq[g] <- NA_real_
      result$p_value[g] <- NA_real_
    }
  }

  # --- Overall chi-square ---------------------------------------------------
  valid <- !is.na(result$chi_sq)
  overall_chi_sq <- sum(result$chi_sq[valid])
  overall_df <- sum(valid) - 1L
  overall_p <- if (overall_df > 0) {
    stats::pchisq(overall_chi_sq, df = overall_df, lower.tail = FALSE)
  } else {
    NA_real_
  }

  attr(result, "overall") <- list(
    chi_sq = overall_chi_sq,
    df = overall_df,
    p_value = overall_p,
    time = time,
    groups = groups,
    total_events = sum(status == 1),
    total_expected = sum(cumhaz)
  )

  class(result) <- c("hzr_deciles", "data.frame")
  result
}

#' Print method for hzr_deciles
#'
#' @param x An `hzr_deciles` object.
#' @param digits Number of significant digits for formatting.
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_deciles <- function(x, digits = 3, ...) {
  ov <- attr(x, "overall")
  cat("Decile-of-risk calibration at time =", ov$time, "\n")
  cat(ov$groups, "groups,", ov$total_events, "observed events,",
      round(ov$total_expected, 1), "expected\n\n")

  # Format for display
  display <- x
  display$expected <- round(display$expected, digits)
  display$observed_rate <- round(display$observed_rate, digits)
  display$expected_rate <- round(display$expected_rate, digits)
  display$chi_sq <- round(display$chi_sq, digits)
  display$p_value <- round(display$p_value, digits)
  display$mean_survival <- round(display$mean_survival, digits)
  display$mean_cumhaz <- round(display$mean_cumhaz, digits)
  class(display) <- "data.frame"
  print(display, row.names = FALSE)

  cat("\nOverall: chi-sq =", round(ov$chi_sq, digits),
      "on", ov$df, "df, p =",
      format.pval(ov$p_value, digits = digits), "\n")

  invisible(x)
}
