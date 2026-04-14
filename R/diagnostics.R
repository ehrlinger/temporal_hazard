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
#' chi-square goodness-of-fit statistic. Subjects censored before the
#' requested horizon are excluded from the observed-vs-expected comparison.
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
#'   \item{events}{Observed event count by the requested horizon
#'     (event indicator = 1 and event time <= \code{time}).}
#'   \item{expected}{Expected event count by the requested horizon,
#'     computed as the sum of individual event probabilities
#'     (\eqn{1 - S(time)}).}
#'   \item{observed_rate}{Observed event rate (events / n).}
#'   \item{expected_rate}{Expected event rate (expected / n).}
#'   \item{chi_sq}{Chi-square contribution: (events - expected)^2 /
#'     expected.}
#'   \item{p_value}{Upper-tail p-value from the chi-square test for
#'     this group (1 df).}
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
  n_obs <- length(status)
  if (length(event_time) != n_obs) {
    stop("'status' and 'event_time' must have the same length.",
         call. = FALSE)
  }
  if (any(!is.finite(event_time)) || any(event_time < 0)) {
    stop("'event_time' must be finite and non-negative.", call. = FALSE)
  }
  if (any(!is.finite(status)) || any(!status %in% c(0, 1))) {
    stop("'status' must be coded as 0/1 with finite values.", call. = FALSE)
  }

  n_model <- length(object$data$time)
  if (n_obs != n_model) {
    stop("'status'/'event_time' lengths must match the fitted data length (",
         n_model, ").", call. = FALSE)
  }

  # --- Predicted cumulative hazard and survival at the target time ----------
  # Compute per-observation predictions at the calibration horizon.
  if (identical(object$spec$dist, "multiphase")) {
    phases <- object$fit$phases
    if (is.null(phases)) phases <- object$spec$phases
    cov_counts <- object$fit$covariate_counts
    x_list <- object$fit$x_list
    cumhaz <- .hzr_multiphase_cumhaz(
      rep(time, n_model), object$fit$theta, phases, cov_counts, x_list,
      per_phase = FALSE
    )
  } else if (!is.null(object$data$x) && ncol(object$data$x) > 0) {
    nd <- as.data.frame(object$data$x)
    nd$time <- time
    cumhaz <- predict(object, newdata = nd, type = "cumulative_hazard")
  } else {
    nd <- data.frame(time = rep(time, n_model))
    cumhaz <- predict(object, newdata = nd, type = "cumulative_hazard")
  }
  survival <- exp(-cumhaz)

  # Exclude subjects censored before the calibration horizon.
  include <- (status == 1) | (event_time >= time)
  n_excluded <- sum(!include)
  if (!any(include)) {
    stop("No observations are available at the requested horizon after ",
         "excluding subjects censored before 'time'.", call. = FALSE)
  }

  status <- status[include]
  event_time <- event_time[include]
  cumhaz <- cumhaz[include]
  survival <- survival[include]
  observed <- as.integer(status == 1 & event_time <= time)
  n_included <- length(observed)

  if (groups > n_included) {
    stop("'groups' must be <= number of included observations at horizon (",
         n_included, ").", call. = FALSE)
  }

  # --- Rank into groups by predicted risk (cumulative hazard) ---------------
  # Higher cumhaz = higher risk; group 1 = lowest risk.
  # Use ntile-style assignment that handles ties gracefully (no duplicate

  # breaks).  When all predictions are identical (intercept-only model),
  # observations are distributed as evenly as possible across groups.
  ranks <- rank(cumhaz, ties.method = "first")
  group <- as.integer(cut(ranks,
                          breaks = seq(0, n_included,
                                       length.out = groups + 1L),
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
    obs_events <- sum(observed[idx] == 1)
    exp_events <- sum(1 - survival[idx])

    result$n[g] <- ng
    result$events[g] <- obs_events
    result$expected[g] <- exp_events
    result$observed_rate[g] <- if (ng > 0) obs_events / ng else NA_real_
    result$expected_rate[g] <- if (ng > 0) exp_events / ng else NA_real_
    result$mean_survival[g] <- if (ng > 0) mean(survival[idx]) else NA_real_
    result$mean_cumhaz[g] <- if (ng > 0) mean(cumhaz[idx]) else NA_real_

    # Per-group chi-square: (O - E)^2 / E
    if (exp_events > 0) {
      result$chi_sq[g] <- (obs_events - exp_events)^2 / exp_events
      # Upper-tail p-value from chi-square with 1 df
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
    total_events = sum(observed == 1),
    total_expected = sum(1 - survival),
    n_included = n_included,
    n_excluded = n_excluded
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
  if (!is.null(ov$n_included) && !is.null(ov$n_excluded)) {
    cat("Included", ov$n_included, "observations (excluded", ov$n_excluded,
        "censored before horizon).\n")
  }
  cat(ov$groups, "groups,", ov$total_events, "observed events,",
      signif(ov$total_expected, digits), "expected\n\n")

  # Format for display
  display <- x
  display$expected <- signif(display$expected, digits)
  display$observed_rate <- signif(display$observed_rate, digits)
  display$expected_rate <- signif(display$expected_rate, digits)
  display$chi_sq <- signif(display$chi_sq, digits)
  display$p_value <- signif(display$p_value, digits)
  display$mean_survival <- signif(display$mean_survival, digits)
  display$mean_cumhaz <- signif(display$mean_cumhaz, digits)
  class(display) <- "data.frame"
  print(display, row.names = FALSE)

  cat("\nOverall: chi-sq =", signif(ov$chi_sq, digits),
      "on", ov$df, "df, p =",
      format.pval(ov$p_value, digits = digits), "\n")

  invisible(x)
}


# =========================================================================
# hzr_gof — Observed vs. expected goodness-of-fit
# =========================================================================

#' Goodness-of-fit: observed vs. predicted events
#'
#' Compare a fitted hazard model against the nonparametric Kaplan-Meier
#' estimate by computing observed and expected (parametric) event counts
#' at each distinct event time.  This is the R equivalent of the SAS
#' `hazplot.sas` macro and implements the conservation-of-events
#' diagnostic.
#'
#' At each observed event time the function computes:
#' \itemize{
#'   \item The Kaplan-Meier survival and cumulative hazard.
#'   \item The parametric survival and cumulative hazard from the fitted
#'     model (and per-phase components for multiphase models).
#'   \item Cumulative observed events vs. cumulative expected events
#'     (sum of individual cumulative hazards for those exiting the risk
#'     set at each time).
#'   \item The running residual (expected minus observed).
#' }
#'
#' Perfect model fit implies the expected and observed event counts track
#' each other (residual near zero).  This is the conservation-of-events
#' principle.
#'
#' @param object A fitted `hazard` object (with `fit = TRUE`).
#' @param time_grid Optional numeric vector of time points at which to
#'   evaluate the parametric model.
#'   If `NULL` (default), uses the sorted unique event times from the
#'   fitted data.
#'
#' @return A data frame with one row per time point and columns:
#' \describe{
#'   \item{time}{Evaluation time.}
#'   \item{n_risk}{Number at risk (Kaplan-Meier).}
#'   \item{n_event}{Number of events at this time.}
#'   \item{n_censor}{Number censored at this time.}
#'   \item{km_surv}{Kaplan-Meier survival estimate.}
#'   \item{km_cumhaz}{Kaplan-Meier cumulative hazard
#'     (\eqn{-\log(\text{km\_surv})}).}
#'   \item{par_surv}{Parametric survival from the fitted model.}
#'   \item{par_cumhaz}{Parametric cumulative hazard.}
#'   \item{cum_observed}{Cumulative observed events to this time.}
#'   \item{cum_expected}{Cumulative expected events (sum of individual
#'     cumulative hazards for observations exiting the risk set).}
#'   \item{residual}{Expected minus observed
#'     (\code{cum_expected - cum_observed}).}
#' }
#'
#' For multiphase models, additional columns are appended for each
#' phase: \code{par_cumhaz_<phase>}.
#'
#' An attribute `"summary"` is attached with scalar diagnostics:
#' total observed events, total expected events, and the final residual.
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
#' gof <- hzr_gof(fit)
#' print(gof)
#'
#' # Plot observed vs expected events
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(gof, aes(x = time)) +
#'     geom_line(aes(y = cum_observed), colour = "#D55E00") +
#'     geom_line(aes(y = cum_expected), colour = "#0072B2") +
#'     labs(x = "Time", y = "Cumulative events") +
#'     theme_minimal()
#' }
#' }
#'
#' @seealso [hzr_deciles()] for decile-of-risk calibration,
#'   [predict.hazard()] for prediction types.
#' @export
hzr_gof <- function(object, time_grid = NULL) {
  if (!inherits(object, "hazard")) {
    stop("'object' must be a fitted hazard object.", call. = FALSE)
  }
  if (is.null(object$fit$theta) ||
      (is.logical(object$fit$converged) && is.na(object$fit$converged))) {
    stop("'object' has no fitted parameters. Refit with fit = TRUE.",
         call. = FALSE)
  }

  # --- Extract observed data ------------------------------------------------
  obs_time   <- object$data$time
  obs_status <- object$data$status
  n_total    <- length(obs_time)

  if (length(obs_status) != n_total) {
    stop("Stored 'time' and 'status' vectors must have the same length.",
         call. = FALSE)
  }
  if (any(!is.finite(obs_time)) || any(obs_time < 0)) {
    stop("Stored event/censoring times must be finite and non-negative.",
         call. = FALSE)
  }
  if (any(!is.finite(obs_status)) || any(!obs_status %in% c(0, 1))) {
    stop("Stored status must be coded as 0/1 with finite values.",
         call. = FALSE)
  }

  # --- Kaplan-Meier via survival::survfit -----------------------------------
  km_fit <- survival::survfit(survival::Surv(obs_time, obs_status) ~ 1)

  # survfit output: time, n.risk, n.event, n.censor, surv
  km_times   <- km_fit$time
  km_n_risk  <- km_fit$n.risk
  km_n_event <- km_fit$n.event
  km_n_censor <- km_fit$n.censor
  km_surv    <- km_fit$surv

  # --- Decide time grid -----------------------------------------------------
  if (is.null(time_grid)) {
    time_grid <- km_times
  }

  # --- Parametric predictions at each time point ----------------------------
  is_multiphase <- (object$spec$dist == "multiphase")

  if (!is.null(object$data$x) && ncol(object$data$x) > 0) {
    # For covariate models, evaluate at covariate means (baseline patient)
    x_means <- colMeans(object$data$x)
    nd <- as.data.frame(t(x_means))
    nd <- nd[rep(1, length(time_grid)), , drop = FALSE]
    nd$time <- time_grid
  } else {
    nd <- data.frame(time = time_grid)
  }

  par_cumhaz <- predict(object, newdata = nd, type = "cumulative_hazard")

  # Phase decomposition for multiphase models
  phase_cumhaz <- NULL
  if (is_multiphase) {
    decomp <- predict(object, newdata = nd, type = "cumulative_hazard",
                      decompose = TRUE)
    # decomp is a matrix; first column is "total", rest are phase names
    phase_cols <- colnames(decomp)[colnames(decomp) != "total"]
    phase_cumhaz <- as.data.frame(decomp[, phase_cols, drop = FALSE])
  }

  par_surv <- exp(-par_cumhaz)

  # --- Interpolate KM at the time grid --------------------------------------
  # Use stepfun-style interpolation for KM (right-continuous)
  km_surv_at_grid <- stats::approx(
    x = c(0, km_times), y = c(1, km_surv),
    xout = time_grid, method = "constant", f = 0, rule = 2
  )$y
  km_cumhaz_at_grid <- -log(pmax(km_surv_at_grid, .Machine$double.xmin))

  # Interpolate n.risk, n.event, n.censor at grid times
  # For event counts, sum events at matching times; 0 otherwise
  km_n_risk_grid <- stats::approx(
    x = c(0, km_times), y = c(n_total, km_n_risk),
    xout = time_grid, method = "constant", f = 0, rule = 2
  )$y
  km_n_event_grid <- rep(0, length(time_grid))
  km_n_censor_grid <- rep(0, length(time_grid))
  for (i in seq_along(km_times)) {
    match_idx <- which(abs(time_grid - km_times[i]) < .Machine$double.eps * 100)
    if (length(match_idx) > 0) {
      km_n_event_grid[match_idx[1]] <- km_n_event[i]
      km_n_censor_grid[match_idx[1]] <- km_n_censor[i]
    }
  }

  # --- Conservation of Events accounting ------------------------------------
  # At each event time, accumulate:
  #   cum_observed: running sum of observed events
  #   cum_expected: running sum of individual cumulative hazards for
  #                 observations exiting the risk set (events + censored)
  #
  # The expected events for observations leaving at time t is:
  #   (n_event + n_censor) * parametric_cumhaz(t)
  # This is the SAS hazplot approach: total * _CUMHAZ at that interval.

  cum_observed <- cumsum(km_n_event_grid)
  interval_expected <- (km_n_event_grid + km_n_censor_grid) * par_cumhaz
  cum_expected <- cumsum(interval_expected)
  residual <- cum_expected - cum_observed

  # --- Assemble result ------------------------------------------------------
  result <- data.frame(
    time         = time_grid,
    n_risk       = km_n_risk_grid,
    n_event      = km_n_event_grid,
    n_censor     = km_n_censor_grid,
    km_surv      = km_surv_at_grid,
    km_cumhaz    = km_cumhaz_at_grid,
    par_surv     = par_surv,
    par_cumhaz   = par_cumhaz,
    cum_observed = cum_observed,
    cum_expected = cum_expected,
    residual     = residual
  )

  # Add phase columns for multiphase
  if (!is.null(phase_cumhaz)) {
    for (ph in names(phase_cumhaz)) {
      result[[paste0("par_cumhaz_", ph)]] <- phase_cumhaz[[ph]]
    }
  }

  # Summary diagnostics
  attr(result, "summary") <- list(
    total_observed = cum_observed[length(cum_observed)],
    total_expected = cum_expected[length(cum_expected)],
    final_residual = residual[length(residual)],
    dist = object$spec$dist,
    n = n_total
  )

  class(result) <- c("hzr_gof", "data.frame")
  result
}

#' Print method for hzr_gof
#'
#' @param x An `hzr_gof` object.
#' @param digits Number of significant digits for formatting.
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_gof <- function(x, digits = 3, ...) {
  s <- attr(x, "summary")
  cat("Goodness-of-fit: observed vs. expected events\n")
  cat("Distribution:", s$dist, " | n =", s$n, "\n\n")
  cat("Total observed events:", s$total_observed, "\n")
  cat("Total expected events:", round(s$total_expected, digits), "\n")
  cat("Final residual (E - O):", round(s$final_residual, digits), "\n")
  cat("Conservation ratio (E/O):",
      round(s$total_expected / max(s$total_observed, 1), digits), "\n")
  cat("\nUse plot columns: time, km_surv, par_surv, cum_observed,",
      "cum_expected, residual\n")
  invisible(x)
}


# =========================================================================
# hzr_kaplan — Kaplan-Meier with exact logit confidence limits
# =========================================================================

#' Kaplan-Meier survival with exact logit confidence limits
#'
#' Compute the product-limit (Kaplan-Meier) survival estimate with
#' logit-transformed confidence limits that respect the \eqn{[0, 1]}
#' boundary.
#' This is the R equivalent of the SAS `kaplan.sas` macro.
#'
#' The standard Greenwood confidence interval can exceed \eqn{[0, 1]} in the
#' tails. The logit-transformed interval avoids this by working on the
#' log-odds scale:
#'
#' \deqn{
#'   \text{CL}_{\text{lower}} = S / \bigl(S + (1-S)\,
#'   \exp(z_\alpha\,\text{SI})\bigr)
#' }
#' \deqn{
#'   \text{CL}_{\text{upper}} = S / \bigl(S + (1-S)\,
#'   \exp(-z_\alpha\,\text{SI})\bigr)
#' }
#'
#' where \eqn{\text{SI} = \sqrt{V_P - 1} / (1 - S)}, \eqn{V_P} is the
#' cumulative Greenwood variance product, and \eqn{z_\alpha} is the
#' normal quantile for the requested confidence level.
#'
#' @param time Numeric vector of follow-up times.
#' @param status Numeric event indicator (1 = event, 0 = censored).
#' @param conf_level Confidence level for the interval (default 0.95).
#'   The SAS default of 0.68268948 corresponds to a 1-SD interval.
#' @param event_only Logical; if `TRUE` (default), only return rows at
#'   event times.
#'   If `FALSE`, include rows at censoring times as well.
#'
#' @return A data frame with one row per time point and columns:
#' \describe{
#'   \item{time}{Event/censoring time.}
#'   \item{n_risk}{Number at risk at start of interval.}
#'   \item{n_event}{Number of events at this time.}
#'   \item{n_censor}{Number censored at this time.}
#'   \item{survival}{Kaplan-Meier survival estimate.}
#'   \item{std_err}{Standard error of survival (Greenwood).}
#'   \item{cl_lower}{Lower confidence limit (logit-transformed).}
#'   \item{cl_upper}{Upper confidence limit (logit-transformed).}
#'   \item{cumhaz}{Cumulative hazard \eqn{= -\log(S)}.}
#'   \item{hazard}{Interval hazard rate
#'     \eqn{= \log(S_{t-1} / S_t) / \Delta t}.}
#'   \item{density}{Probability density estimate
#'     \eqn{= (S_{t-1} - S_t) / \Delta t}.}
#'   \item{life}{Restricted mean survival time (area under curve to
#'     this time).}
#' }
#'
#' @examples
#' data(cabgkul)
#' km <- hzr_kaplan(cabgkul$int_dead, cabgkul$dead)
#' head(km)
#'
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(km, aes(time)) +
#'     geom_step(aes(y = survival * 100)) +
#'     geom_ribbon(aes(ymin = cl_lower * 100, ymax = cl_upper * 100),
#'                 stat = "identity", alpha = 0.2) +
#'     labs(x = "Months", y = "Survival (%)") +
#'     theme_minimal()
#' }
#' }
#'
#' @seealso [hzr_gof()] for parametric vs. nonparametric comparison.
#' @export
hzr_kaplan <- function(time, status, conf_level = 0.95,
                       event_only = TRUE) {
  if (!is.numeric(time) || !is.numeric(status)) {
    stop("'time' and 'status' must be numeric vectors.", call. = FALSE)
  }
  if (length(time) != length(status)) {
    stop("'time' and 'status' must have the same length.", call. = FALSE)
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.", call. = FALSE)
  }

  z_alpha <- stats::qnorm(0.5 + 0.5 * conf_level)

  # Use survival::survfit for the core KM computation
  km_fit <- survival::survfit(survival::Surv(time, status) ~ 1)

  n_times <- length(km_fit$time)

  # --- Build output vectors at each event/censoring time --------------------
  km_time    <- km_fit$time
  km_n_risk  <- km_fit$n.risk
  km_n_event <- km_fit$n.event
  km_n_censor <- km_fit$n.censor
  km_surv    <- km_fit$surv

  # --- Greenwood variance product and exact logit CL -----------------------
  # VAR_PROD = prod_i [1/(n_i - d_i) - 1/n_i + 1] for all event times
  # The product is cumulative: each time contributes a multiplicative factor.
  var_prod <- rep(1.0, n_times)
  for (i in seq_len(n_times)) {
    ni <- km_n_risk[i]
    di <- km_n_event[i]
    if (ni > 0 && (ni - di) > 0 && di > 0) {
      factor_i <- 1.0 / (ni - di) - 1.0 / ni + 1.0
      if (i == 1L) {
        var_prod[i] <- factor_i
      } else {
        var_prod[i] <- var_prod[i - 1L] * factor_i
      }
    } else if (i > 1L) {
      var_prod[i] <- var_prod[i - 1L]
    }
  }

  # Standard error (Greenwood)
  var_exact <- km_surv^2 * pmax(var_prod - 1, 0)
  std_err <- sqrt(var_exact)

  # Logit-transformed confidence limits (SAS SI_EXACT formula)
  si_exact <- rep(0.0, n_times)
  idx_valid <- km_surv < 1 & km_surv > 0
  si_exact[idx_valid] <- sqrt(pmax(var_prod[idx_valid] - 1, 0)) /
    (1 - km_surv[idx_valid])

  cl_lower <- km_surv / (km_surv + (1 - km_surv) *
                            exp(z_alpha * si_exact))
  cl_upper <- km_surv / (km_surv + (1 - km_surv) *
                            exp(-z_alpha * si_exact))

  # Edge cases: if S = 1 or S = 0, CL = S

  cl_lower[km_surv >= 1] <- 1
  cl_upper[km_surv >= 1] <- 1
  cl_lower[km_surv <= 0] <- 0
  cl_upper[km_surv <= 0] <- 0

  # --- Cumulative hazard -----------------------------------------------------
  cumhaz <- -log(pmax(km_surv, .Machine$double.xmin))

  # --- Interval hazard and density ------------------------------------------
  lag_surv <- c(1, km_surv[-n_times])
  lag_time <- c(0, km_time[-n_times])
  delta_t <- km_time - lag_time

  hazard <- rep(NA_real_, n_times)
  density <- rep(NA_real_, n_times)
  idx_haz <- km_n_event > 0 & delta_t > 0 & km_surv > 0
  hazard[idx_haz] <- log(lag_surv[idx_haz] / km_surv[idx_haz]) /
    delta_t[idx_haz]
  density[idx_haz] <- (lag_surv[idx_haz] - km_surv[idx_haz]) /
    delta_t[idx_haz]

  # --- Restricted mean survival (life integral) -----------------------------
  # Modified trapezoidal rule: LIFE(t) = LIFE(t-1) + dt*(3S(t) - S(t-1))/2
  life <- rep(0, n_times)
  lag_life <- 0
  for (i in seq_len(n_times)) {
    life[i] <- lag_life + delta_t[i] * (3 * km_surv[i] - lag_surv[i]) / 2
    lag_life <- life[i]
  }

  # --- Assemble result ------------------------------------------------------
  result <- data.frame(
    time      = km_time,
    n_risk    = km_n_risk,
    n_event   = km_n_event,
    n_censor  = km_n_censor,
    survival  = km_surv,
    std_err   = std_err,
    cl_lower  = cl_lower,
    cl_upper  = cl_upper,
    cumhaz    = cumhaz,
    hazard    = hazard,
    density   = density,
    life      = life
  )

  if (event_only) {
    result <- result[result$n_event > 0, , drop = FALSE]
    rownames(result) <- NULL
  }

  class(result) <- c("hzr_kaplan", "data.frame")
  result
}

#' Print method for hzr_kaplan
#'
#' @param x An `hzr_kaplan` object.
#' @param digits Number of significant digits.
#' @param n Maximum rows to print (default 20).
#' @param ... Additional arguments (ignored).
#' @export
print.hzr_kaplan <- function(x, digits = 4, n = 20, ...) {
  cat("Kaplan-Meier estimate with logit confidence limits\n")
  total_events <- sum(x$n_event)
  cat("Events:", total_events, " | Time points:", nrow(x), "\n")
  if (nrow(x) > 0) {
    cat("Survival range:",
        round(min(x$survival), digits), "to",
        round(max(x$survival), digits), "\n")
    cat("RMST at last event:", round(x$life[nrow(x)], digits), "\n\n")
  }
  display <- x
  display$survival <- round(display$survival, digits)
  display$std_err <- round(display$std_err, digits)
  display$cl_lower <- round(display$cl_lower, digits)
  display$cl_upper <- round(display$cl_upper, digits)
  display$cumhaz <- round(display$cumhaz, digits)
  display$hazard <- round(display$hazard, digits)
  display$density <- round(display$density, digits)
  display$life <- round(display$life, digits)
  class(display) <- "data.frame"
  if (nrow(display) > n) {
    print(utils::head(display, n), row.names = FALSE)
    cat("... (", nrow(display) - n, " more rows)\n")
  } else {
    print(display, row.names = FALSE)
  }
  invisible(x)
}
