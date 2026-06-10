#' @keywords internal
NULL

# hessian-invert.R -- Numerically stable Hessian inversion with diagnostics.

# Reciprocal-condition warning threshold (~1.5e-8). Shared by .hzr_safe_solve()
# and the summary() diagnostic note so they never drift apart.
.hzr_rcond_tol <- .Machine$double.eps^0.5

#' Stable Hessian inversion with conditioning diagnostics
#'
#' Inverts a negative-log-likelihood Hessian into a variance-covariance
#' matrix, hardening against the ill-conditioning that arises at high
#' parameter counts.  Symmetrizes the input, checks the reciprocal condition
#' number, inverts via Cholesky (with a \code{solve()} fallback for non-PD
#' Hessians), and guards non-positive variances.  Emits a named warning for
#' every degenerate path.
#'
#' @param H Square numeric Hessian of the negative log-likelihood.
#' @param tol Reciprocal-condition warning threshold.
#' @return A list with:
#'   \code{vcov} (the variance-covariance matrix, or \code{NA} on failure;
#'   diagonals with non-positive variance, and their rows/cols, set to
#'   \code{NA}); \code{rcond} (reciprocal condition number of the symmetrized
#'   Hessian, \code{NA} if unavailable); \code{pd} (\code{TRUE} if the Hessian
#'   was positive-definite, \code{FALSE} if inverted via fallback, \code{NA}
#'   if not invertible).
#' @noRd
.hzr_safe_solve <- function(H, tol = .hzr_rcond_tol) {
  # (1) Non-finite / non-matrix guard
  if (is.null(H) || !is.matrix(H) || anyNA(H) || any(!is.finite(H))) {
    warning("Hessian contains non-finite entries; standard errors unavailable")
    return(list(vcov = NA, rcond = NA_real_, pd = NA))
  }

  # (2) Symmetrize (numDeriv Hessians are only symmetric to Richardson tol)
  H <- (H + t(H)) / 2

  # (3) Conditioning check
  rc <- tryCatch(rcond(H), error = function(e) NA_real_)
  if (is.na(rc)) {
    warning("Hessian conditioning could not be assessed; standard errors may be unreliable")
  } else if (rc < tol) {
    warning(sprintf(
      "Hessian is ill-conditioned (rcond = %.3g); standard errors may be unreliable",
      rc
    ))
  }

  # (4) Stable inversion: Cholesky (PD) with solve() fallback (non-PD)
  ch <- tryCatch(chol(H), error = function(e) NULL)
  if (!is.null(ch)) {
    pd <- TRUE
    vcov <- chol2inv(ch)
  } else {
    pd <- FALSE
    vcov <- tryCatch(solve(H), error = function(e) NULL)
    if (is.null(vcov)) {
      warning("Hessian not invertible; standard errors unavailable")
      return(list(vcov = NA, rcond = rc, pd = NA))
    }
    warning("Hessian is not positive-definite at the optimum; standard errors may be unreliable")
  }
  dimnames(vcov) <- dimnames(H)  # chol2inv() drops names

  # (5) Non-positive-variance guard
  d <- diag(vcov)
  bad <- !is.finite(d) | d <= 0
  if (any(bad)) {
    warning("Non-positive variance estimates; the optimum may not be a proper maximum")
    vcov[bad, ] <- NA_real_
    vcov[, bad] <- NA_real_
  }

  list(vcov = vcov, rcond = rc, pd = pd)
}
