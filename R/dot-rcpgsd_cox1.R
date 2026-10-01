#' Cox Regression Estimate for a Single Binary Covariate
#'
#' @description
#' Internal helper returning the maximum partial likelihood estimate of the
#' log hazard ratio (experimental versus control) in a Cox proportional
#' hazards model with a single binary covariate, for data without tied event
#' times. The score function is
#' \eqn{U(\beta) = \sum_j \{d_{T,j} - n_{T,j} e^{\beta} /
#' (n_{T,j} e^{\beta} + n_{C,j})\}}, summed over event times at which both
#' groups are at risk (other event times do not contribute). It is solved by
#' Newton's method starting from the Peto estimate, with steps limited to 1
#' and bisection within the bracketing interval as a safeguard.
#'
#' When no event occurs at a time when both groups are at risk, the partial
#' likelihood is flat and \code{NA} is returned. When all such events occur
#' in the control group (or all in the experimental group), the partial
#' likelihood is monotone and \code{-Inf} (or \code{Inf}) is returned.
#'
#' @param rs List returned by \code{.rcpgsd_riskset()}.
#' @param tol Numeric scalar. Convergence tolerance for the estimate.
#' @param max_iter Integer. Maximum number of iterations.
#'
#' @return Numeric scalar. Estimated log hazard ratio.
#'
#' @noRd
.rcpgsd_cox1 <- function(rs, tol = 1e-10, max_iter = 100L) {
  info_ev <- rs$n_t > 0 & rs$n_c > 0
  if (!any(info_ev)) {
    return(NA_real_)
  }
  n_t <- rs$n_t[info_ev]
  n_c <- rs$n_c[info_ev]
  o_t <- sum(rs$d_t[info_ev])
  m   <- length(n_t)
  if (o_t == 0) {
    return(-Inf)
  }
  if (o_t == m) {
    return(Inf)
  }

  p0 <- n_t / (n_t + n_c)
  b  <- (o_t - sum(p0)) / sum(p0 * (1 - p0))
  lo <- -Inf
  hi <- Inf
  for (iter in seq_len(max_iter)) {
    e    <- n_t * exp(b)
    p    <- e / (e + n_c)
    u    <- o_t - sum(p)
    info <- sum(p * (1 - p))
    if (u > 0) lo <- b else hi <- b
    step <- if (info > 0) u / info else sign(u)
    step <- max(-1, min(1, step))
    b_new <- b + step
    if (is.finite(lo) && is.finite(hi) && (b_new <= lo || b_new >= hi)) {
      b_new <- (lo + hi) / 2
    }
    if (abs(b_new - b) < tol) {
      return(b_new)
    }
    b <- b_new
  }
  b
}
