#' Log Hazard Ratio Estimate from a Risk Set Summary
#'
#' @description
#' Internal helper returning the estimated log hazard ratio (experimental
#' versus control) from the event-time summary of \code{.rcpgsd_riskset()},
#' by Cox regression (maximum partial likelihood, \code{.rcpgsd_cox1()}),
#' by the method of Peto, \eqn{(O_T - E_T) / V}, or by the method of Pike,
#' \eqn{\log\{(O_T / E_T) / (O_C / E_C)\}} (\code{.rcpgsd_logrank()}). The
#' same estimator is used for the overall population and for each region in
#' \code{rcpgsd_sim()}.
#'
#' For all three estimators, \code{NA} is returned when no event occurs at a
#' time at which both groups are at risk. The Cox and Pike estimates are
#' infinite when all events occur in one group (\code{-Inf} when they all
#' occur in the control group); the Peto estimate is always finite.
#'
#' @param rs List returned by \code{.rcpgsd_riskset()}.
#' @param estimator Character. \code{"cox"}, \code{"peto"} or \code{"pike"}.
#'
#' @return Numeric scalar. Estimated log hazard ratio.
#'
#' @noRd
.rcpgsd_log_hr <- function(rs, estimator = "cox") {
  if (!any(rs$n_t > 0 & rs$n_c > 0)) {
    return(NA_real_)
  }
  if (estimator == "cox") {
    return(.rcpgsd_cox1(rs))
  }
  est <- .rcpgsd_logrank(rs, estimator)$log_hr
  if (is.nan(est)) NA_real_ else est
}
