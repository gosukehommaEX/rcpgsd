#' Log-Rank Statistic and Pike or Peto Hazard Ratio Estimate
#'
#' @description
#' Internal helper computing, from a risk set summary, the observed and
#' expected numbers of events in the experimental group, the hypergeometric
#' variance, the log-rank Z-statistic oriented so that positive values favour
#' the experimental group, and the log hazard ratio estimate by the method of
#' Pike or of Peto.
#'
#' With \eqn{O_T}, \eqn{E_T}, \eqn{O_C = d - O_T}, \eqn{E_C = d - E_T} and
#' variance \eqn{V}, the Pike estimate is
#' \eqn{\log \widehat{\mathrm{HR}} = \log(O_T / E_T) - \log(O_C / E_C)} and
#' the Peto estimate is \eqn{(O_T - E_T) / V}. The log-rank statistic is
#' \eqn{Z = (E_T - O_T) / \sqrt{V}}.
#'
#' @param rs List returned by \code{.rcpgsd_riskset()}.
#' @param estimator Character. \code{"pike"} or \code{"peto"}.
#'
#' @return A list with elements \code{d}, \code{o_t}, \code{e_t}, \code{v},
#'   \code{z} and \code{log_hr}.
#'
#' @noRd
.rcpgsd_logrank <- function(rs, estimator = "pike") {
  d <- length(rs$d_t)
  if (d == 0L) {
    return(list(d = 0L, o_t = 0, e_t = 0, v = 0, z = NA_real_,
                log_hr = NA_real_))
  }
  n_risk <- rs$n_t + rs$n_c
  o_t <- sum(rs$d_t)
  e_t <- sum(rs$n_t / n_risk)
  v   <- sum(rs$n_t * rs$n_c / n_risk^2)
  z   <- if (v > 0) (e_t - o_t) / sqrt(v) else NA_real_
  log_hr <- if (estimator == "pike") {
    log(o_t / e_t) - log((d - o_t) / (d - e_t))
  } else {
    if (v > 0) (o_t - e_t) / v else NA_real_
  }
  list(d = d, o_t = o_t, e_t = e_t, v = v, z = z, log_hr = log_hr)
}
