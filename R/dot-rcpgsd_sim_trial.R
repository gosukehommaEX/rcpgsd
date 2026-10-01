#' Generate Patient-Level Data for One Simulated Trial
#'
#' @description
#' Internal helper that generates entry times, treatment assignments and
#' event times for one multi-regional trial. Within region \eqn{s} and
#' accrual interval \eqn{j}, \code{n_accrual[s, j]} patients enter uniformly
#' on the interval (piecewise uniform accrual). Within each region,
#' \code{round(n_s r / (1 + r))} patients are randomly assigned to the
#' experimental group and the others to the control group. Event times are
#' exponential with hazard \code{lambda_c[s]} in the control group and
#' \code{lambda_c[s] * hr[s]} in the experimental group.
#'
#' @param n_accrual Numeric \code{S}-by-\code{J} matrix of numbers of
#'   patients by region and accrual interval.
#' @param accrual_time Numeric vector of length \code{J + 1} of accrual
#'   interval boundaries.
#' @param lambda_c Numeric vector of length \code{S}. Control hazards.
#' @param hr Numeric vector of length \code{S}. Hazard ratios.
#' @param r Numeric scalar. Allocation ratio.
#'
#' @return A list with elements \code{region}, \code{entry}, \code{trt} and
#'   \code{t_event} (time from entry to event), one element per patient.
#'
#' @importFrom stats runif rexp
#'
#' @noRd
.rcpgsd_sim_trial <- function(n_accrual, accrual_time, lambda_c, hr, r) {
  S <- nrow(n_accrual)
  J <- ncol(n_accrual)
  region <- vector("list", S)
  entry  <- vector("list", S)
  trt    <- vector("list", S)
  for (s in seq_len(S)) {
    e_s <- unlist(lapply(seq_len(J), function(j) {
      stats::runif(n_accrual[s, j], accrual_time[j], accrual_time[j + 1L])
    }))
    n_s <- length(e_s)
    n_t <- round(n_s * r / (1 + r))
    a_s <- rep(c(1L, 0L), c(n_t, n_s - n_t))
    region[[s]] <- rep(s, n_s)
    entry[[s]]  <- e_s
    trt[[s]]    <- a_s[sample.int(n_s)]
  }
  region <- unlist(region)
  trt    <- unlist(trt)
  haz    <- lambda_c[region] * hr[region]^trt
  list(region  = region,
       entry   = unlist(entry),
       trt     = trt,
       t_event = stats::rexp(length(region), haz))
}
