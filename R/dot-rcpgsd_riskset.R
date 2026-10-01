#' Risk Set Summary at Event Times for Two Groups
#'
#' @description
#' Internal helper that returns, for each observed event time, the numbers of
#' patients at risk in the experimental and control groups and the group of
#' the patient with the event. Event times are assumed to be distinct
#' (continuous survival times without ties).
#'
#' @param time Numeric vector of observed times.
#' @param status Integer vector of event indicators (1 = event,
#'   0 = censored).
#' @param trt Integer vector of group indicators (1 = experimental,
#'   0 = control).
#'
#' @return A list with elements \code{n_t} and \code{n_c} (numbers at risk
#'   in the experimental and control groups just before each event time) and
#'   \code{d_t} (1 if the event occurred in the experimental group), each of
#'   length equal to the number of events.
#'
#' @noRd
.rcpgsd_riskset <- function(time, status, trt) {
  n <- length(time)
  if (n == 0L) {
    return(list(n_t = numeric(0), n_c = numeric(0), d_t = numeric(0)))
  }
  o      <- order(time)
  status <- status[o]
  trt    <- trt[o]
  n_t    <- rev(cumsum(rev(trt)))
  n_c    <- rev(seq_len(n)) - n_t
  ev     <- status == 1L
  list(n_t = as.numeric(n_t[ev]),
       n_c = as.numeric(n_c[ev]),
       d_t = as.numeric(trt[ev]))
}
