#' Expected Numbers of Events and Patients by Region at a Calendar Time
#'
#' @description
#' Internal helper returning, for each region, the expected number of events
#' and the expected number of enrolled patients by calendar time \eqn{t},
#' under piecewise uniform accrual and exponential event times without
#' drop-out. For \eqn{n} patients entering uniformly on \eqn{(a, b)} with
#' hazard \eqn{\lambda}, the expected number of events by time \eqn{t > a}
#' is
#' \deqn{\frac{n}{b - a} \left[(c - a) - \frac{e^{-\lambda (t - c)} -
#' e^{-\lambda (t - a)}}{\lambda}\right], \quad c = \min(b, t).}
#' Within each region the proportion of patients in the experimental group
#' is \code{round(n_s r / (1 + r)) / n_s}, as in \code{\link{rcpgsd_sim}}.
#'
#' @param t Numeric scalar. Calendar time.
#' @param n_accrual Numeric \code{S}-by-\code{J} matrix of numbers of
#'   patients by region and accrual interval.
#' @param accrual_time Numeric vector of length \code{J + 1} of accrual
#'   interval boundaries.
#' @param lambda_c Numeric vector of length \code{S}. Control hazards.
#' @param hr Numeric vector of length \code{S}. Hazard ratios.
#' @param r Numeric scalar. Allocation ratio.
#'
#' @return A list with elements \code{events} and \code{enrolled}, numeric
#'   vectors of length \code{S}.
#'
#' @noRd
.rcpgsd_expected_events <- function(t, n_accrual, accrual_time, lambda_c,
                                    hr, r) {
  S   <- nrow(n_accrual)
  J   <- ncol(n_accrual)
  n_s <- rowSums(n_accrual)
  p_t <- round(n_s * r / (1 + r)) / n_s
  events   <- numeric(S)
  enrolled <- numeric(S)
  for (s in seq_len(S)) {
    haz  <- c(lambda_c[s] * hr[s], lambda_c[s])
    prop <- c(p_t[s], 1 - p_t[s])
    for (j in seq_len(J)) {
      a <- accrual_time[j]
      b <- accrual_time[j + 1L]
      n <- n_accrual[s, j]
      if (t <= a || n == 0) next
      cc <- min(b, t)
      enrolled[s] <- enrolled[s] + n * (cc - a) / (b - a)
      for (g in 1:2) {
        lam <- haz[g]
        ev  <- ((cc - a) - (exp(-lam * (t - cc)) - exp(-lam * (t - a))) /
                  lam) / (b - a)
        events[s] <- events[s] + n * prop[g] * ev
      }
    }
  }
  list(events = events, enrolled = enrolled)
}
