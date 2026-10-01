#' Expected Numbers of Events and Patients by Region at a Calendar Time
#'
#' @description
#' Internal helper returning, for each region, the expected number of events
#' and the expected number of enrolled patients by calendar time \eqn{t},
#' under piecewise uniform accrual, exponential event times and exponential
#' drop-out. For \eqn{n} patients entering uniformly on \eqn{(a, b)} with
#' event hazard \eqn{\lambda} and drop-out hazard \eqn{\eta}, and with
#' \eqn{\theta = \lambda + \eta}, the expected number of events by time
#' \eqn{t > a} is
#' \deqn{\frac{\lambda}{\theta} \frac{n}{b - a} \left[(c - a) -
#' \frac{e^{-\theta (t - c)} - e^{-\theta (t - a)}}{\theta}\right],
#' \quad c = \min(b, t).}
#' Without drop-out (\eqn{\eta = 0}) the factor \eqn{\lambda / \theta}
#' equals 1.
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
#' @param lambda_d Numeric \code{S}-by-2 matrix. Drop-out hazards in the
#'   control (first column) and experimental (second column) groups.
#'   Default is no drop-out.
#'
#' @return A list with elements \code{events} and \code{enrolled}, numeric
#'   vectors of length \code{S}.
#'
#' @noRd
.rcpgsd_expected_events <- function(t, n_accrual, accrual_time, lambda_c,
                                    hr, r,
                                    lambda_d = matrix(0, nrow(n_accrual),
                                                      2L)) {
  S   <- nrow(n_accrual)
  J   <- ncol(n_accrual)
  n_s <- rowSums(n_accrual)
  p_t <- round(n_s * r / (1 + r)) / n_s
  events   <- numeric(S)
  enrolled <- numeric(S)
  for (s in seq_len(S)) {
    haz  <- c(lambda_c[s] * hr[s], lambda_c[s])
    eta  <- c(lambda_d[s, 2L], lambda_d[s, 1L])
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
        th  <- lam + eta[g]
        ev  <- lam / th * ((cc - a) - (exp(-th * (t - cc)) -
                                         exp(-th * (t - a))) / th) / (b - a)
        events[s] <- events[s] + n * prop[g] * ev
      }
    }
  }
  list(events = events, enrolled = enrolled)
}
