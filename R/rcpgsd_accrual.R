#' Regional Accrual with Delayed Starts under a Fixed Overall Accrual Plan
#'
#' @description
#' Constructs the piecewise uniform accrual by region (\code{accrual_time}
#' and \code{n_accrual} of \code{\link{rcpgsd_events}} and
#' \code{\link{rcpgsd_sim}}) when some regions start enrolling later than
#' others, while the overall accrual plan of the protocol (total number of
#' patients, accrual period and overall accrual rate) is kept unchanged.
#' The delay in a region is compensated for by the regions that start at
#' time 0, so that the overall accrual over time, and hence the expected
#' analysis times of an event-driven trial with hazards common to all
#' regions, are those of the protocol.
#'
#' @details
#' Let \eqn{T} be the end of the overall accrual period, \eqn{R(t)} the
#' overall accrual rate of the protocol (piecewise constant), and
#' \eqn{n_s} and \eqn{\tau_s} the number of patients and the start of
#' enrolment in region \eqn{s}. The accrual is constructed as follows.
#' \itemize{
#'   \item A delayed region (\eqn{\tau_s > 0}) enrols its \eqn{n_s}
#'     patients uniformly between \eqn{\tau_s} and \eqn{T}, that is, at the
#'     rate \eqn{n_s / (T - \tau_s)}.
#'   \item The regions that start at time 0 share the remaining overall
#'     accrual rate in proportion to their numbers of patients, so that
#'     each of them also enrols exactly its \eqn{n_s} patients.
#' }
#' This is possible if and only if, at every time \eqn{t}, the delayed
#' regions that have started do not need more than the overall accrual
#' rate:
#' \deqn{\sum_{s: 0 < \tau_s < t} \frac{n_s}{T - \tau_s} \le R(t).}
#' For a single delayed region and a constant overall accrual rate
#' \eqn{N / T}, this reduces to \eqn{\tau_s \le (1 - n_s / N) T}. When the
#' condition is violated, the function stops with an error that reports the
#' first time interval in which the overall accrual cannot absorb the
#' delayed regions.
#'
#' The numbers of patients by region and interval are generally not
#' integers. The returned \code{n_accrual} is rounded so that the
#' cumulative number of patients in each region is rounded at every
#' interval boundary; the regional totals are therefore exactly
#' \code{n_region}, and the cumulative numbers differ from the exact values
#' by at most 0.5 patient. The exact values are returned as
#' \code{n_accrual_exact}.
#'
#' @param n_region Numeric vector of length \code{S}. Planned numbers of
#'   patients by region (positive integers); region 1 is the region of
#'   interest. Their sum must equal the total number of patients of the
#'   overall accrual plan.
#' @param start Numeric scalar or vector of length \code{S}. Start of
#'   enrolment in each region (non-negative and smaller than the end of the
#'   overall accrual period). At least one region must start at time 0.
#'   Default is \code{0} (no delay).
#' @param accrual_time Numeric vector of length \code{J + 1}. Boundaries of
#'   the intervals of the overall accrual plan, starting at 0.
#' @param n_overall Numeric vector of length \code{J}. Number of patients in
#'   the overall population entering during each interval. Exactly one of
#'   \code{n_overall} and \code{rate_overall} must be given.
#' @param rate_overall Numeric vector of length \code{J}. Overall accrual
#'   rate (patients per unit of time) in each interval.
#'
#' @return An object of class \code{"rcpgsd_accrual"}, a list containing:
#' \describe{
#'   \item{\code{accrual_time}}{Boundaries of the accrual intervals, the
#'     union of the boundaries of the overall plan and the start times.}
#'   \item{\code{n_accrual}}{\code{S}-by-\code{M} integer matrix of numbers
#'     of patients by region and interval, to be passed to
#'     \code{\link{rcpgsd_events}} or \code{\link{rcpgsd_sim}} together
#'     with \code{accrual_time}.}
#'   \item{\code{n_accrual_exact}}{The same numbers before rounding.}
#'   \item{\code{rate}}{\code{S}-by-\code{M} matrix of accrual rates by
#'     region and interval.}
#'   \item{\code{rate_overall}}{Overall accrual rate in each interval of
#'     \code{accrual_time}.}
#'   \item{\code{n_region}, \code{start}}{Inputs.}
#'   \item{\code{plan}}{List with the overall accrual plan
#'     (\code{accrual_time}, \code{n_overall}, \code{rate_overall}).}
#' }
#'
#' @seealso \code{\link{rcpgsd_events}}, \code{\link{rcpgsd_sim}},
#'   \code{\link{print.rcpgsd_accrual}}
#'
#' @examples
#' # 600 patients at 50 patients per month over 12 months; Japan (region 1,
#' # 60 patients) starts 6 months late, EU and US start at month 0
#' acc <- rcpgsd_accrual(n_region = c(60, 270, 270), start = c(6, 0, 0),
#'                       accrual_time = c(0, 12), rate_overall = 50)
#' print(acc)
#'
#' # Expected analysis times are those of the plan without delay, but the
#' # proportion of Japanese events is smaller
#' ev  <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
#'                      lambda_c = log(2) / 12,
#'                      accrual_time = acc$accrual_time,
#'                      n_accrual = acc$n_accrual)
#' acc0 <- rcpgsd_accrual(n_region = c(60, 270, 270),
#'                        accrual_time = c(0, 12), rate_overall = 50)
#' ev0 <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
#'                      lambda_c = log(2) / 12,
#'                      accrual_time = acc0$accrual_time,
#'                      n_accrual = acc0$n_accrual)
#' rbind(time_delay = ev$analysis$time, time_no_delay = ev0$analysis$time,
#'       fs1_delay = ev$fs[1, ], fs1_no_delay = ev0$fs[1, ])
#'
#' # A delay of 11 months cannot be absorbed (11 > (1 - 60 / 600) * 12)
#' try(rcpgsd_accrual(n_region = c(60, 270, 270), start = c(11, 0, 0),
#'                    accrual_time = c(0, 12), rate_overall = 50))
#'
#' @export
rcpgsd_accrual <- function(n_region,
                           start        = 0,
                           accrual_time,
                           n_overall    = NULL,
                           rate_overall = NULL) {

  # ========== Argument checks ==========
  if (!is.numeric(n_region) || length(n_region) < 2L ||
      any(!is.finite(n_region)) || any(n_region < 1) ||
      any(n_region != round(n_region))) {
    stop("'n_region' must be a vector of at least two positive integers.",
         call. = FALSE)
  }
  S <- length(n_region)

  if (!is.numeric(accrual_time) || length(accrual_time) < 2L ||
      any(!is.finite(accrual_time)) || accrual_time[1L] != 0 ||
      any(diff(accrual_time) <= 0)) {
    stop("'accrual_time' must be a strictly increasing numeric vector ",
         "starting at 0.", call. = FALSE)
  }
  J     <- length(accrual_time) - 1L
  len   <- diff(accrual_time)
  T_end <- accrual_time[J + 1L]

  if (is.null(n_overall) == is.null(rate_overall)) {
    stop("Exactly one of 'n_overall' and 'rate_overall' must be given.",
         call. = FALSE)
  }
  if (!is.null(n_overall)) {
    if (!is.numeric(n_overall) || length(n_overall) != J ||
        any(!is.finite(n_overall)) || any(n_overall < 0)) {
      stop(sprintf(paste0("'n_overall' must be a non-negative numeric ",
                          "vector of length %d (one element per interval ",
                          "of 'accrual_time')."), J), call. = FALSE)
    }
    rate_overall <- n_overall / len
  } else {
    if (!is.numeric(rate_overall) || length(rate_overall) != J ||
        any(!is.finite(rate_overall)) || any(rate_overall < 0)) {
      stop(sprintf(paste0("'rate_overall' must be a non-negative numeric ",
                          "vector of length %d (one element per interval ",
                          "of 'accrual_time')."), J), call. = FALSE)
    }
    n_overall <- rate_overall * len
  }
  n_total <- sum(n_overall)
  if (abs(sum(n_region) - n_total) > 1e-8 * max(1, n_total)) {
    stop(sprintf(paste0("The sum of 'n_region' (%s) must equal the total ",
                        "number of patients of the overall accrual plan ",
                        "(%s)."),
                 format(sum(n_region)), format(signif(n_total, 10))),
         call. = FALSE)
  }

  if (!is.numeric(start) || !(length(start) %in% c(1L, S)) ||
      any(!is.finite(start)) || any(start < 0) || any(start >= T_end)) {
    stop(sprintf(paste0("'start' must be a numeric scalar or a vector of ",
                        "length %d with values in [0, %s)."), S,
                 format(T_end)), call. = FALSE)
  }
  start   <- rep_len(as.numeric(start), S)
  delayed <- start > 0
  if (all(delayed)) {
    stop("At least one region must start enrolling at time 0 to keep the ",
         "overall accrual plan.", call. = FALSE)
  }

  # ========== Rates on the merged grid ==========
  grid <- sort(unique(c(accrual_time, start)))
  M    <- length(grid) - 1L
  mid  <- (grid[-1L] + grid[-(M + 1L)]) / 2
  r_ov <- rate_overall[findInterval(mid, accrual_time)]

  rate <- matrix(0, S, M)
  for (s in which(delayed)) {
    rate[s, mid > start[s]] <- n_region[s] / (T_end - start[s])
  }
  need   <- colSums(rate)
  excess <- need - r_ov
  bad    <- which(excess > 1e-10 * max(1, max(r_ov)))
  if (length(bad) > 0L) {
    j   <- bad[1L]
    act <- which(delayed & start < mid[j])
    one <- which(delayed)
    hint <- if (length(one) == 1L && J == 1L) {
      sprintf(paste0(" With a single delayed region and a constant overall ",
                     "accrual rate, the start must not exceed (1 - %s / %s) ",
                     "* %s = %s."),
              format(n_region[one]), format(signif(n_total, 10)),
              format(T_end),
              format(signif((1 - n_region[one] / n_total) * T_end, 6)))
    } else {
      ""
    }
    stop(sprintf(paste0("The delayed regions cannot be compensated for by ",
                        "the regions that start at time 0: between times ",
                        "%s and %s, region(s) %s enrolling uniformly until ",
                        "the end of accrual need %s patients per unit of ",
                        "time, but the overall accrual rate is %s.%s"),
                 format(signif(grid[j], 6)), format(signif(grid[j + 1L], 6)),
                 paste(act, collapse = ", "), format(signif(need[j], 6)),
                 format(signif(r_ov[j], 6)), hint),
         call. = FALSE)
  }

  # Regions starting at time 0 share the remaining overall rate
  absorb <- which(!delayed)
  w      <- n_region[absorb] / sum(n_region[absorb])
  rate[absorb, ] <- outer(w, pmax(r_ov - need, 0))

  # ========== Numbers of patients ==========
  n_exact <- sweep(rate, 2L, diff(grid), "*")
  n_int   <- n_exact
  for (s in seq_len(S)) {
    n_int[s, ] <- diff(c(0, round(cumsum(n_exact[s, ]))))
  }

  structure(
    list(accrual_time    = grid,
         n_accrual       = n_int,
         n_accrual_exact = n_exact,
         rate            = rate,
         rate_overall    = r_ov,
         n_region        = as.numeric(n_region),
         start           = start,
         plan            = list(accrual_time = as.numeric(accrual_time),
                                n_overall    = n_overall,
                                rate_overall = rate_overall)),
    class = "rcpgsd_accrual"
  )
}
