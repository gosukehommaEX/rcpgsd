#' Expected Analysis Times and Regional Event Proportions
#'
#' @description
#' Computes, for an event-driven multi-regional trial, the expected calendar
#' time of each analysis and the expected number and proportion of events in
#' each region at that time, from the accrual pattern by region and
#' exponential event times. The resulting proportions can be passed to
#' \code{\link{rcpgsd}} or \code{\link{rcpgsd_get_fs}} as \code{fs}, so that
#' delayed enrolment in a region is reflected in the regional consistency
#' probabilities.
#'
#' @details
#' Patients enter by piecewise uniform accrual (\code{n_accrual[s, j]}
#' patients of region \eqn{s} uniformly within accrual interval \eqn{j}) and
#' have exponential event times with hazard \code{lambda_c[s]} in the
#' control group and \code{lambda_c[s] * hr[s]} in the experimental group.
#' Patients may drop out at exponential times with the hazards given by
#' \code{lambda_d}; a patient who drops out before the event is censored at
#' the drop-out time, so only a fraction \eqn{\lambda / (\lambda + \eta)}
#' of the patients with event hazard \eqn{\lambda} and drop-out hazard
#' \eqn{\eta} eventually has an event. These are the assumptions of
#' \code{\link{rcpgsd_sim}}.
#' The expected number of events in each region by calendar time \eqn{t} is
#' obtained in closed form, and the time of analysis \eqn{k} is defined as
#' the time at which the expected total number of events equals
#' \code{events[k]} (found with \code{stats::uniroot()}). The regional
#' proportions are the expected regional numbers of events at that time
#' divided by \code{events[k]}.
#'
#' These are expectations at the expected analysis time. They differ
#' slightly from the mean regional proportions at the random analysis times
#' of a simulated trial (\code{fs_mean} of \code{\link{rcpgsd_sim}}), but the
#' difference is usually small.
#'
#' @inheritParams rcpgsd_sim
#'
#' @return An object of class \code{"rcpgsd_events"}, a list containing:
#' \describe{
#'   \item{\code{analysis}}{Data frame with one row per analysis and columns
#'     \code{stage}, \code{events}, \code{time} (expected calendar time) and
#'     \code{enrolled} (expected number of enrolled patients).}
#'   \item{\code{events_region}}{\code{S}-by-\code{K} matrix of expected
#'     numbers of events by region and analysis.}
#'   \item{\code{enrolled_region}}{\code{S}-by-\code{K} matrix of expected
#'     numbers of enrolled patients by region and analysis.}
#'   \item{\code{fs}}{\code{S}-by-\code{K} matrix of expected regional event
#'     proportions, with columns summing to 1.}
#'   \item{\code{events}, \code{hr}, \code{lambda_c}, \code{accrual_time},
#'     \code{n_accrual}, \code{lambda_d}, \code{r}}{Inputs
#'     (\code{lambda_d} as an \code{S}-by-2 matrix with the control group
#'     in the first column).}
#' }
#'
#' @seealso \code{\link{rcpgsd}}, \code{\link{rcpgsd_get_fs}},
#'   \code{\link{rcpgsd_sim}}, \code{\link{print.rcpgsd_events}}
#'
#' @examples
#' # 600 patients over 18 months; Japan (region 1) starts after month 6
#' n_acc <- rbind(c(0, 60), c(100, 170), c(100, 170))
#' ev <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
#'                     lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
#'                     n_accrual = n_acc)
#' print(ev)
#'
#' # Regional consistency probabilities with the expected proportions
#' design <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
#' res <- rcpgsd(design, events = c(150, 300, 450), fs = ev$fs, hr = 0.72)
#' res$overall
#'
#' # A drop-out probability of 0.05 per year in both groups delays the
#' # analyses
#' ev_d <- rcpgsd_events(events = c(150, 300, 450), hr = 0.72,
#'                       lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
#'                       n_accrual = n_acc, lambda_d = -log(1 - 0.05) / 12)
#' rbind(no_dropout = ev$analysis$time, dropout = ev_d$analysis$time)
#'
#' @importFrom stats uniroot
#'
#' @export
rcpgsd_events <- function(events,
                          hr,
                          lambda_c,
                          accrual_time,
                          n_accrual,
                          lambda_d = 0,
                          r        = 1) {
  K   <- length(events)
  inp <- .rcpgsd_check_sim_inputs(events, hr, lambda_c, accrual_time,
                                  n_accrual, r, K, lambda_d)
  S <- nrow(inp$n_accrual)
  if (inp$events[K] >= sum(inp$n_accrual)) {
    stop("The number of events at the final analysis must be smaller than ",
         "the total number of patients.", call. = FALSE)
  }

  expected <- function(t) {
    .rcpgsd_expected_events(t, inp$n_accrual, inp$accrual_time,
                            inp$lambda_c, inp$hr, inp$r, inp$lambda_d)
  }

  time     <- numeric(K)
  ev_reg   <- matrix(0, S, K)
  enr_reg  <- matrix(0, S, K)
  for (k in seq_len(K)) {
    d_k <- inp$events[k]
    f <- function(t) sum(expected(t)$events) - d_k
    hi <- max(inp$accrual_time)
    while (f(hi) < 0) {
      hi <- 2 * hi
    }
    time[k] <- stats::uniroot(f, c(0, hi), tol = 1e-10)$root
    e_k <- expected(time[k])
    ev_reg[, k]  <- e_k$events
    enr_reg[, k] <- e_k$enrolled
  }
  fs <- sweep(ev_reg, 2L, colSums(ev_reg), "/")

  structure(
    list(analysis        = data.frame(stage    = seq_len(K),
                                      events   = inp$events,
                                      time     = time,
                                      enrolled = colSums(enr_reg)),
         events_region   = ev_reg,
         enrolled_region = enr_reg,
         fs              = fs,
         events          = inp$events,
         hr              = inp$hr,
         lambda_c        = inp$lambda_c,
         accrual_time    = inp$accrual_time,
         n_accrual       = inp$n_accrual,
         lambda_d        = inp$lambda_d,
         r               = inp$r),
    class = "rcpgsd_events"
  )
}
