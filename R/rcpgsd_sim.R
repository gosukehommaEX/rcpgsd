#' Patient-Level Simulation of Regional Consistency under a Group Sequential
#' Design
#'
#' @description
#' Estimates by Monte Carlo simulation the same quantities as
#' \code{\link{rcpgsd}} for an event-driven multi-regional trial, using
#' patient-level time-to-event data. The simulation does not rely on the
#' normal approximation of the regional log hazard ratio estimates, so its
#' results can be compared with the closed-form results of
#' \code{\link{rcpgsd}} to assess the accuracy of the approximation, in
#' particular when the number of events in the region of interest is small.
#'
#' @details
#' \strong{Data generation.} Patients enter the trial by piecewise uniform
#' accrual: \code{n_accrual[s, j]} patients of region \eqn{s} enter uniformly
#' within the interval from \code{accrual_time[j]} to
#' \code{accrual_time[j + 1]}. Setting some entries of \code{n_accrual} to
#' zero represents delayed enrolment in a region. Within each region,
#' \code{round(n_s r / (1 + r))} of the \eqn{n_s} patients are randomly
#' assigned to the experimental group. Event times are exponential with
#' hazard \code{lambda_c[s]} in the control group and
#' \code{lambda_c[s] * hr[s]} in the experimental group. There is no
#' drop-out.
#'
#' \strong{Event-driven analyses.} Analysis \eqn{k} takes place at the
#' calendar time at which the cumulative number of events in the trial
#' reaches \code{events[k]}. Patients who have entered by then are included,
#' and those without an event are censored at that time.
#'
#' \strong{Analysis at each analysis time.}
#' \itemize{
#'   \item The efficacy and futility boundaries are applied to the
#'     unstratified log-rank statistic, oriented so that positive values
#'     favour the experimental group.
#'   \item The overall hazard ratio used in the consistency criteria is
#'     estimated by the method of Pike (\code{overall_est = "pike"},
#'     default), \eqn{(O_T / E_T) / (O_C / E_C)}, or of Peto
#'     (\code{overall_est = "peto"}), \eqn{\exp\{(O_T - E_T) / V\}}, where
#'     \eqn{O}, \eqn{E} and \eqn{V} are the observed and expected numbers of
#'     events and the variance of the log-rank statistic.
#'   \item The hazard ratio in each region is estimated by a Cox proportional
#'     hazards model fitted to the data of that region only, with treatment
#'     as the single covariate. Without tied event times this is the exact
#'     maximum partial likelihood estimate (identical to
#'     \code{survival::coxph()}).
#' }
#' When all events of a region occur in one group, the regional estimate is
#' infinite (hazard ratio of 0 or infinity) and is used as such in the
#' consistency criteria. When no estimate is available (no events, or no
#' event while both groups are at risk), the consistency criterion is
#' regarded as not met. The frequencies of these cases are reported.
#'
#' The consistency criteria and the reported probabilities are those of
#' \code{\link{rcpgsd}}. To compare with the closed form, the mean observed
#' regional event proportions \code{fs_mean} can be passed to
#' \code{\link{rcpgsd}} (see Examples).
#'
#' @inheritParams rcpgsd
#' @param events Numeric vector of length \code{K}. Cumulative total numbers
#'   of events that trigger the analyses (positive integers, strictly
#'   increasing).
#' @param hr Numeric scalar or vector of length \code{S}. True hazard
#'   ratio(s) (experimental versus control).
#' @param lambda_c Numeric scalar or vector of length \code{S}. Hazard of the
#'   control group (per unit of time), common to all regions or
#'   region-specific.
#' @param accrual_time Numeric vector of length \code{J + 1}. Boundaries of
#'   the accrual intervals, starting at 0.
#' @param n_accrual Numeric \code{S}-by-\code{J} matrix. Number of patients
#'   entering in region \code{s} (row) during accrual interval \code{j}
#'   (column); region 1 is the region of interest. A vector of length
#'   \code{S} is accepted when there is a single accrual interval.
#' @param overall_est Character. Estimator of the overall hazard ratio used
#'   in the consistency criteria: \code{"pike"} (default) or \code{"peto"}.
#' @param nsim Positive integer. Number of simulated trials. Default is
#'   \code{10000}.
#' @param seed Integer or \code{NULL}. Random seed. Default is \code{1}.
#'
#' @return An object of class \code{c("rcpgsd_sim", "rcpgsd")}, a list
#'   containing:
#' \describe{
#'   \item{\code{stage}}{Data frame with the same columns as for
#'     \code{\link{rcpgsd}} (with \code{events_r1} the mean number of events
#'     in region 1) and, in addition, \code{se_rcp}, the Monte Carlo
#'     standard error of \code{rcp}.}
#'   \item{\code{overall}}{Data frame with columns \code{method},
#'     \code{power}, \code{p_joint} and \code{rcp}.}
#'   \item{\code{analysis}}{Data frame with one row per analysis and columns
#'     \code{stage}, \code{events}, \code{time_mean} and \code{time_sd}
#'     (calendar time of the analysis), \code{events_r1_mean},
#'     \code{p_no_est_r1} (proportion of trials without an estimate in
#'     region 1) and \code{p_inf_r1} (proportion of trials with an infinite
#'     estimate in region 1).}
#'   \item{\code{fs_mean}}{\code{S}-by-\code{K} matrix of mean regional event
#'     proportions.}
#'   \item{\code{p_no_est}, \code{p_inf}}{\code{S}-by-\code{K} matrices of the
#'     proportions of trials without an estimate and with an infinite
#'     estimate, by region and analysis.}
#'   \item{\code{design}, \code{events}, \code{hr}, \code{lambda_c},
#'     \code{accrual_time}, \code{n_accrual}, \code{r}, \code{method},
#'     \code{scale}, \code{gamma_M1}, \code{gamma_M2}, \code{overall_est},
#'     \code{nsim}}{Settings used.}
#' }
#' Because the object also inherits from class \code{"rcpgsd"}, the
#' \code{\link{plot.rcpgsd}} method can be used.
#'
#' @importFrom stats sd
#'
#' @seealso \code{\link{rcpgsd}}, \code{\link{print.rcpgsd_sim}},
#'   \code{\link{plot.rcpgsd}}
#'
#' @examples
#' # Three analyses; Japan (region 1) starts enrolling after 6 months
#' design <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
#' n_acc  <- rbind(c(0, 60), c(100, 170), c(100, 170))
#' sim <- rcpgsd_sim(design, events = c(150, 300, 450), hr = 0.72,
#'                   lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
#'                   n_accrual = n_acc, nsim = 200)
#' print(sim)
#'
#' \donttest{
#' # Comparison with the closed form using the mean regional proportions
#' sim <- rcpgsd_sim(design, events = c(150, 300, 450), hr = 0.72,
#'                   lambda_c = log(2) / 12, accrual_time = c(0, 6, 18),
#'                   n_accrual = n_acc, nsim = 5000)
#' cf <- rcpgsd(design, events = c(150, 300, 450), fs = sim$fs_mean,
#'              hr = 0.72)
#' cbind(simulation = sim$stage$rcp, closed_form = cf$stage$rcp)
#' plot(sim)
#' }
#'
#' @export
rcpgsd_sim <- function(design,
                       events,
                       hr,
                       lambda_c,
                       accrual_time,
                       n_accrual,
                       r           = 1,
                       method      = c("M1", "M2"),
                       scale       = c("loghr", "hr"),
                       gamma_M1    = 0.5,
                       gamma_M2    = 0,
                       overall_est = c("pike", "peto"),
                       futility    = c("obey", "ignore"),
                       nsim        = 10000,
                       seed        = 1) {

  # ========== Argument checks ==========
  method      <- match.arg(method, c("M1", "M2"), several.ok = TRUE)
  method      <- intersect(c("M1", "M2"), method)
  scale       <- match.arg(scale)
  overall_est <- match.arg(overall_est)
  futility    <- match.arg(futility)

  bnd <- .rcpgsd_bounds(design)
  K   <- bnd$K
  inp <- .rcpgsd_check_sim_inputs(events, hr, lambda_c, accrual_time,
                                  n_accrual, r, K)
  S   <- nrow(inp$n_accrual)

  gamma_M1 <- .rcpgsd_check_gamma(gamma_M1, K, "gamma_M1", scale)
  gamma_M2 <- .rcpgsd_check_gamma(gamma_M2, K, "gamma_M2", scale)
  if (!is.numeric(nsim) || length(nsim) != 1L || nsim < 1) {
    stop("'nsim' must be a positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)

  fut <- .rcpgsd_futility(bnd, futility)
  info_frac <- inp$events / inp$events[K]
  if (!is.null(bnd$timing) && max(abs(info_frac - bnd$timing)) > 0.01) {
    warning("The information fractions implied by 'events' differ from ",
            "the planned information fractions of 'design'; the boundaries ",
            "of 'design' are used without recalculation.", call. = FALSE)
  }

  # ========== Simulation ==========
  if (!is.null(seed)) set.seed(seed)
  z       <- matrix(NA_real_, nsim, K)
  est_all <- matrix(NA_real_, nsim, K)
  t_anal  <- matrix(NA_real_, nsim, K)
  est_reg <- array(NA_real_, c(nsim, S, K))
  n_ev    <- array(0, c(nsim, S, K))
  ev_idx  <- as.integer(inp$events)

  for (i in seq_len(nsim)) {
    pt  <- .rcpgsd_sim_trial(inp$n_accrual, inp$accrual_time, inp$lambda_c,
                             inp$hr, inp$r)
    cal <- pt$entry + pt$t_event
    t_k <- sort(cal, partial = ev_idx)[ev_idx]
    for (k in seq_len(K)) {
      inc <- pt$entry < t_k[k]
      obs <- pmin(pt$t_event[inc], t_k[k] - pt$entry[inc])
      st  <- as.integer(cal[inc] <= t_k[k])
      tr  <- pt$trt[inc]
      rg  <- pt$region[inc]
      lr  <- .rcpgsd_logrank(.rcpgsd_riskset(obs, st, tr), overall_est)
      z[i, k]       <- lr$z
      est_all[i, k] <- -lr$log_hr
      t_anal[i, k]  <- t_k[k]
      for (s in seq_len(S)) {
        sel <- rg == s
        n_ev[i, s, k]    <- sum(st[sel])
        est_reg[i, s, k] <- -.rcpgsd_cox1(
          .rcpgsd_riskset(obs[sel], st[sel], tr[sel]))
      }
    }
  }

  # ========== Consistency and stopping ==========
  consistent <- function(m, k) {
    if (m == "M1") {
      x1 <- est_reg[, 1L, k]
      v  <- est_all[, k]
      g  <- gamma_M1[k]
      out <- if (scale == "loghr") {
        x1 > g * v
      } else {
        (1 - exp(-x1)) > g * (1 - exp(-v))
      }
    } else {
      thr <- if (scale == "loghr") gamma_M2[k] else -log(1 - gamma_M2[k])
      out <- rowSums(matrix(est_reg[, , k], nrow = nsim) > thr) == S
    }
    out[is.na(out)] <- FALSE
    out
  }

  tab <- lapply(method, function(m) {
    alive <- rep(TRUE, nsim)
    rows  <- vector("list", K)
    for (k in seq_len(K)) {
      eff  <- alive & !is.na(z[, k]) & z[, k] >= bnd$upper[k]
      futl <- if (k < K) {
        alive & !is.na(z[, k]) & z[, k] <= fut$lower[k]
      } else {
        rep(FALSE, nsim)
      }
      cons <- consistent(m, k)
      rows[[k]] <- data.frame(method     = m,
                              stage      = k,
                              gamma      = if (m == "M1") gamma_M1[k] else
                                gamma_M2[k],
                              p_eff      = mean(eff),
                              p_fut      = mean(futl),
                              p_joint    = mean(eff & cons),
                              rcp_uncond = mean(cons),
                              stringsAsFactors = FALSE)
      alive <- alive & !eff & !futl
    }
    do.call(rbind, rows)
  })
  tab <- do.call(rbind, tab)

  n_ev_mean <- matrix(apply(n_ev, c(2L, 3L), mean), S, K)
  fs_mean   <- sweep(n_ev_mean, 2L, inp$events, "/")
  p_no_est  <- matrix(apply(is.na(est_reg), c(2L, 3L), mean), S, K)
  p_inf     <- matrix(apply(is.infinite(est_reg), c(2L, 3L), mean), S, K)

  rcp    <- ifelse(tab$p_eff > 0, tab$p_joint / tab$p_eff, NA_real_)
  se_rcp <- ifelse(tab$p_eff > 0,
                   sqrt(rcp * (1 - rcp) / (nsim * tab$p_eff)), NA_real_)
  stage <- data.frame(method     = tab$method,
                      stage      = tab$stage,
                      events     = inp$events[tab$stage],
                      info_frac  = info_frac[tab$stage],
                      events_r1  = n_ev_mean[1L, tab$stage],
                      upper      = bnd$upper[tab$stage],
                      lower      = fut$lower[tab$stage],
                      gamma      = tab$gamma,
                      p_eff      = tab$p_eff,
                      p_fut      = tab$p_fut,
                      p_joint    = tab$p_joint,
                      rcp        = rcp,
                      rcp_uncond = tab$rcp_uncond,
                      se_rcp     = se_rcp,
                      stringsAsFactors = FALSE)
  rownames(stage) <- NULL

  overall <- do.call(rbind, lapply(method, function(m) {
    sub   <- stage[stage$method == m, , drop = FALSE]
    power <- sum(sub$p_eff)
    pj    <- sum(sub$p_joint)
    data.frame(method  = m,
               power   = power,
               p_joint = pj,
               rcp     = if (power > 0) pj / power else NA_real_,
               stringsAsFactors = FALSE)
  }))

  analysis <- data.frame(stage          = seq_len(K),
                         events         = inp$events,
                         time_mean      = colMeans(t_anal),
                         time_sd        = apply(t_anal, 2L, stats::sd),
                         events_r1_mean = n_ev_mean[1L, ],
                         p_no_est_r1    = p_no_est[1L, ],
                         p_inf_r1       = p_inf[1L, ])

  structure(
    list(stage        = stage,
         overall      = overall,
         analysis     = analysis,
         fs_mean      = fs_mean,
         p_no_est     = p_no_est,
         p_inf        = p_inf,
         design       = list(K            = K,
                             upper        = bnd$upper,
                             lower        = fut$lower,
                             lower_design = bnd$lower,
                             binding      = bnd$binding,
                             futility     = fut$label,
                             alpha        = bnd$alpha,
                             timing       = bnd$timing,
                             source       = bnd$source),
         events       = inp$events,
         hr           = inp$hr,
         lambda_c     = inp$lambda_c,
         accrual_time = inp$accrual_time,
         n_accrual    = inp$n_accrual,
         r            = inp$r,
         method       = method,
         scale        = scale,
         gamma_M1     = gamma_M1,
         gamma_M2     = gamma_M2,
         overall_est  = overall_est,
         nsim         = nsim),
    class = c("rcpgsd_sim", "rcpgsd")
  )
}
