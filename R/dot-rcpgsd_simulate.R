#' Stagewise Probabilities by Monte Carlo Simulation
#'
#' @description
#' Internal helper that simulates the canonical joint distribution of the
#' regional log hazard ratio estimates (independent increments within each
#' region, independence across regions) and estimates the same quantities
#' as \code{.rcpgsd_formula()}. The consistency criteria are applied in
#' their original (nonlinear) form, so the simulation can be used to verify
#' both the exact and the delta-method results on the hazard ratio scale.
#'
#' @param st List returned by \code{.rcpgsd_setup()}.
#' @param upper,lower Numeric vectors of efficacy and futility boundaries
#'   used in the calculation.
#' @param method Character vector, a subset of \code{c("M1", "M2")}.
#' @param scale Character. \code{"loghr"} or \code{"hr"}.
#' @param gamma_M1,gamma_M2 Numeric vectors of length \code{K}.
#' @param nsim Positive integer. Number of simulated trials.
#' @param seed Integer or \code{NULL}. Random seed.
#'
#' @return A data frame with the same columns as \code{.rcpgsd_formula()}.
#'
#' @importFrom stats rnorm
#'
#' @noRd
.rcpgsd_simulate <- function(st, upper, lower, method, scale,
                             gamma_M1, gamma_M2, nsim, seed) {
  if (!is.null(seed)) set.seed(seed)
  K <- st$K
  S <- st$S

  score     <- matrix(0, nsim, S)
  alive     <- rep(TRUE, nsim)
  info_prev <- rep(0, S)
  res <- vector("list", K)

  for (k in seq_len(K)) {
    # Independent score increments in each region
    d_info <- st$info_region[, k] - info_prev
    for (s in seq_len(S)) {
      if (d_info[s] > 0) {
        score[, s] <- score[, s] +
          stats::rnorm(nsim, st$delta[s] * d_info[s], sqrt(d_info[s]))
      }
    }
    info_prev <- st$info_region[, k]
    est   <- sweep(score, 2L, st$info_region[, k], "/")
    v_hat <- as.numeric(est %*% st$weight[, k])
    z_k   <- v_hat * sqrt(st$info_total[k])

    eff <- alive & z_k >= upper[k]
    fut <- if (k < K) alive & z_k <= lower[k] else rep(FALSE, nsim)

    cons <- list()
    if ("M1" %in% method) {
      g1 <- gamma_M1[k]
      cons$M1 <- if (scale == "loghr") {
        est[, 1L] > g1 * v_hat
      } else {
        (1 - exp(-est[, 1L])) > g1 * (1 - exp(-v_hat))
      }
    }
    if ("M2" %in% method) {
      thr <- if (scale == "loghr") gamma_M2[k] else -log(1 - gamma_M2[k])
      cons$M2 <- rowSums(est > thr) == S
    }
    res[[k]] <- do.call(rbind, lapply(method, function(m) {
      data.frame(method     = m,
                 stage      = k,
                 gamma      = if (m == "M1") gamma_M1[k] else gamma_M2[k],
                 p_eff      = mean(eff),
                 p_fut      = mean(fut),
                 p_joint    = mean(eff & cons[[m]]),
                 rcp_uncond = mean(cons[[m]]),
                 stringsAsFactors = FALSE)
    }))
    alive <- alive & !eff & !fut
  }
  out <- do.call(rbind, res)
  out <- out[order(match(out$method, method), out$stage), , drop = FALSE]
  rownames(out) <- NULL
  out
}
