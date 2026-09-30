#' Unconditional Regional Consistency Probability at an Analysis
#'
#' @description
#' Internal helper returning \eqn{\Pr(C_k)}, the probability that the
#' regional consistency criterion is met with the data available at
#' analysis \eqn{k}, irrespective of the stopping boundaries. For a
#' single-stage design this is the unconditional regional consistency
#' probability of Kawai et al. (2008), Ikeda and Bretz (2010) and, for the
#' hazard ratio scale, Hayashi and Itoh (2018). The criteria are those of
#' \code{.rcpgsd_prob_joint()}.
#'
#' For Method 1 with \code{scale = "hr"} and \code{hr_method = "exact"},
#' the probability is evaluated with \code{stats::integrate()} as
#' \eqn{\Pr(\hat{\delta}_{1,k} \ge x_u) + \int_{-\infty}^{x_u}
#' \phi(x) \Pr(\hat{\delta}_k < V(x) \mid \hat{\delta}_{1,k} = x) dx}.
#'
#' @param st List returned by \code{.rcpgsd_setup()}.
#' @param k Integer. Analysis index.
#' @param method Character. \code{"M1"} or \code{"M2"}.
#' @param scale Character. \code{"loghr"} or \code{"hr"}.
#' @param gamma Numeric scalar. Consistency threshold.
#' @param hr_method Character. \code{"exact"} or \code{"delta"}.
#'
#' @return Numeric scalar.
#'
#' @importFrom stats dnorm pnorm integrate
#'
#' @noRd
.rcpgsd_prob_uncond <- function(st, k, method, scale, gamma,
                                hr_method = "exact") {
  K <- st$K
  S <- st$S
  row_A <- function(s) {
    a <- numeric(S * K)
    a[(s - 1L) * K + k] <- 1
    a
  }

  if (method == "M2") {
    thr <- if (scale == "loghr") gamma else -log(1 - gamma)
    mv  <- .rcpgsd_lin(st, t(vapply(seq_len(S), row_A, numeric(S * K))))
    return(.rcpgsd_pmvn(rep(thr, S), rep(Inf, S), mv$mean, mv$sigma))
  }

  A1 <- row_A(1L)
  Vk <- st$L_V[k, ]
  if (scale == "loghr" || hr_method == "delta") {
    if (scale == "loghr") {
      D      <- A1 - gamma * Vk
      offset <- 0
    } else {
      mu_V   <- sum(Vk * st$mu)
      denom  <- 1 - gamma + gamma * exp(-mu_V)
      a      <- gamma * exp(-mu_V) / denom
      D      <- A1 - a * Vk
      offset <- -log(denom) - a * mu_V
    }
    mv <- .rcpgsd_lin(st, D)
    return(.rcpgsd_pmvn(offset, Inf, mv$mean, mv$sigma))
  }

  x_u <- -log(1 - gamma)
  mv  <- .rcpgsd_lin(st, rbind(A1, Vk))
  m_X  <- mv$mean[1L]
  sd_X <- sqrt(mv$sigma[1L, 1L])
  m_V  <- mv$mean[2L]
  c_XV <- mv$sigma[1L, 2L]
  sd_V_cond <- sqrt(mv$sigma[2L, 2L] - c_XV^2 / mv$sigma[1L, 1L])
  p1 <- stats::pnorm(x_u, m_X, sd_X, lower.tail = FALSE)
  if (gamma == 0) {
    return(p1)
  }
  integrand <- function(x) {
    arg <- exp(-x) - 1 + gamma
    v_x <- ifelse(arg <= 0, Inf, -log(pmax(arg, 1e-300) / gamma))
    m_c <- m_V + c_XV / sd_X^2 * (x - m_X)
    stats::dnorm(x, m_X, sd_X) * stats::pnorm((v_x - m_c) / sd_V_cond)
  }
  lower_int <- m_X - 12 * sd_X
  if (lower_int >= x_u) {
    return(p1)
  }
  p2 <- stats::integrate(integrand, lower = lower_int, upper = x_u,
                         rel.tol = 1e-10, abs.tol = 1e-13,
                         subdivisions = 1000L)$value
  min(1, p1 + p2)
}
