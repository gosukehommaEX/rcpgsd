#' Joint Probability of Efficacy Stopping and Regional Consistency
#'
#' @description
#' Internal helper returning
#' \eqn{\Pr(E_k \cap C_k)}, where
#' \eqn{E_k = \{l_j < Z_j < u_j, j < k; Z_k \ge u_k\}} is the event that the
#' trial stops for efficacy at analysis \eqn{k} and \eqn{C_k} is the
#' regional consistency event evaluated with the data available at analysis
#' \eqn{k}.
#'
#' Consistency criteria (region 1 is the region of interest, and
#' \eqn{\hat{\delta} = -\log \widehat{\mathrm{HR}}}):
#' \itemize{
#'   \item Method 1, \code{scale = "loghr"}:
#'     \eqn{\hat{\delta}_{1,k} > \gamma \hat{\delta}_k}.
#'   \item Method 1, \code{scale = "hr"}:
#'     \eqn{1 - \widehat{\mathrm{HR}}_{1,k} >
#'     \gamma (1 - \widehat{\mathrm{HR}}_k)}. With
#'     \code{hr_method = "exact"} the probability is evaluated by conditioning
#'     on \eqn{\hat{\delta}_{1,k}} and integrating numerically
#'     (Gauss-Legendre quadrature) over the part of the range where the
#'     criterion constrains \eqn{\hat{\delta}_k}; the remaining part is a
#'     multivariate normal probability. With \code{hr_method = "delta"} the
#'     criterion is linearized around the mean of \eqn{\hat{\delta}_k} by the
#'     delta method (Hayashi and Itoh, 2018), which gives a multivariate
#'     normal probability.
#'   \item Method 2, \code{scale = "loghr"}:
#'     \eqn{\hat{\delta}_{s,k} > \gamma} for all \eqn{s}.
#'   \item Method 2, \code{scale = "hr"}:
#'     \eqn{1 - \widehat{\mathrm{HR}}_{s,k} > \gamma} for all \eqn{s}.
#' }
#'
#' @param st List returned by \code{.rcpgsd_setup()}.
#' @param upper,lower Numeric vectors of efficacy and futility boundaries.
#' @param k Integer. Analysis index.
#' @param method Character. \code{"M1"} or \code{"M2"}.
#' @param scale Character. \code{"loghr"} or \code{"hr"}.
#' @param gamma Numeric scalar. Consistency threshold used at analysis
#'   \code{k}.
#' @param hr_method Character. \code{"exact"} or \code{"delta"}; used only
#'   for Method 1 with \code{scale = "hr"}.
#' @param n_nodes Integer. Number of Gauss-Legendre nodes.
#'
#' @return Numeric scalar.
#'
#' @importFrom stats dnorm
#'
#' @noRd
.rcpgsd_prob_joint <- function(st, upper, lower, k, method, scale, gamma,
                               hr_method = "exact", n_nodes = 64L) {
  K    <- st$K
  S    <- st$S
  prev <- seq_len(k - 1L)
  z_lo <- c(lower[prev], upper[k])
  z_up <- c(upper[prev], Inf)
  L_Zk <- st$L_Z[seq_len(k), , drop = FALSE]

  # Selection row for the stage-k estimate in region s
  row_A <- function(s) {
    a <- numeric(S * K)
    a[(s - 1L) * K + k] <- 1
    a
  }

  # ---------- Method 2 ----------
  if (method == "M2") {
    thr <- if (scale == "loghr") gamma else -log(1 - gamma)
    L   <- rbind(L_Zk, t(vapply(seq_len(S), row_A, numeric(S * K))))
    mv  <- .rcpgsd_lin(st, L)
    return(.rcpgsd_pmvn(c(z_lo, rep(thr, S)), c(z_up, rep(Inf, S)),
                        mv$mean, mv$sigma))
  }

  # ---------- Method 1: linear criteria ----------
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
    mv <- .rcpgsd_lin(st, rbind(L_Zk, D))
    return(.rcpgsd_pmvn(c(z_lo, offset), c(z_up, Inf), mv$mean, mv$sigma))
  }

  # ---------- Method 1: exact criterion on the hazard ratio scale ----------
  # Consistency holds if hat_delta_1 >= x_u, or if x_s <= hat_delta_1 < x_u
  # and hat_delta_k < V(hat_delta_1), where
  # V(x) = -log((exp(-x) - 1 + gamma) / gamma).
  x_u <- -log(1 - gamma)
  mv  <- .rcpgsd_lin(st, rbind(L_Zk, A1))
  p1  <- .rcpgsd_pmvn(c(z_lo, x_u), c(z_up, Inf), mv$mean, mv$sigma)
  if (gamma == 0 || !is.finite(upper[k])) {
    return(p1)
  }

  sqrt_I <- sqrt(st$info_total[k])
  c_k    <- upper[k]
  x_s    <- -log(1 - gamma + gamma * exp(-c_k / sqrt_I))

  m_Z   <- mv$mean[seq_len(k)]
  C_Z   <- mv$sigma[seq_len(k), seq_len(k), drop = FALSE]
  c_ZX  <- mv$sigma[seq_len(k), k + 1L]
  m_X   <- mv$mean[k + 1L]
  v_X   <- mv$sigma[k + 1L, k + 1L]
  sd_X  <- sqrt(v_X)
  C_cond <- C_Z - tcrossprod(c_ZX) / v_X

  a_int <- max(x_s, m_X - 10 * sd_X)
  b_int <- min(x_u, m_X + 10 * sd_X)
  if (a_int >= b_int) {
    return(p1)
  }
  gl <- .rcpgsd_gauss_legendre(n_nodes, a_int, b_int)
  vals <- vapply(seq_along(gl$x), function(i) {
    x   <- gl$x[i]
    arg <- exp(-x) - 1 + gamma
    u_k <- if (arg <= 0) Inf else -log(arg / gamma) * sqrt_I
    if (u_k <= c_k) {
      return(0)
    }
    up_i <- z_up
    up_i[k] <- u_k
    stats::dnorm(x, m_X, sd_X) *
      .rcpgsd_pmvn(z_lo, up_i, m_Z + c_ZX / v_X * (x - m_X), C_cond)
  }, numeric(1L))
  p1 + sum(gl$w * vals)
}
