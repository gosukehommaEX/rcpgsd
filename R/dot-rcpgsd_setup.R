#' Joint Distribution of Regional Log Hazard Ratio Estimates
#'
#' @description
#' Internal helper that builds the canonical joint distribution of the
#' regional treatment effect estimates across all regions and analyses.
#'
#' The treatment effect is \eqn{\delta_s = -\log(\mathrm{HR}_s)}, so that a
#' positive value indicates a benefit. With \eqn{d_{s,k}} events in region
#' \eqn{s} at analysis \eqn{k} and allocation ratio \eqn{r}, the regional
#' Fisher information is \eqn{I_{s,k} = r d_{s,k} / (1 + r)^2}. The regional
#' estimates satisfy
#' \eqn{\hat{\delta}_{s,k} \sim N(\delta_s, 1 / I_{s,k})} with
#' \eqn{\mathrm{Cov}(\hat{\delta}_{s,k}, \hat{\delta}_{s,k'}) =
#' 1 / I_{s,\max(k,k')}} (independent increments) and independence across
#' regions. The overall estimate is
#' \eqn{\hat{\delta}_k = \sum_s w_{s,k} \hat{\delta}_{s,k}} with
#' \eqn{w_{s,k} = I_{s,k} / I_k} and \eqn{I_k = \sum_s I_{s,k}}, and the
#' overall Wald statistic is \eqn{Z_k = \hat{\delta}_k \sqrt{I_k}}.
#'
#' The stacked vector of regional estimates is ordered by region and then by
#' analysis, that is, element \eqn{(s - 1) K + k} corresponds to
#' \eqn{\hat{\delta}_{s,k}}.
#'
#' @param events Numeric vector of cumulative total events (length \code{K}).
#' @param fs Numeric \code{S}-by-\code{K} matrix of regional event
#'   proportions.
#' @param hr Numeric vector of true hazard ratios (length \code{S}).
#' @param r Numeric scalar. Allocation ratio (experimental to control).
#'
#' @return A list with elements \code{K}, \code{S}, \code{delta},
#'   \code{info_region} (\code{S}-by-\code{K}), \code{info_total}
#'   (length \code{K}), \code{weight} (\code{S}-by-\code{K}), \code{mu}
#'   and \code{sigma} (mean and covariance of the stacked regional
#'   estimates), \code{L_Z} and \code{L_V} (\code{K}-row matrices mapping the
#'   stacked vector to \eqn{Z_k} and \eqn{\hat{\delta}_k}).
#'
#' @noRd
.rcpgsd_setup <- function(events, fs, hr, r) {
  K <- length(events)
  S <- nrow(fs)
  info_region <- sweep(fs, 2L, events, "*") * r / (1 + r)^2
  info_total  <- colSums(info_region)
  weight      <- sweep(info_region, 2L, info_total, "/")
  delta       <- -log(hr)

  n     <- S * K
  mu    <- rep(delta, each = K)
  sigma <- matrix(0, n, n)
  idx_max <- outer(seq_len(K), seq_len(K), pmax)
  for (s in seq_len(S)) {
    pos <- (s - 1L) * K + seq_len(K)
    sigma[pos, pos] <- matrix(1 / info_region[s, idx_max], K, K)
  }

  L_V <- matrix(0, K, n)
  for (k in seq_len(K)) {
    for (s in seq_len(S)) {
      L_V[k, (s - 1L) * K + k] <- weight[s, k]
    }
  }
  L_Z <- L_V * sqrt(info_total)

  list(K           = K,
       S           = S,
       delta       = delta,
       info_region = info_region,
       info_total  = info_total,
       weight      = weight,
       mu          = mu,
       sigma       = sigma,
       L_Z         = L_Z,
       L_V         = L_V)
}
