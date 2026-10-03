#' Multivariate Normal Rectangle Probability
#'
#' @description
#' Internal wrapper around \code{mvtnorm::pmvnorm()} that returns
#' \eqn{\Pr(a < X < b)} for \eqn{X \sim N(\mu, \Sigma)}. The problem is
#' standardized to a correlation matrix before integration. One-dimensional
#' problems are evaluated with \code{stats::pnorm()} using the upper tail
#' when this is numerically more accurate. For two or more dimensions the
#' Genz-Bretz algorithm is used with a fixed seed, so that repeated calls
#' return identical values, and with absolute and relative error
#' tolerances of 1e-10 and 1e-5, so that small probabilities (e.g., early
#' stopping under the null hypothesis) are evaluated accurately. Singular
#' covariance matrices, which arise for
#' Method 2 because the overall statistic is a linear combination of all
#' regional estimates, are handled by \code{mvtnorm::pmvnorm()}.
#'
#' @param lower Numeric vector of lower limits.
#' @param upper Numeric vector of upper limits.
#' @param mean Numeric mean vector.
#' @param sigma Numeric covariance matrix.
#'
#' @return Numeric scalar in \eqn{[0, 1]}.
#'
#' @importFrom mvtnorm pmvnorm GenzBretz
#' @importFrom stats pnorm
#'
#' @noRd
.rcpgsd_pmvn <- function(lower, upper, mean, sigma) {
  mean  <- as.numeric(mean)
  sigma <- as.matrix(sigma)
  if (any(lower >= upper)) {
    return(0)
  }
  sds <- sqrt(diag(sigma))
  lo  <- (lower - mean) / sds
  up  <- (upper - mean) / sds
  if (length(mean) == 1L) {
    if (lo > 0) {
      p <- stats::pnorm(lo, lower.tail = FALSE) -
        stats::pnorm(up, lower.tail = FALSE)
    } else {
      p <- stats::pnorm(up) - stats::pnorm(lo)
    }
    return(min(1, max(0, p)))
  }
  corr <- sigma / outer(sds, sds)
  corr <- (corr + t(corr)) / 2
  diag(corr) <- 1
  p <- mvtnorm::pmvnorm(lower     = lo,
                        upper     = up,
                        corr      = corr,
                        algorithm = mvtnorm::GenzBretz(maxpts = 500000,
                                                       abseps = 1e-10,
                                                       releps = 1e-5),
                        keepAttr  = FALSE,
                        seed      = 1L)
  min(1, max(0, as.numeric(p)))
}
