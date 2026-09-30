#' Gauss-Legendre Quadrature Nodes and Weights
#'
#' @description
#' Internal helper computing Gauss-Legendre nodes and weights on a finite
#' interval by the Golub-Welsch algorithm (eigen-decomposition of the
#' symmetric Jacobi matrix of the Legendre polynomials).
#'
#' @param n Positive integer. Number of nodes.
#' @param a Numeric scalar. Lower limit of the interval.
#' @param b Numeric scalar. Upper limit of the interval.
#'
#' @return A list with elements \code{x} (nodes) and \code{w} (weights),
#'   such that \code{sum(w * f(x))} approximates the integral of \code{f}
#'   over \eqn{(a, b)}.
#'
#' @noRd
.rcpgsd_gauss_legendre <- function(n, a = -1, b = 1) {
  n <- as.integer(n)
  if (n == 1L) {
    return(list(x = (a + b) / 2, w = b - a))
  }
  i    <- seq_len(n - 1L)
  beta <- i / sqrt(4 * i^2 - 1)
  J    <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- beta
  J[cbind(i + 1L, i)] <- beta
  eig  <- eigen(J, symmetric = TRUE)
  t    <- eig$values
  w    <- 2 * eig$vectors[1L, ]^2
  ord  <- order(t)
  list(x = (b - a) / 2 * t[ord] + (a + b) / 2,
       w = (b - a) / 2 * w[ord])
}
