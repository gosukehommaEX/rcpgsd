#' Mean and Covariance of Linear Combinations of Regional Estimates
#'
#' @description
#' Internal helper returning the mean vector and covariance matrix of
#' \eqn{L \hat{\delta}}, where \eqn{\hat{\delta}} is the stacked vector of
#' regional estimates defined in \code{.rcpgsd_setup()}.
#'
#' @param st List returned by \code{.rcpgsd_setup()}.
#' @param L Numeric matrix with \code{S * K} columns.
#'
#' @return A list with elements \code{mean} and \code{sigma}.
#'
#' @noRd
.rcpgsd_lin <- function(st, L) {
  L <- matrix(L, ncol = length(st$mu))
  list(mean  = as.numeric(L %*% st$mu),
       sigma = L %*% st$sigma %*% t(L))
}
