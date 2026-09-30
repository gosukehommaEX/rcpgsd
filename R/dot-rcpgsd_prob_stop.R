#' Stagewise Stopping Probability
#'
#' @description
#' Internal helper returning the probability that the trial stops at
#' analysis \code{k} for efficacy, that is,
#' \eqn{\Pr(l_j < Z_j < u_j, j < k; Z_k \ge u_k)}, or for futility, that
#' is, \eqn{\Pr(l_j < Z_j < u_j, j < k; Z_k \le l_k)}.
#'
#' @param st List returned by \code{.rcpgsd_setup()}.
#' @param upper Numeric vector of efficacy boundaries (length \code{K}).
#' @param lower Numeric vector of futility boundaries (length \code{K}); use
#'   \code{-Inf} where no futility stopping applies.
#' @param k Integer. Analysis index.
#' @param type Character. \code{"efficacy"} or \code{"futility"}.
#'
#' @return Numeric scalar.
#'
#' @noRd
.rcpgsd_prob_stop <- function(st, upper, lower, k,
                              type = c("efficacy", "futility")) {
  type <- match.arg(type)
  prev <- seq_len(k - 1L)
  if (type == "efficacy") {
    lo <- c(lower[prev], upper[k])
    up <- c(upper[prev], Inf)
  } else {
    if (!is.finite(lower[k])) {
      return(0)
    }
    lo <- c(lower[prev], -Inf)
    up <- c(upper[prev], lower[k])
  }
  mv <- .rcpgsd_lin(st, st$L_Z[seq_len(k), , drop = FALSE])
  .rcpgsd_pmvn(lo, up, mv$mean, mv$sigma)
}
