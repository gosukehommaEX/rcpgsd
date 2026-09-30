#' Stagewise Probabilities by the Closed-Form Approach
#'
#' @description
#' Internal helper that evaluates, for each consistency method and each
#' analysis, the probability of stopping for efficacy, the probability of
#' stopping for futility, the joint probability of efficacy stopping and
#' regional consistency, and the unconditional regional consistency
#' probability, using multivariate normal probabilities.
#'
#' @param st List returned by \code{.rcpgsd_setup()}.
#' @param upper,lower Numeric vectors of efficacy and futility boundaries
#'   used in the calculation.
#' @param method Character vector, a subset of \code{c("M1", "M2")}.
#' @param scale Character. \code{"loghr"} or \code{"hr"}.
#' @param gamma_M1,gamma_M2 Numeric vectors of length \code{K}.
#' @param hr_method Character. \code{"exact"} or \code{"delta"}.
#' @param n_nodes Integer. Number of Gauss-Legendre nodes.
#'
#' @return A data frame with one row per method and analysis and columns
#'   \code{method}, \code{stage}, \code{gamma}, \code{p_eff}, \code{p_fut},
#'   \code{p_joint} and \code{rcp_uncond}.
#'
#' @noRd
.rcpgsd_formula <- function(st, upper, lower, method, scale,
                            gamma_M1, gamma_M2, hr_method, n_nodes) {
  K     <- st$K
  p_eff <- vapply(seq_len(K), function(k) {
    .rcpgsd_prob_stop(st, upper, lower, k, "efficacy")
  }, numeric(1L))
  p_fut <- vapply(seq_len(K), function(k) {
    if (k == K) 0 else .rcpgsd_prob_stop(st, upper, lower, k, "futility")
  }, numeric(1L))

  out <- lapply(method, function(m) {
    gam <- if (m == "M1") gamma_M1 else gamma_M2
    p_joint <- vapply(seq_len(K), function(k) {
      if (p_eff[k] <= 0) {
        return(0)
      }
      .rcpgsd_prob_joint(st, upper, lower, k, m, scale, gam[k],
                         hr_method, n_nodes)
    }, numeric(1L))
    p_uncond <- vapply(seq_len(K), function(k) {
      .rcpgsd_prob_uncond(st, k, m, scale, gam[k], hr_method)
    }, numeric(1L))
    data.frame(method     = m,
               stage      = seq_len(K),
               gamma      = gam,
               p_eff      = p_eff,
               p_fut      = p_fut,
               p_joint    = pmin(p_joint, p_eff),
               rcp_uncond = p_uncond,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
