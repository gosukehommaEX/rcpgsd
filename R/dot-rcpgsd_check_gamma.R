#' Validate Consistency Thresholds
#'
#' @description
#' Internal helper that validates a consistency threshold supplied as a
#' scalar or as a vector with one element per analysis, and returns it as a
#' vector of length \code{K}.
#'
#' @param gamma Numeric scalar or vector of length \code{K}.
#' @param K Integer. Number of analyses.
#' @param name Character. Argument name used in error messages
#'   (\code{"gamma_M1"} or \code{"gamma_M2"}).
#' @param scale Character. \code{"loghr"} or \code{"hr"}.
#'
#' @return Numeric vector of length \code{K}.
#'
#' @noRd
.rcpgsd_check_gamma <- function(gamma, K, name, scale) {
  if (!is.numeric(gamma) || !(length(gamma) %in% c(1L, K)) ||
      any(!is.finite(gamma))) {
    stop(sprintf("'%s' must be a numeric scalar or a vector of length %d.",
                 name, K), call. = FALSE)
  }
  if (any(gamma < 0)) {
    stop(sprintf("'%s' must be non-negative.", name), call. = FALSE)
  }
  if (name == "gamma_M1" && any(gamma >= 1)) {
    stop("'gamma_M1' must be smaller than 1.", call. = FALSE)
  }
  if (name == "gamma_M2" && scale == "hr" && any(gamma >= 1)) {
    stop("'gamma_M2' must be smaller than 1 when scale = \"hr\".",
         call. = FALSE)
  }
  rep_len(as.numeric(gamma), K)
}
