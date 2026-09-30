#' Futility Boundaries Used in the Calculation
#'
#' @description
#' Internal helper that determines the futility boundaries used in the
#' probability calculation. Binding futility boundaries are always applied.
#' Non-binding futility boundaries are applied when
#' \code{futility = "obey"} and replaced by \code{-Inf} when
#' \code{futility = "ignore"}.
#'
#' @param bnd List returned by \code{.rcpgsd_bounds()}.
#' @param futility Character. \code{"obey"} or \code{"ignore"}.
#'
#' @return A list with elements \code{lower} (numeric vector of length
#'   \code{K}) and \code{label} (character description).
#'
#' @noRd
.rcpgsd_futility <- function(bnd, futility) {
  if (!bnd$has_futility) {
    return(list(lower = bnd$lower, label = "none"))
  }
  if (bnd$binding) {
    if (futility == "ignore") {
      warning("The futility boundaries of 'design' are binding and are ",
              "therefore obeyed; futility = \"ignore\" was not applied.",
              call. = FALSE)
    }
    return(list(lower = bnd$lower, label = "binding"))
  }
  if (futility == "obey") {
    list(lower = bnd$lower, label = "non-binding (obeyed)")
  } else {
    list(lower = rep(-Inf, bnd$K), label = "non-binding (ignored)")
  }
}
