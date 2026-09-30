#' Extract Stopping Boundaries from a Group Sequential Design
#'
#' @description
#' Internal helper that converts a group sequential design into a common
#' representation on the standardized Z-scale. Supported inputs are a
#' \code{gsDesign} object (package \pkg{gsDesign}), a
#' \code{TrialDesignGroupSequential} or \code{TrialDesignInverseNormal}
#' object (package \pkg{rpact}), a list with elements \code{upper} and
#' optionally \code{lower}, \code{binding}, \code{alpha} and \code{timing},
#' or a numeric vector of efficacy boundaries.
#'
#' Efficacy is assumed to correspond to large values of the Z-statistic.
#' Lower boundaries at or below \code{-6} (the value used by \pkg{rpact}
#' for "no futility boundary") are treated as \code{-Inf}. The
#' lower boundary at the final analysis is always set to \code{-Inf},
#' because only efficacy is evaluated at the final analysis.
#'
#' @param design A design object or list (see Description).
#'
#' @return A list with elements \code{K} (number of analyses), \code{upper}
#'   (efficacy boundaries), \code{lower} (futility boundaries, length
#'   \code{K}), \code{binding} (logical), \code{has_futility} (logical),
#'   \code{alpha} (one-sided significance level, or \code{NA}),
#'   \code{timing} (planned information fractions, or \code{NULL}) and
#'   \code{source} (character label of the input type).
#'
#' @noRd
.rcpgsd_bounds <- function(design) {
  if (inherits(design, "gsDesign")) {
    upper     <- as.numeric(design$upper$bound)
    K         <- length(upper)
    test_type <- design$test.type
    if (is.null(test_type) || test_type == 1 || is.null(design$lower$bound)) {
      lower <- rep(-Inf, K)
    } else {
      lower <- as.numeric(design$lower$bound)
    }
    binding <- !is.null(test_type) && test_type %in% c(2, 3, 5)
    alpha   <- if (is.null(design$alpha)) NA_real_ else as.numeric(design$alpha)
    timing  <- if (is.null(design$timing)) NULL else as.numeric(design$timing)
    source  <- "gsDesign"
  } else if (inherits(design, c("TrialDesignGroupSequential",
                                "TrialDesignInverseNormal"))) {
    upper <- as.numeric(design$criticalValues)
    K     <- length(upper)
    sided <- design$sided
    if (!is.null(sided) && sided == 2) {
      lower <- -upper
      alpha <- as.numeric(design$alpha) / 2
      binding <- TRUE
    } else {
      fut <- design$futilityBounds
      if (is.null(fut) || K == 1L) {
        lower <- rep(-Inf, K)
      } else {
        lower <- c(as.numeric(fut), -Inf)
      }
      alpha   <- as.numeric(design$alpha)
      binding <- isTRUE(design$bindingFutility)
    }
    timing <- as.numeric(design$informationRates)
    source <- "rpact"
  } else if (inherits(design, "TrialDesign")) {
    stop("Only rpact designs of class 'TrialDesignGroupSequential' or ",
         "'TrialDesignInverseNormal' are supported.", call. = FALSE)
  } else if (is.list(design)) {
    if (is.null(design$upper)) {
      stop("A list 'design' must contain an element 'upper'.", call. = FALSE)
    }
    upper <- as.numeric(design$upper)
    K     <- length(upper)
    if (is.null(design$lower)) {
      lower <- rep(-Inf, K)
    } else {
      lower <- as.numeric(design$lower)
      if (length(lower) == K - 1L) lower <- c(lower, -Inf)
      if (length(lower) != K) {
        stop("'design$lower' must have length K or K - 1.", call. = FALSE)
      }
    }
    binding <- if (is.null(design$binding)) FALSE else isTRUE(design$binding)
    alpha   <- if (is.null(design$alpha)) NA_real_ else as.numeric(design$alpha)
    timing  <- if (is.null(design$timing)) NULL else as.numeric(design$timing)
    source  <- "list"
  } else if (is.numeric(design)) {
    upper   <- as.numeric(design)
    K       <- length(upper)
    lower   <- rep(-Inf, K)
    binding <- FALSE
    alpha   <- NA_real_
    timing  <- NULL
    source  <- "numeric"
  } else {
    stop("'design' must be a gsDesign object, an rpact design object, ",
         "a list with element 'upper', or a numeric vector.", call. = FALSE)
  }

  if (K < 1L || any(is.na(upper))) {
    stop("The efficacy boundaries must be non-missing.", call. = FALSE)
  }
  lower[is.na(lower) | lower <= -6] <- -Inf
  lower[K] <- -Inf
  if (K > 1L && any(lower[-K] >= upper[-K])) {
    stop("Each futility boundary must be smaller than the corresponding ",
         "efficacy boundary.", call. = FALSE)
  }
  has_futility <- any(is.finite(lower))
  if (!has_futility) binding <- FALSE
  if (!is.null(timing) && length(timing) != K) timing <- NULL

  list(K            = K,
       upper        = upper,
       lower        = lower,
       binding      = binding,
       has_futility = has_futility,
       alpha        = alpha,
       timing       = timing,
       source       = source)
}
