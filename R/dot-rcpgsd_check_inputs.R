#' Validate and Standardize Trial Inputs
#'
#' @description
#' Internal helper that validates the numbers of events, the regional event
#' proportions, the true hazard ratios, and the allocation ratio, and returns
#' them in a standardized form.
#'
#' @param events Numeric vector of cumulative numbers of events at the
#'   \code{K} analyses.
#' @param fs Numeric vector of length \code{S} or numeric \code{S}-by-\code{K}
#'   matrix of regional event proportions.
#' @param hr Numeric scalar or vector of length \code{S} of true hazard
#'   ratios. \code{NULL} is allowed and returns \code{NULL}.
#' @param r Numeric scalar. Allocation ratio.
#' @param K Integer. Number of analyses of the design.
#'
#' @return A list with elements \code{events}, \code{fs} (an
#'   \code{S}-by-\code{K} matrix), \code{hr} (a vector of length \code{S}, or
#'   \code{NULL}) and \code{r}.
#'
#' @noRd
.rcpgsd_check_inputs <- function(events, fs, hr, r, K) {
  if (!is.numeric(events) || length(events) != K || any(!is.finite(events)) ||
      any(events <= 0)) {
    stop(sprintf("'events' must be a positive numeric vector of length %d ",
                 K), "(one element per analysis of 'design').", call. = FALSE)
  }
  if (K > 1L && any(diff(events) <= 0)) {
    stop("'events' must be strictly increasing.", call. = FALSE)
  }

  if (!is.numeric(fs) || any(!is.finite(fs))) {
    stop("'fs' must be a numeric vector or matrix.", call. = FALSE)
  }
  if (is.matrix(fs)) {
    if (ncol(fs) != K) {
      stop(sprintf("A matrix 'fs' must have %d columns (one per analysis).",
                   K), call. = FALSE)
    }
    fs_mat <- fs
  } else {
    fs_mat <- matrix(fs, nrow = length(fs), ncol = K)
  }
  S <- nrow(fs_mat)
  if (S < 2L) {
    stop("'fs' must describe at least two regions.", call. = FALSE)
  }
  if (any(fs_mat <= 0)) {
    stop("All elements of 'fs' must be positive.", call. = FALSE)
  }
  if (any(abs(colSums(fs_mat) - 1) > 1e-8)) {
    stop("The regional proportions in 'fs' must sum to 1 at each analysis.",
         call. = FALSE)
  }
  events_region <- sweep(fs_mat, 2L, events, "*")
  if (K > 1L && any(events_region[, -1L, drop = FALSE] -
                    events_region[, -K, drop = FALSE] < -1e-8)) {
    stop("The implied regional numbers of events 'fs * events' must be ",
         "non-decreasing across analyses in every region.", call. = FALSE)
  }

  hr_vec <- NULL
  if (!is.null(hr)) {
    if (!is.numeric(hr) || any(!is.finite(hr)) || any(hr <= 0) ||
        !(length(hr) %in% c(1L, S))) {
      stop(sprintf(paste0("'hr' must be a positive numeric scalar or a ",
                          "vector of length %d."), S), call. = FALSE)
    }
    hr_vec <- rep_len(as.numeric(hr), S)
  }

  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= 0) {
    stop("'r' must be a single positive number.", call. = FALSE)
  }

  list(events = as.numeric(events),
       fs     = fs_mat,
       hr     = hr_vec,
       r      = as.numeric(r))
}
