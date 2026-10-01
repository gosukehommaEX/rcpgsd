#' Validate Inputs of the Patient-Level Simulation
#'
#' @description
#' Internal helper that validates the numbers of events, the accrual
#' specification, the hazards and the allocation ratio for
#' \code{\link{rcpgsd_sim}}, and returns them in a standardized form.
#'
#' @param events Numeric vector of cumulative numbers of events (length
#'   \code{K}).
#' @param hr Numeric scalar or vector of true hazard ratios.
#' @param lambda_c Numeric scalar or vector of control hazards.
#' @param accrual_time Numeric vector of accrual interval boundaries.
#' @param n_accrual Numeric matrix (or vector when there is one accrual
#'   interval) of numbers of patients by region and accrual interval.
#' @param r Numeric scalar. Allocation ratio.
#' @param K Integer. Number of analyses of the design.
#'
#' @return A list with elements \code{events}, \code{hr},
#'   \code{lambda_c} (vectors of length \code{S}), \code{accrual_time},
#'   \code{n_accrual} (an \code{S}-by-\code{J} matrix) and \code{r}.
#'
#' @noRd
.rcpgsd_check_sim_inputs <- function(events, hr, lambda_c, accrual_time,
                                     n_accrual, r, K) {
  if (!is.numeric(events) || length(events) != K || any(!is.finite(events)) ||
      any(events < 1) || any(events != round(events))) {
    stop(sprintf(paste0("'events' must be a vector of %d positive integers ",
                        "(one element per analysis of 'design')."), K),
         call. = FALSE)
  }
  if (K > 1L && any(diff(events) <= 0)) {
    stop("'events' must be strictly increasing.", call. = FALSE)
  }

  if (!is.numeric(accrual_time) || length(accrual_time) < 2L ||
      any(!is.finite(accrual_time)) || accrual_time[1L] != 0 ||
      any(diff(accrual_time) <= 0)) {
    stop("'accrual_time' must be a strictly increasing numeric vector ",
         "starting at 0.", call. = FALSE)
  }
  J <- length(accrual_time) - 1L

  if (!is.numeric(n_accrual) || any(!is.finite(n_accrual))) {
    stop("'n_accrual' must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(n_accrual)) {
    if (J != 1L) {
      stop("'n_accrual' must be a matrix with one column per accrual ",
           "interval when there are several accrual intervals.",
           call. = FALSE)
    }
    n_accrual <- matrix(n_accrual, ncol = 1L)
  }
  if (ncol(n_accrual) != J) {
    stop(sprintf(paste0("'n_accrual' must have %d columns (one per accrual ",
                        "interval defined by 'accrual_time')."), J),
         call. = FALSE)
  }
  S <- nrow(n_accrual)
  if (S < 2L) {
    stop("'n_accrual' must describe at least two regions.", call. = FALSE)
  }
  if (any(n_accrual < 0) || any(n_accrual != round(n_accrual))) {
    stop("'n_accrual' must contain non-negative integers.", call. = FALSE)
  }
  if (any(rowSums(n_accrual) == 0)) {
    stop("Every region must enrol at least one patient.", call. = FALSE)
  }
  if (sum(n_accrual) < events[K]) {
    stop("The total number of patients must be at least the number of ",
         "events at the final analysis.", call. = FALSE)
  }

  if (!is.numeric(hr) || any(!is.finite(hr)) || any(hr <= 0) ||
      !(length(hr) %in% c(1L, S))) {
    stop(sprintf(paste0("'hr' must be a positive numeric scalar or a ",
                        "vector of length %d."), S), call. = FALSE)
  }
  if (!is.numeric(lambda_c) || any(!is.finite(lambda_c)) ||
      any(lambda_c <= 0) || !(length(lambda_c) %in% c(1L, S))) {
    stop(sprintf(paste0("'lambda_c' must be a positive numeric scalar or a ",
                        "vector of length %d."), S), call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r <= 0) {
    stop("'r' must be a single positive number.", call. = FALSE)
  }

  list(events       = as.numeric(events),
       hr           = rep_len(as.numeric(hr), S),
       lambda_c     = rep_len(as.numeric(lambda_c), S),
       accrual_time = as.numeric(accrual_time),
       n_accrual    = n_accrual,
       r            = as.numeric(r))
}
