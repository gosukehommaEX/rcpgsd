#' Validate Inputs of the Patient-Level Simulation
#'
#' @description
#' Internal helper that validates the numbers of events, the accrual
#' specification, the event and drop-out hazards and the allocation ratio
#' for \code{\link{rcpgsd_sim}} and \code{\link{rcpgsd_events}}, and
#' returns them in a standardized form.
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
#' @param lambda_d Numeric scalar, vector of length 2, or \code{S}-by-2
#'   matrix of drop-out hazards (see \code{\link{rcpgsd_sim}}).
#'
#' @return A list with elements \code{events}, \code{hr},
#'   \code{lambda_c} (vectors of length \code{S}), \code{accrual_time},
#'   \code{n_accrual} (an \code{S}-by-\code{J} matrix), \code{lambda_d}
#'   (an \code{S}-by-2 matrix with the control group in the first column)
#'   and \code{r}.
#'
#' @noRd
.rcpgsd_check_sim_inputs <- function(events, hr, lambda_c, accrual_time,
                                     n_accrual, r, K, lambda_d = 0) {
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

  # Drop-out hazards as an S-by-2 matrix (control, experimental)
  msg_d <- sprintf(paste0("'lambda_d' must be a non-negative numeric ",
                          "scalar, a vector of length 2 (control, ",
                          "experimental), or a %d-by-2 matrix."), S)
  if (!is.numeric(lambda_d) || length(lambda_d) == 0L ||
      any(!is.finite(lambda_d)) || any(lambda_d < 0)) {
    stop(msg_d, call. = FALSE)
  }
  if (is.matrix(lambda_d)) {
    if (nrow(lambda_d) != S || ncol(lambda_d) != 2L) {
      stop(msg_d, call. = FALSE)
    }
    lambda_d <- matrix(as.numeric(lambda_d), S, 2L)
  } else if (length(lambda_d) == 1L) {
    lambda_d <- matrix(as.numeric(lambda_d), S, 2L)
  } else if (length(lambda_d) == 2L) {
    lambda_d <- matrix(as.numeric(lambda_d), S, 2L, byrow = TRUE)
  } else {
    stop(msg_d, call. = FALSE)
  }

  hr_s <- rep_len(as.numeric(hr), S)
  lc_s <- rep_len(as.numeric(lambda_c), S)

  # With drop-out, only a fraction lambda / (lambda + lambda_d) of the
  # patients of each group eventually has an event
  if (any(lambda_d > 0)) {
    n_s  <- rowSums(n_accrual)
    n_t  <- round(n_s * r / (1 + r))
    l_t  <- lc_s * hr_s
    e_max <- sum((n_s - n_t) * lc_s / (lc_s + lambda_d[, 1L]) +
                   n_t * l_t / (l_t + lambda_d[, 2L]))
    if (events[K] >= e_max) {
      stop(sprintf(paste0("The number of events at the final analysis ",
                          "must be smaller than the expected total number ",
                          "of events with drop-out (%.1f)."), e_max),
           call. = FALSE)
    }
  }

  list(events       = as.numeric(events),
       hr           = hr_s,
       lambda_c     = lc_s,
       accrual_time = as.numeric(accrual_time),
       n_accrual    = n_accrual,
       lambda_d     = lambda_d,
       r            = as.numeric(r))
}
