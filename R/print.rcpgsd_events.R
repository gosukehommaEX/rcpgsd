#' Print Method for rcpgsd_events Objects
#'
#' @description
#' Prints the expected analysis times and the expected regional numbers and
#' proportions of events computed by \code{\link{rcpgsd_events}}.
#'
#' @param x An object of class \code{"rcpgsd_events"}.
#' @param digits Integer. Number of decimal places for proportions. Default
#'   is \code{4}.
#' @param ... Further arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rcpgsd_events}}
#'
#' @examples
#' ev <- rcpgsd_events(events = c(100, 200), hr = 0.7, lambda_c = 0.05,
#'                     accrual_time = c(0, 12), n_accrual = c(40, 260))
#' print(ev)
#'
#' @export
print.rcpgsd_events <- function(x, digits = 4, ...) {
  S <- nrow(x$n_accrual)
  K <- length(x$events)

  cat("\nExpected Analysis Times and Regional Event Proportions\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("   Patients          : (%s) by region, total = %d\n",
              paste(rowSums(x$n_accrual), collapse = ", "),
              as.integer(sum(x$n_accrual))))
  cat(sprintf("   Accrual intervals : %s\n",
              paste(format(x$accrual_time, trim = TRUE), collapse = " | ")))
  cat(sprintf("   Control hazard    : (%s)\n",
              paste(signif(x$lambda_c, 4), collapse = ", ")))
  cat(sprintf("   True HR           : (%s)\n",
              paste(signif(x$hr, 4), collapse = ", ")))
  if (is.null(x$lambda_d) || all(x$lambda_d == 0)) {
    cat("   Drop-out hazard   : none\n")
  } else {
    cat(sprintf("   Drop-out hazard   : control (%s), experimental (%s)\n",
                paste(signif(x$lambda_d[, 1L], 4), collapse = ", "),
                paste(signif(x$lambda_d[, 2L], 4), collapse = ", ")))
  }
  cat(sprintf("   Allocation ratio  : r = %s\n", format(x$r)))
  cat(strrep("-", 70), "\n")

  an <- x$analysis
  df <- data.frame(Stage    = an$stage,
                   Events   = format(an$events),
                   Time     = formatC(an$time, format = "f", digits = 2),
                   Enrolled = formatC(an$enrolled, format = "f", digits = 1),
                   stringsAsFactors = FALSE)
  cat("Expected analysis times:\n\n")
  print(df, row.names = FALSE, right = TRUE)

  ev <- formatC(x$events_region, format = "f", digits = 1)
  fs <- formatC(x$fs, format = "f", digits = digits)
  tab <- matrix(paste0(ev, " (", fs, ")"), S, K)
  dimnames(tab) <- list(paste("Region", seq_len(S)),
                        paste("Stage", seq_len(K)))
  cat("\nExpected events by region (proportion):\n\n")
  print(noquote(tab), right = TRUE)
  cat("\nUse x$fs as 'fs' in rcpgsd() or rcpgsd_get_fs().\n\n")
  invisible(x)
}
