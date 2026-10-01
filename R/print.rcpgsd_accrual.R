#' Print Method for rcpgsd_accrual Objects
#'
#' @description
#' Prints the regional accrual constructed by \code{\link{rcpgsd_accrual}}:
#' the start of enrolment and the number of patients in each region, and the
#' numbers of patients and accrual rates by region and interval.
#'
#' @param x An object of class \code{"rcpgsd_accrual"}.
#' @param digits Integer. Number of significant digits for accrual rates.
#'   Default is \code{4}.
#' @param ... Further arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rcpgsd_accrual}}
#'
#' @examples
#' acc <- rcpgsd_accrual(n_region = c(60, 270, 270), start = c(6, 0, 0),
#'                       accrual_time = c(0, 12), rate_overall = 50)
#' print(acc)
#'
#' @export
print.rcpgsd_accrual <- function(x, digits = 4, ...) {
  S  <- length(x$n_region)
  tm <- format(x$accrual_time, trim = TRUE)
  M  <- length(tm) - 1L
  iv <- paste0("[", tm[-(M + 1L)], ", ", tm[-1L], ")")

  cat("\nRegional Accrual under the Overall Accrual Plan\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("   Overall plan      : %s patients, accrual intervals %s\n",
              format(sum(x$n_region)),
              paste(format(x$plan$accrual_time, trim = TRUE),
                    collapse = " | ")))
  cat(sprintf("   Overall rate      : (%s) per unit of time\n",
              paste(signif(x$plan$rate_overall, digits), collapse = ", ")))
  cat(sprintf("   Patients by region: (%s)\n",
              paste(format(x$n_region, trim = TRUE), collapse = ", ")))
  cat(sprintf("   Start by region   : (%s)\n",
              paste(format(x$start, trim = TRUE), collapse = ", ")))
  cat(strrep("-", 70), "\n")

  tab <- matrix(format(x$n_accrual, trim = TRUE), S, M,
                dimnames = list(paste("Region", seq_len(S)), iv))
  tab <- rbind(tab, Overall = format(colSums(x$n_accrual), trim = TRUE))
  cat("Number of patients by interval (n_accrual):\n\n")
  print(noquote(tab), right = TRUE)

  rt <- matrix(format(signif(x$rate, digits), trim = TRUE), S, M,
               dimnames = list(paste("Region", seq_len(S)), iv))
  rt <- rbind(rt, Overall = format(signif(x$rate_overall, digits),
                                   trim = TRUE))
  cat("\nAccrual rate by interval (patients per unit of time):\n\n")
  print(noquote(rt), right = TRUE)
  cat("\nUse x$accrual_time and x$n_accrual in rcpgsd_events() or",
      "rcpgsd_sim().\n\n")
  invisible(x)
}
