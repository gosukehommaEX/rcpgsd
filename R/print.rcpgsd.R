#' Print Method for rcpgsd Objects
#'
#' @description
#' Prints the settings, the stagewise probabilities and the overall
#' regional consistency probabilities computed by \code{\link{rcpgsd}}.
#'
#' @param x An object of class \code{"rcpgsd"}.
#' @param digits Integer. Number of decimal places for probabilities.
#'   Default is \code{4}.
#' @param ... Further arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rcpgsd}}
#'
#' @examples
#' res <- rcpgsd(list(upper = c(2.963, 1.969)), events = c(200, 400),
#'               fs = c(0.15, 0.85), hr = 0.7)
#' print(res)
#'
#' @export
print.rcpgsd <- function(x, digits = 4, ...) {
  fmt <- function(v) {
    ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
  }
  fmt_b <- function(v) {
    ifelse(is.finite(v), formatC(v, format = "f", digits = 3),
           ifelse(v > 0, "Inf", "-Inf"))
  }
  S <- nrow(x$fs)

  cat("\nRegional Consistency Probability under Group Sequential Design\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("   Approach          : %s\n",
              if (x$approach == "formula") {
                "closed form (multivariate normal)"
              } else {
                sprintf("Monte Carlo simulation (nsim = %d)", x$nsim)
              }))
  scale_lab <- if (x$scale == "loghr") {
    "log hazard ratio"
  } else if ("M1" %in% x$method && x$approach == "formula") {
    sprintf("hazard reduction (1 - HR), Method 1 by %s",
            if (x$hr_method == "exact") "numerical integration" else
              "delta method")
  } else {
    "hazard reduction (1 - HR)"
  }
  cat(sprintf("   Consistency scale : %s\n", scale_lab))
  cat(sprintf("   Analyses          : K = %d (boundaries from %s)\n",
              x$design$K, x$design$source))
  cat(sprintf("   Futility bounds   : %s\n", x$design$futility))
  cat(sprintf("   Regions           : S = %d (region 1 = region of interest)\n",
              S))
  cat(sprintf("   True HR           : (%s)\n",
              paste(formatC(x$hr, format = "g", digits = 4), collapse = ", ")))
  cat(sprintf("   Allocation ratio  : r = %s\n", format(x$r)))
  cat(sprintf("   Total events      : (%s)\n",
              paste(format(x$events), collapse = ", ")))
  cat(strrep("-", 70), "\n")

  lab <- c(M1 = "Method 1", M2 = "Method 2")
  for (m in x$method) {
    sub <- x$stage[x$stage$method == m, , drop = FALSE]
    ov  <- x$overall[x$overall$method == m, , drop = FALSE]
    gname <- if (m == "M1") "gamma_M1" else "gamma_M2"
    cat(sprintf("%s (%s):\n\n", lab[[m]], gname))
    df <- data.frame(Stage        = sub$stage,
                     Events       = format(round(sub$events, 1)),
                     InfoFrac     = formatC(sub$info_frac, format = "f",
                                            digits = 3),
                     Events_R1    = formatC(sub$events_r1, format = "f",
                                            digits = 1),
                     Upper        = fmt_b(sub$upper),
                     Lower        = fmt_b(sub$lower),
                     Gamma        = formatC(sub$gamma, format = "f",
                                            digits = 4),
                     Pr_eff       = fmt(sub$p_eff),
                     Pr_eff_cons  = fmt(sub$p_joint),
                     RCP_cond     = fmt(sub$rcp),
                     RCP_uncond   = fmt(sub$rcp_uncond),
                     stringsAsFactors = FALSE)
    colnames(df) <- c("Stage", "Events", "InfoFrac", "Events(R1)",
                      "Upper", "Lower", gname, "Pr(eff)",
                      "Pr(eff&cons)", "RCP|eff", "RCP(uncond)")
    print(df, row.names = FALSE, right = TRUE)
    cat(sprintf(paste0("\n   Overall: Pr(efficacy) = %s, ",
                       "Pr(efficacy & consistency) = %s, RCP = %s\n\n"),
                fmt(ov$power), fmt(ov$p_joint), fmt(ov$rcp)))
  }
  cat("RCP|eff    : RCP conditional on stopping for efficacy at the stage\n")
  cat("RCP(uncond): Pr(consistency) at the stage, ignoring the boundaries\n")
  cat("\n")
  invisible(x)
}
