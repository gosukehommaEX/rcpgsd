#' Print Method for rcpgsd_fs Objects
#'
#' @description
#' Prints the required proportion and number of events in the region of
#' interest computed by \code{\link{rcpgsd_get_fs}}, together with the
#' stagewise probabilities at that proportion.
#'
#' @param x An object of class \code{"rcpgsd_fs"}.
#' @param digits Integer. Number of decimal places for probabilities and
#'   proportions. Default is \code{4}.
#' @param ... Further arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rcpgsd_get_fs}}
#'
#' @examples
#' req <- rcpgsd_get_fs(list(upper = c(2.963, 1.969)), events = c(200, 400),
#'                      fs = c(0.15, 0.85), hr = 0.7, method = "M1")
#' print(req)
#'
#' @export
print.rcpgsd_fs <- function(x, digits = 4, ...) {
  fmt <- function(v) {
    ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
  }
  type_lab <- switch(x$type,
                     overall       = "overall RCP given success",
                     stagewise     = "smallest RCP given efficacy stopping",
                     unconditional = "unconditional Pr(consistency) at the final analysis")
  K <- x$design$K

  cat("\nRequired Proportion of Events in the Region of Interest\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("   Criterion         : %s >= %s\n", type_lab,
              format(x$target)))
  cat(sprintf("   Consistency scale : %s\n",
              if (x$scale == "loghr") "log hazard ratio" else
                "hazard reduction (1 - HR)"))
  cat(sprintf("   Analyses          : K = %d (boundaries from %s)\n", K,
              x$design$source))
  cat(sprintf("   Futility bounds   : %s\n", x$design$futility))
  cat(sprintf("   Total events      : (%s)\n",
              paste(format(x$events, trim = TRUE), collapse = ", ")))
  cat(sprintf("   True HR           : (%s)\n",
              paste(signif(x$hr, 4), collapse = ", ")))
  cat(sprintf("   Allocation ratio  : r = %s\n", format(x$r)))
  cat(strrep("-", 70), "\n")

  res <- x$result
  df <- data.frame(Method   = ifelse(res$method == "M1", "Method 1",
                                     "Method 2"),
                   fs_R1    = fmt(res$fs_r1),
                   Events   = ifelse(is.na(res$events_r1), "NA",
                                     formatC(res$events_r1, format = "f",
                                             digits = 2)),
                   Ceiling  = ifelse(is.na(res$events_r1_ceiling), "NA",
                                     format(res$events_r1_ceiling)),
                   Value    = fmt(res$value_ceiling),
                   stringsAsFactors = FALSE)
  colnames(df) <- c("Method", "fs(R1) final", "Events(R1) final",
                    "Rounded up", "Criterion at rounded")
  cat("Required proportion and number of events in region 1:\n\n")
  print(df, row.names = FALSE, right = TRUE)
  if (any(res$at_lower)) {
    cat("\n   Note: the target is already met at the lower limit of",
        "'interval' for", paste(res$method[res$at_lower], collapse = ", "),
        "\n")
  }

  for (m in x$method) {
    f <- x$fit[[m]]
    if (is.null(f)) next
    sub <- f$stage
    cat(sprintf("\n%s at the required proportion:\n\n",
                if (m == "M1") "Method 1" else "Method 2"))
    dfm <- data.frame(Stage      = sub$stage,
                      Events     = format(sub$events, trim = TRUE),
                      fs_R1      = fmt(f$fs[1L, sub$stage]),
                      Events_R1  = formatC(sub$events_r1, format = "f",
                                           digits = 1),
                      Pr_eff     = fmt(sub$p_eff),
                      RCP_cond   = fmt(sub$rcp),
                      RCP_uncond = fmt(sub$rcp_uncond),
                      stringsAsFactors = FALSE)
    colnames(dfm) <- c("Stage", "Events", "fs(R1)", "Events(R1)",
                       "Pr(eff)", "RCP|eff", "RCP(uncond)")
    print(dfm, row.names = FALSE, right = TRUE)
    cat(sprintf("\n   Overall RCP given success = %s\n",
                fmt(f$overall$rcp)))
  }
  cat("\n")
  invisible(x)
}
