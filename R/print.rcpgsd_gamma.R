#' Print Method for rcpgsd_gamma Objects
#'
#' @description
#' Prints the reference null RCPs and the adjusted consistency thresholds
#' computed by \code{\link{rcpgsd_get_gamma}}.
#'
#' @param x An object of class \code{"rcpgsd_gamma"}.
#' @param digits Integer. Number of decimal places. Default is \code{4}.
#' @param ... Further arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rcpgsd_get_gamma}}
#'
#' @examples
#' adj <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
#'                         events = c(200, 400), fs = c(0.15, 0.85))
#' print(adj)
#'
#' @export
print.rcpgsd_gamma <- function(x, digits = 4, ...) {
  fmt <- function(v) {
    ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
  }
  fmt_b <- function(v) {
    ifelse(is.finite(v), formatC(v, format = "f", digits = 3),
           ifelse(v > 0, "Inf", "-Inf"))
  }

  cat("\nAdjusted Consistency Thresholds under Group Sequential Design\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("   Consistency scale : %s\n",
              if (x$scale == "loghr") "log hazard ratio" else
                "hazard reduction (1 - HR)"))
  cat(sprintf("   Analyses          : K = %d (boundaries from %s)\n",
              x$design$K, x$design$source))
  cat(sprintf("   Futility bounds   : %s\n", x$design$futility))
  cat(sprintf("   Regional props.   : fs = (%s) at the final analysis\n",
              paste(format(x$fs[, x$design$K]), collapse = ", ")))
  cat(sprintf("   Adjustment type   : %s%s\n", x$type,
              if (x$allow_relax) " (relaxation allowed)" else
                " (thresholds only tightened)"))
  for (m in x$method) {
    base <- if (m == "M1") x$gamma_M1 else x$gamma_M2
    cat(sprintf("   %-17s : unadjusted threshold = %s, reference RCP0 = %s\n",
                if (m == "M1") "Method 1" else "Method 2",
                formatC(base, format = "f", digits = digits),
                fmt(x$ref[[m]])))
  }
  cat(strrep("-", 70), "\n")

  res <- x$result
  df <- data.frame(Stage    = res$stage,
                   Events   = format(round(res$events, 1)),
                   Upper    = fmt_b(res$upper),
                   Lower    = fmt_b(res$lower),
                   Pr_eff0  = formatC(res$p_eff0, format = "e", digits = 3),
                   stringsAsFactors = FALSE)
  nm <- c("Stage", "Events", "Upper", "Lower", "Pr(eff|H0)")
  for (m in x$method) {
    df[[paste0("u", m)]] <- fmt(res[[paste0("rcp0_", m, "_unadj")]])
    df[[paste0("g", m)]] <- fmt(res[[paste0("gamma_", m, "_adj")]])
    df[[paste0("a", m)]] <- fmt(res[[paste0("rcp0_", m, "_adj")]])
    nm <- c(nm, paste0("RCP0_", m, "(unadj)"), paste0("gamma_", m, "_adj"),
            paste0("RCP0_", m, "(adj)"))
  }
  colnames(df) <- nm
  cat("Null RCP given efficacy stop (RCP0) and adjusted thresholds:\n\n")
  print(df, row.names = FALSE, right = TRUE)

  cat("\nOverall null RCP given success:\n")
  for (i in seq_len(nrow(x$overall))) {
    cat(sprintf("   %s: unadjusted = %s, adjusted = %s\n",
                x$overall$method[i], fmt(x$overall$rcp0_unadj[i]),
                fmt(x$overall$rcp0_adj[i])))
  }
  cat("\n")
  invisible(x)
}
