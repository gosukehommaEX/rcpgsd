#' Print Method for rcpgsd_sim Objects
#'
#' @description
#' Prints the settings, the analysis summaries and the stagewise and overall
#' regional consistency probabilities estimated by \code{\link{rcpgsd_sim}}.
#'
#' @param x An object of class \code{"rcpgsd_sim"}.
#' @param digits Integer. Number of decimal places for probabilities.
#'   Default is \code{4}.
#' @param ... Further arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{rcpgsd_sim}}
#'
#' @examples
#' sim <- rcpgsd_sim(list(upper = c(2.963, 1.969)), events = c(100, 200),
#'                   hr = 0.7, lambda_c = 0.05, accrual_time = c(0, 12),
#'                   n_accrual = c(40, 260), nsim = 100)
#' print(sim)
#'
#' @export
print.rcpgsd_sim <- function(x, digits = 4, ...) {
  fmt <- function(v) {
    ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
  }
  fmt_b <- function(v) {
    ifelse(is.finite(v), formatC(v, format = "f", digits = 3),
           ifelse(v > 0, "Inf", "-Inf"))
  }
  S <- nrow(x$n_accrual)

  cat("\nPatient-Level Simulation of Regional Consistency under",
      "Group Sequential Design\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("   Simulated trials  : nsim = %d\n", x$nsim))
  cat(sprintf("   Consistency scale : %s\n",
              if (x$scale == "loghr") "log hazard ratio" else
                "hazard reduction (1 - HR)"))
  est_lab <- c(cox = "Cox regression", peto = "Peto's method",
               pike = "Pike's method")
  cat(sprintf("   HR estimator      : %s (overall and by region)\n",
              est_lab[[x$estimator]]))
  cat("   Test statistic    : log-rank (unstratified)\n")
  cat(sprintf("   Analyses          : K = %d (boundaries from %s)\n",
              x$design$K, x$design$source))
  cat(sprintf("   Futility bounds   : %s\n", x$design$futility))
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
  df_an <- data.frame(Stage     = an$stage,
                      Events    = format(an$events),
                      Time_mean = formatC(an$time_mean, format = "f",
                                          digits = 2),
                      Time_sd   = formatC(an$time_sd, format = "f",
                                          digits = 2),
                      Events_R1 = formatC(an$events_r1_mean, format = "f",
                                          digits = 1),
                      Prop_R1   = formatC(x$fs_mean[1L, ], format = "f",
                                          digits = 4),
                      No_est_R1 = fmt(an$p_no_est_r1),
                      Inf_R1    = fmt(an$p_inf_r1),
                      stringsAsFactors = FALSE)
  colnames(df_an) <- c("Stage", "Events", "Time(mean)", "Time(SD)",
                       "Events(R1)", "fs(R1)", "Pr(no est R1)",
                       "Pr(inf est R1)")
  cat("Analyses (means over simulated trials):\n\n")
  print(df_an, row.names = FALSE, right = TRUE)
  if (!is.null(an$p_not_reached) && any(an$p_not_reached > 0)) {
    cat(sprintf(paste0("\n   Not reached because of drop-out: ",
                       "Pr = (%s) by stage\n"),
                paste(fmt(an$p_not_reached), collapse = ", ")))
  }
  cat("\n")

  lab <- c(M1 = "Method 1", M2 = "Method 2")
  for (m in x$method) {
    sub <- x$stage[x$stage$method == m, , drop = FALSE]
    ov  <- x$overall[x$overall$method == m, , drop = FALSE]
    gname <- if (m == "M1") "gamma_M1" else "gamma_M2"
    cat(sprintf("%s (%s):\n\n", lab[[m]], gname))
    df <- data.frame(Stage       = sub$stage,
                     Upper       = fmt_b(sub$upper),
                     Lower       = fmt_b(sub$lower),
                     Gamma       = formatC(sub$gamma, format = "f",
                                           digits = 4),
                     Pr_eff      = fmt(sub$p_eff),
                     Pr_eff_cons = fmt(sub$p_joint),
                     RCP_cond    = fmt(sub$rcp),
                     SE          = fmt(sub$se_rcp),
                     RCP_uncond  = fmt(sub$rcp_uncond),
                     stringsAsFactors = FALSE)
    colnames(df) <- c("Stage", "Upper", "Lower", gname, "Pr(eff)",
                      "Pr(eff&cons)", "RCP|eff", "SE", "RCP(uncond)")
    print(df, row.names = FALSE, right = TRUE)
    cat(sprintf(paste0("\n   Overall: Pr(efficacy) = %s, ",
                       "Pr(efficacy & consistency) = %s, RCP = %s\n\n"),
                fmt(ov$power), fmt(ov$p_joint), fmt(ov$rcp)))
  }
  cat("RCP|eff    : RCP conditional on stopping for efficacy at the stage\n")
  cat("SE         : Monte Carlo standard error of RCP|eff\n")
  cat("RCP(uncond): Pr(consistency) at the stage, ignoring the boundaries\n")
  cat(sprintf("Regions    : S = %d (region 1 = region of interest)\n", S))
  cat("\n")
  invisible(x)
}
