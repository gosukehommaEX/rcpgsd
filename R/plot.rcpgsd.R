#' Plot Method for rcpgsd Objects
#'
#' @description
#' Visualizes the results of \code{\link{rcpgsd}} by analysis, with one
#' panel per consistency method.
#' \itemize{
#'   \item \code{what = "rcp"}: the RCP conditional on stopping for efficacy
#'     at each analysis (solid line) and, optionally, the unconditional
#'     probability of consistency at each analysis (dashed line) and the
#'     overall RCP given success (dotted horizontal line).
#'   \item \code{what = "prob"}: the probability of stopping for efficacy
#'     and the joint probability of stopping for efficacy and meeting the
#'     consistency criterion at each analysis (bars).
#' }
#'
#' @param x An object of class \code{"rcpgsd"}.
#' @param ... Further arguments (not used).
#' @param what Character. \code{"rcp"} (default) or \code{"prob"}.
#' @param x_axis Character. Variable on the horizontal axis:
#'   \code{"stage"} (default), \code{"info_frac"} or \code{"events"}.
#'   Bars (\code{what = "prob"}) always use \code{"stage"}.
#' @param show_uncond Logical. Show the unconditional probability of
#'   consistency when \code{what = "rcp"}. Default is \code{TRUE}.
#' @param show_overall Logical. Show the overall RCP when
#'   \code{what = "rcp"}. Default is \code{TRUE}.
#' @param colours Character vector of two colour codes. Default uses the
#'   Okabe-Ito palette.
#' @param point_size Numeric scalar. Size of the points. Default \code{2.5}.
#' @param line_width Numeric scalar. Width of the lines. Default \code{0.8}.
#' @param base_size Numeric scalar. Base font size for
#'   \code{ggplot2::theme_bw()}. Default \code{11}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{rcpgsd}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_col
#' @importFrom ggplot2 facet_wrap scale_colour_manual scale_fill_manual
#' @importFrom ggplot2 scale_linetype_manual scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous labs theme_bw theme element_blank
#' @importFrom ggplot2 position_dodge
#' @importFrom rlang .data
#'
#' @examples
#' res <- rcpgsd(list(upper = c(3.7103, 2.5114, 1.9930)),
#'               events = c(150, 300, 450), fs = c(0.1, 0.45, 0.45),
#'               hr = 0.72)
#' plot(res)
#' plot(res, what = "prob")
#'
#' @export
plot.rcpgsd <- function(x,
                        ...,
                        what         = c("rcp", "prob"),
                        x_axis       = c("stage", "info_frac", "events"),
                        show_uncond  = TRUE,
                        show_overall = TRUE,
                        colours      = NULL,
                        point_size   = 2.5,
                        line_width   = 0.8,
                        base_size    = 11) {
  what   <- match.arg(what)
  x_axis <- match.arg(x_axis)
  if (is.null(colours)) colours <- c("#0072B2", "#D55E00")
  colours <- rep_len(colours, 2L)

  lab_method <- c(M1 = "Method 1", M2 = "Method 2")
  st <- x$stage
  st$method_lab <- factor(lab_method[st$method],
                          levels = unname(lab_method[x$method]))

  if (what == "prob") {
    df <- rbind(
      data.frame(stage = st$stage, method_lab = st$method_lab,
                 series = "Pr(efficacy)", value = st$p_eff),
      data.frame(stage = st$stage, method_lab = st$method_lab,
                 series = "Pr(efficacy & consistency)", value = st$p_joint)
    )
    df$series <- factor(df$series, levels = c("Pr(efficacy)",
                                              "Pr(efficacy & consistency)"))
    p <- ggplot2::ggplot(df, ggplot2::aes(x    = .data$stage,
                                          y    = .data$value,
                                          fill = .data$series)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7) +
      ggplot2::facet_wrap(~ method_lab) +
      ggplot2::scale_fill_manual(values = colours, name = NULL) +
      ggplot2::scale_x_continuous(breaks = unique(st$stage)) +
      ggplot2::labs(x = "Analysis", y = "Stagewise probability")
  } else {
    x_val <- switch(x_axis,
                    stage     = st$stage,
                    info_frac = st$info_frac,
                    events    = st$events)
    x_lab <- switch(x_axis,
                    stage     = "Analysis",
                    info_frac = "Information fraction",
                    events    = "Total number of events")
    df <- data.frame(xv = x_val, method_lab = st$method_lab,
                     series = "RCP given efficacy stop", value = st$rcp)
    if (show_uncond) {
      df <- rbind(df, data.frame(xv = x_val, method_lab = st$method_lab,
                                 series = "Unconditional Pr(consistency)",
                                 value = st$rcp_uncond))
    }
    df$series <- factor(df$series,
                        levels = c("RCP given efficacy stop",
                                   "Unconditional Pr(consistency)"))
    df <- df[!is.na(df$value), , drop = FALSE]

    p <- ggplot2::ggplot(df, ggplot2::aes(x        = .data$xv,
                                          y        = .data$value,
                                          colour   = .data$series,
                                          linetype = .data$series)) +
      ggplot2::geom_line(linewidth = line_width) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::facet_wrap(~ method_lab) +
      ggplot2::scale_colour_manual(values = colours, name = NULL,
                                   drop = FALSE) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed"),
                                     name = NULL, drop = FALSE) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = x_lab, y = "Probability")
    if (x_axis == "stage") {
      p <- p + ggplot2::scale_x_continuous(breaks = unique(st$stage))
    }
    if (show_overall) {
      ov <- x$overall
      ov$method_lab <- factor(lab_method[ov$method],
                              levels = unname(lab_method[x$method]))
      ov <- ov[!is.na(ov$rcp), , drop = FALSE]
      p <- p + ggplot2::geom_hline(data = ov,
                                   ggplot2::aes(yintercept = .data$rcp),
                                   inherit.aes = FALSE,
                                   linetype = "dotted", colour = "grey30")
    }
  }

  p + ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position  = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}
