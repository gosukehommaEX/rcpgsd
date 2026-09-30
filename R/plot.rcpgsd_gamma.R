#' Plot Method for rcpgsd_gamma Objects
#'
#' @description
#' Visualizes the results of \code{\link{rcpgsd_get_gamma}} by analysis,
#' with one column per consistency method. The upper row shows the adjusted
#' thresholds together with the unadjusted threshold (dashed line). The
#' lower row shows the null RCP given efficacy stopping with unadjusted and
#' adjusted thresholds together with the reference value (dashed line).
#'
#' @param x An object of class \code{"rcpgsd_gamma"}.
#' @param ... Further arguments (not used).
#' @param colours Character vector of two colour codes. Default uses the
#'   Okabe-Ito palette.
#' @param point_size Numeric scalar. Size of the points. Default \code{2.5}.
#' @param line_width Numeric scalar. Width of the lines. Default \code{0.8}.
#' @param base_size Numeric scalar. Base font size for
#'   \code{ggplot2::theme_bw()}. Default \code{11}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{rcpgsd_get_gamma}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline facet_grid
#' @importFrom ggplot2 scale_colour_manual scale_x_continuous labs theme_bw
#' @importFrom ggplot2 theme element_blank
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' adj <- rcpgsd_get_gamma(list(upper = c(3.7103, 2.5114, 1.9930),
#'                              alpha = 0.025),
#'                         events = c(262, 523, 785),
#'                         fs = c(0.1, 0.45, 0.45))
#' plot(adj)
#' }
#'
#' @export
plot.rcpgsd_gamma <- function(x,
                              ...,
                              colours    = NULL,
                              point_size = 2.5,
                              line_width = 0.8,
                              base_size  = 11) {
  if (is.null(colours)) colours <- c("#0072B2", "#D55E00")
  colours <- rep_len(colours, 2L)

  lab_method <- c(M1 = "Method 1", M2 = "Method 2")
  lab_panel  <- c("Threshold", "Null RCP given efficacy stop")
  res <- x$result

  pts <- list()
  ref <- list()
  for (m in x$method) {
    base <- if (m == "M1") x$gamma_M1 else x$gamma_M2
    pts[[length(pts) + 1L]] <- data.frame(
      stage = res$stage, method = lab_method[[m]], panel = lab_panel[1L],
      series = "Adjusted", value = res[[paste0("gamma_", m, "_adj")]])
    pts[[length(pts) + 1L]] <- data.frame(
      stage = res$stage, method = lab_method[[m]], panel = lab_panel[2L],
      series = "Unadjusted", value = res[[paste0("rcp0_", m, "_unadj")]])
    pts[[length(pts) + 1L]] <- data.frame(
      stage = res$stage, method = lab_method[[m]], panel = lab_panel[2L],
      series = "Adjusted", value = res[[paste0("rcp0_", m, "_adj")]])
    ref[[length(ref) + 1L]] <- data.frame(
      method = lab_method[[m]], panel = lab_panel[1L], value = base)
    ref[[length(ref) + 1L]] <- data.frame(
      method = lab_method[[m]], panel = lab_panel[2L], value = x$ref[[m]])
  }
  df  <- do.call(rbind, pts)
  ref <- do.call(rbind, ref)
  df  <- df[!is.na(df$value), , drop = FALSE]
  lev_m <- unname(lab_method[x$method])
  df$method   <- factor(df$method, levels = lev_m)
  ref$method  <- factor(ref$method, levels = lev_m)
  df$panel    <- factor(df$panel, levels = lab_panel)
  ref$panel   <- factor(ref$panel, levels = lab_panel)
  df$series   <- factor(df$series, levels = c("Unadjusted", "Adjusted"))

  ggplot2::ggplot(df, ggplot2::aes(x      = .data$stage,
                                   y      = .data$value,
                                   colour = .data$series)) +
    ggplot2::geom_hline(data = ref,
                        ggplot2::aes(yintercept = .data$value),
                        inherit.aes = FALSE,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_grid(panel ~ method, scales = "free_y") +
    ggplot2::scale_colour_manual(values = colours, name = NULL,
                                 drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = res$stage) +
    ggplot2::labs(x = "Analysis", y = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position  = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}
