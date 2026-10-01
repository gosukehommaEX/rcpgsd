#' Rescale Regional Event Proportions for a Given Share of Region 1
#'
#' @description
#' Internal helper that rescales a template matrix of regional event
#' proportions so that region 1 has proportion \code{x} at the final
#' analysis. The proportion of region 1 at analysis \eqn{k} is
#' \eqn{x f_{1,k} / f_{1,K}}, which keeps the shape of the template across
#' analyses (e.g., delayed enrolment). The proportions of the other regions
#' are rescaled so that their relative sizes are those of the template and
#' each column sums to 1.
#'
#' @param x Numeric scalar. Proportion of events in region 1 at the final
#'   analysis.
#' @param tmpl Numeric \code{S}-by-\code{K} template matrix of regional event
#'   proportions.
#'
#' @return Numeric \code{S}-by-\code{K} matrix of regional event proportions.
#'
#' @noRd
.rcpgsd_fs_scale <- function(x, tmpl) {
  K  <- ncol(tmpl)
  f1 <- x * tmpl[1L, ] / tmpl[1L, K]
  others <- sweep(tmpl[-1L, , drop = FALSE], 2L, (1 - f1) / (1 - tmpl[1L, ]),
                  "*")
  rbind(f1, others, deparse.level = 0)
}
