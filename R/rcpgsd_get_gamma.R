#' Adjusted Consistency Thresholds under a Group Sequential Design
#'
#' @description
#' Computes adjusted consistency thresholds for Method 1 and Method 2 such
#' that the regional consistency probability (RCP) under the null
#' hypothesis (\eqn{\mathrm{HR}_s = 1} in all regions), conditional on
#' stopping for efficacy, equals its value in the corresponding fixed
#' design. This extends the adjustment of Homma (2026), developed for the
#' Bonferroni procedure with multiple primary endpoints, to group sequential
#' designs.
#'
#' @details
#' \strong{Motivation.} A trial that stops for efficacy at an early analysis
#' has crossed a boundary that is larger than the fixed-design critical
#' value (e.g., O'Brien-Fleming-type boundaries). Conditioning on this event
#' selects large overall effect estimates, and the regional consistency
#' criterion is then met more easily even when no treatment effect exists.
#' The null RCP at early analyses is therefore inflated relative to the
#' fixed design, in the same way as it is inflated by the Bonferroni-adjusted
#' critical value in Homma (2026).
#'
#' \strong{Reference value.} The reference null RCP for Method \eqn{m} is
#' \eqn{\phi^{(m)}_{\mathrm{Ref}} = \Pr(C \mid Z > z_{1-\alpha};
#' \mathrm{HR} = 1)} for a single-stage design with the final number of
#' events \code{events[K]}, the final regional proportions \code{fs[, K]},
#' the one-sided significance level \eqn{\alpha} of \code{design} and the
#' unadjusted thresholds \code{gamma_M1} and \code{gamma_M2}. On the log
#' hazard ratio scale the reference value does not depend on the number of
#' events. A user-specified reference value can be supplied through
#' \code{rcp_ref}.
#'
#' \strong{Adjustment.} With \code{type = "stagewise"}, a separate threshold
#' \eqn{\gamma_{m,k}} is obtained for each analysis \eqn{k} by solving
#' \eqn{\Pr(C_k \mid E_k; \mathrm{HR} = 1) = \phi^{(m)}_{\mathrm{Ref}}}
#' with \code{stats::uniroot()}. With \code{type = "overall"}, a single
#' threshold common to all analyses is obtained by solving
#' \eqn{\sum_k \Pr(E_k \cap C_k) / \sum_k \Pr(E_k) =
#' \phi^{(m)}_{\mathrm{Ref}}} under \eqn{\mathrm{HR} = 1}. Given that the
#' overall effect estimate is positive on the stopping event, the null RCP
#' is decreasing in the threshold, so the solution is unique.
#'
#' At some analyses (typically the final analysis) the null RCP with the
#' unadjusted threshold may be slightly below the reference value. With
#' \code{allow_relax = FALSE} (default) the unadjusted threshold is kept in
#' this case, so thresholds are only ever made stricter. With
#' \code{allow_relax = TRUE} the threshold is solved within the full search
#' interval and may be smaller than the unadjusted threshold.
#'
#' @inheritParams rcpgsd
#' @param design Group sequential design (see \code{\link{rcpgsd}}). The
#'   one-sided significance level must be available from the design object
#'   (element \code{alpha} for a list) unless \code{rcp_ref} is supplied.
#' @param gamma_M1 Numeric scalar. Unadjusted effect retention threshold for
#'   Method 1. Default is \code{0.5}.
#' @param gamma_M2 Numeric scalar. Unadjusted threshold for Method 2.
#'   Default is \code{0}.
#' @param type Character. \code{"stagewise"} (default) for analysis-specific
#'   thresholds or \code{"overall"} for a single threshold.
#' @param allow_relax Logical. Whether the adjusted threshold may be smaller
#'   than the unadjusted threshold. Default is \code{FALSE}.
#' @param rcp_ref Optional numeric vector of reference null RCPs, named
#'   \code{"M1"} and/or \code{"M2"} (or unnamed, in the order of
#'   \code{method}). \code{NULL} (default) computes the fixed-design
#'   reference.
#' @param interval_M1,interval_M2 Numeric vectors of length 2. Search
#'   intervals for the adjusted thresholds. Defaults are \code{c(0, 0.999)}.
#' @param tol Numeric scalar. Tolerance passed to \code{stats::uniroot()}.
#'   Default is \code{1e-6}.
#'
#' @return An object of class \code{"rcpgsd_gamma"}, a list containing:
#' \describe{
#'   \item{\code{result}}{Data frame with one row per analysis and columns
#'     \code{stage}, \code{events}, \code{info_frac}, \code{upper},
#'     \code{lower}, \code{p_eff0} (probability of stopping for efficacy
#'     under the null hypothesis), and, for each method requested,
#'     \code{rcp0_M1_unadj}, \code{gamma_M1_adj}, \code{rcp0_M1_adj}
#'     (and the same for \code{M2}).}
#'   \item{\code{overall}}{Data frame with one row per method and columns
#'     \code{method}, \code{rcp0_unadj} and \code{rcp0_adj} (overall null
#'     RCP given success, with unadjusted and adjusted thresholds).}
#'   \item{\code{ref}}{Named numeric vector of reference null RCPs.}
#'   \item{\code{design}, \code{events}, \code{fs}, \code{r},
#'     \code{method}, \code{scale}, \code{hr_method}, \code{gamma_M1},
#'     \code{gamma_M2}, \code{type}, \code{allow_relax}}{Settings used.}
#' }
#' The adjusted thresholds can be passed to \code{\link{rcpgsd}} through
#' \code{gamma_M1 = x$result$gamma_M1_adj} and
#' \code{gamma_M2 = x$result$gamma_M2_adj}.
#'
#' @references
#' Homma G (2026). A cautionary note on regional consistency evaluation with
#' multiple primary endpoints in multiregional clinical trials.
#' \emph{Pharmaceutical Statistics}, 25, e70111.
#'
#' @seealso \code{\link{rcpgsd}}, \code{\link{print.rcpgsd_gamma}},
#'   \code{\link{plot.rcpgsd_gamma}}
#'
#' @importFrom stats qnorm uniroot
#'
#' @examples
#' # Two analyses, Method 1 only
#' adj1 <- rcpgsd_get_gamma(list(upper = c(2.963, 1.969), alpha = 0.025),
#'                          events = c(200, 400), fs = c(0.15, 0.85),
#'                          method = "M1")
#' print(adj1)
#'
#' \donttest{
#' # Three analyses with O'Brien-Fleming-type boundaries, both methods
#' design <- list(upper = c(3.7103, 2.5114, 1.9930), alpha = 0.025)
#' adj <- rcpgsd_get_gamma(design, events = c(262, 523, 785),
#'                         fs = c(0.1, 0.45, 0.45))
#' print(adj)
#'
#' # RCP under the alternative with adjusted thresholds
#' res <- rcpgsd(design, events = c(262, 523, 785),
#'               fs = c(0.1, 0.45, 0.45), hr = exp(-0.2),
#'               gamma_M1 = adj$result$gamma_M1_adj,
#'               gamma_M2 = adj$result$gamma_M2_adj)
#' print(res)
#' }
#'
#' @export
rcpgsd_get_gamma <- function(design,
                             events,
                             fs,
                             r           = 1,
                             method      = c("M1", "M2"),
                             scale       = c("loghr", "hr"),
                             gamma_M1    = 0.5,
                             gamma_M2    = 0,
                             type        = c("stagewise", "overall"),
                             hr_method   = c("exact", "delta"),
                             futility    = c("obey", "ignore"),
                             allow_relax = FALSE,
                             rcp_ref     = NULL,
                             interval_M1 = c(0, 0.999),
                             interval_M2 = c(0, 0.999),
                             tol         = 1e-6,
                             n_nodes     = 64L) {

  # ========== Argument checks ==========
  method    <- match.arg(method, c("M1", "M2"), several.ok = TRUE)
  method    <- intersect(c("M1", "M2"), method)
  scale     <- match.arg(scale)
  type      <- match.arg(type)
  hr_method <- match.arg(hr_method)
  futility  <- match.arg(futility)
  if (!is.logical(allow_relax) || length(allow_relax) != 1L ||
      is.na(allow_relax)) {
    stop("'allow_relax' must be TRUE or FALSE.", call. = FALSE)
  }

  bnd <- .rcpgsd_bounds(design)
  K   <- bnd$K
  inp <- .rcpgsd_check_inputs(events, fs, NULL, r, K)
  S   <- nrow(inp$fs)

  if (!is.numeric(gamma_M1) || length(gamma_M1) != 1L) {
    stop("'gamma_M1' must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(gamma_M2) || length(gamma_M2) != 1L) {
    stop("'gamma_M2' must be a numeric scalar.", call. = FALSE)
  }
  .rcpgsd_check_gamma(gamma_M1, K, "gamma_M1", scale)
  .rcpgsd_check_gamma(gamma_M2, K, "gamma_M2", scale)
  base <- c(M1 = gamma_M1, M2 = gamma_M2)

  intervals <- list(M1 = interval_M1, M2 = interval_M2)
  for (m in method) {
    iv <- intervals[[m]]
    if (!is.numeric(iv) || length(iv) != 2L || iv[1L] >= iv[2L] ||
        iv[1L] < 0 || iv[2L] >= 1 ||
        base[[m]] < iv[1L] || base[[m]] > iv[2L]) {
      stop(sprintf(paste0("'interval_%s' must be an increasing numeric ",
                          "vector of length 2 within [0, 1) that contains ",
                          "the unadjusted threshold."), m), call. = FALSE)
    }
  }

  fut <- .rcpgsd_futility(bnd, futility)
  info_frac <- inp$events / inp$events[K]
  if (!is.null(bnd$timing) && max(abs(info_frac - bnd$timing)) > 0.01) {
    warning("The information fractions implied by 'events' differ from ",
            "the planned information fractions of 'design'; the boundaries ",
            "of 'design' are used without recalculation.", call. = FALSE)
  }

  # ========== Reference null RCP (fixed design) ==========
  if (is.null(rcp_ref)) {
    if (is.na(bnd$alpha)) {
      stop("The one-sided significance level could not be obtained from ",
           "'design'; supply 'design$alpha' or 'rcp_ref'.", call. = FALSE)
    }
    st_ref <- .rcpgsd_setup(inp$events[K], inp$fs[, K, drop = FALSE],
                            rep(1, S), inp$r)
    z_ref  <- stats::qnorm(1 - bnd$alpha)
    ref <- vapply(method, function(m) {
      .rcpgsd_prob_joint(st_ref, z_ref, -Inf, 1L, m, scale, base[[m]],
                         hr_method, n_nodes) /
        .rcpgsd_prob_stop(st_ref, z_ref, -Inf, 1L, "efficacy")
    }, numeric(1L))
  } else {
    if (!is.numeric(rcp_ref) || length(rcp_ref) != length(method) ||
        any(rcp_ref <= 0 | rcp_ref >= 1)) {
      stop("'rcp_ref' must contain one value in (0, 1) per method.",
           call. = FALSE)
    }
    if (is.null(names(rcp_ref))) {
      names(rcp_ref) <- method
    }
    if (!setequal(names(rcp_ref), method)) {
      stop("The names of 'rcp_ref' must match 'method'.", call. = FALSE)
    }
    ref <- rcp_ref[method]
  }
  names(ref) <- method

  # ========== Null model ==========
  st0    <- .rcpgsd_setup(inp$events, inp$fs, rep(1, S), inp$r)
  upper  <- bnd$upper
  lower  <- fut$lower
  p_eff0 <- vapply(seq_len(K), function(k) {
    .rcpgsd_prob_stop(st0, upper, lower, k, "efficacy")
  }, numeric(1L))
  p_eff_total <- sum(p_eff0)

  joint0 <- function(m, k, g) {
    if (p_eff0[k] <= 0) {
      return(0)
    }
    .rcpgsd_prob_joint(st0, upper, lower, k, m, scale, g, hr_method,
                       n_nodes)
  }

  # Solve f(g) = 0 for a decreasing function f on the search interval
  solve_gamma <- function(f, m, label) {
    g0 <- base[[m]]
    f0 <- f(g0)
    if (!allow_relax && f0 <= 0) {
      return(g0)
    }
    lo  <- if (allow_relax) intervals[[m]][1L] else g0
    hi  <- intervals[[m]][2L]
    flo <- if (lo == g0) f0 else f(lo)
    fhi <- f(hi)
    if (flo < 0 || fhi > 0) {
      warning(sprintf(paste0("No adjusted threshold for %s (%s) was found ",
                             "within the search interval; NA is returned."),
                      m, label), call. = FALSE)
      return(NA_real_)
    }
    if (flo == 0) return(lo)
    if (fhi == 0) return(hi)
    stats::uniroot(f, c(lo, hi), f.lower = flo, f.upper = fhi,
                   tol = tol)$root
  }

  result <- data.frame(stage     = seq_len(K),
                       events    = inp$events,
                       info_frac = info_frac,
                       upper     = upper,
                       lower     = lower,
                       p_eff0    = p_eff0)
  overall <- vector("list", length(method))
  names(overall) <- method

  for (m in method) {
    j_unadj <- vapply(seq_len(K), function(k) joint0(m, k, base[[m]]),
                      numeric(1L))
    rcp_unadj <- ifelse(p_eff0 > 0, j_unadj / p_eff0, NA_real_)

    if (type == "stagewise") {
      g_adj <- vapply(seq_len(K), function(k) {
        if (p_eff0[k] <= 0) {
          return(base[[m]])
        }
        f <- function(g) joint0(m, k, g) / p_eff0[k] - ref[[m]]
        solve_gamma(f, m, sprintf("stage %d", k))
      }, numeric(1L))
    } else {
      f <- function(g) {
        sum(vapply(seq_len(K), function(k) joint0(m, k, g),
                   numeric(1L))) / p_eff_total - ref[[m]]
      }
      g_adj <- rep(solve_gamma(f, m, "overall"), K)
    }

    j_adj <- vapply(seq_len(K), function(k) {
      if (is.na(g_adj[k])) NA_real_ else joint0(m, k, g_adj[k])
    }, numeric(1L))
    rcp_adj <- ifelse(p_eff0 > 0, j_adj / p_eff0, NA_real_)

    result[[paste0("rcp0_", m, "_unadj")]] <- rcp_unadj
    result[[paste0("gamma_", m, "_adj")]]  <- g_adj
    result[[paste0("rcp0_", m, "_adj")]]   <- rcp_adj

    overall[[m]] <- data.frame(method     = m,
                               rcp0_unadj = sum(j_unadj) / p_eff_total,
                               rcp0_adj   = sum(j_adj) / p_eff_total,
                               stringsAsFactors = FALSE)
  }

  overall <- do.call(rbind, overall)
  rownames(overall) <- NULL

  structure(
    list(result      = result,
         overall     = overall,
         ref         = ref,
         design      = list(K        = K,
                            upper    = upper,
                            lower    = lower,
                            binding  = bnd$binding,
                            futility = fut$label,
                            alpha    = bnd$alpha,
                            source   = bnd$source),
         events      = inp$events,
         fs          = inp$fs,
         r           = inp$r,
         method      = method,
         scale       = scale,
         hr_method   = hr_method,
         gamma_M1    = gamma_M1,
         gamma_M2    = gamma_M2,
         type        = type,
         allow_relax = allow_relax),
    class = "rcpgsd_gamma"
  )
}
