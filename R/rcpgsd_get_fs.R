#' Required Proportion of Events in the Region of Interest
#'
#' @description
#' Computes the smallest proportion of events in region 1 (the region of
#' interest, e.g., Japan) at the final analysis such that a regional
#' consistency probability (RCP) reaches a target value under a group
#' sequential design. This is the inverse problem of \code{\link{rcpgsd}}
#' and corresponds to the regional sample size determination of Ikeda and
#' Bretz (2010) and Hayashi and Itoh (2018) for fixed designs.
#'
#' @details
#' \strong{Regional proportions.} The argument \code{fs} is a template. For
#' a candidate proportion \eqn{x} of region 1 at the final analysis, the
#' proportion of region 1 at analysis \eqn{k} is \eqn{x f_{1,k} / f_{1,K}},
#' so that the shape of the template across analyses (e.g., delayed
#' enrolment obtained with \code{\link{rcpgsd_events}}) is kept. The other
#' regions share the remaining events in the proportions of the template.
#'
#' \strong{Criterion} (\code{type}):
#' \itemize{
#'   \item \code{"overall"}: the overall RCP given success,
#'     \eqn{\sum_k \Pr(E_k \cap C_k) / \sum_k \Pr(E_k)}.
#'   \item \code{"stagewise"}: the smallest RCP given efficacy stopping,
#'     \eqn{\min_k \Pr(C_k \mid E_k)}, over the analyses at which stopping
#'     for efficacy is possible.
#'   \item \code{"unconditional"}: the unconditional probability of
#'     consistency at the final analysis, \eqn{\Pr(C_K)}. For a single-stage
#'     design this is the criterion used by Hayashi and Itoh (2018) and, on
#'     the log hazard ratio scale, by Ikeda and Bretz (2010).
#' }
#'
#' \strong{Search.} The criterion is evaluated at \code{interval[1]} and
#' then at doubled values of \eqn{x} until the target is reached, and the
#' solution is refined with \code{stats::uniroot()} within the last
#' bracket. This finds the smallest solution when the criterion increases
#' with \eqn{x} below the solution, which holds for Method 1 and, for
#' Method 2, as long as region 1 is not the largest region. If the target
#' is already met at \code{interval[1]}, that value is returned. If it is
#' not reached within the search interval, \code{NA} is returned with a
#' warning.
#'
#' The consistency thresholds are kept fixed during the search.
#'
#' @inheritParams rcpgsd
#' @param fs Numeric vector of length \code{S}, or numeric
#'   \code{S}-by-\code{K} matrix. Template of regional event proportions
#'   (see Details). The proportion of region 1 is replaced during the
#'   search, and only its shape across analyses is used.
#' @param target Numeric scalar in \eqn{(0, 1)}. Target value of the
#'   criterion. Default is \code{0.8}.
#' @param type Character. Criterion to be met: \code{"overall"} (default),
#'   \code{"stagewise"} or \code{"unconditional"} (see Details).
#' @param interval Numeric vector of length 2. Search interval for the
#'   proportion of events in region 1 at the final analysis. Default is
#'   \code{c(0.01, 0.99)}; the upper limit is reduced if necessary so that
#'   the proportion of region 1 stays below 1 at every analysis and the
#'   numbers of events in the other regions do not decrease across
#'   analyses.
#' @param tol Numeric scalar. Tolerance passed to \code{stats::uniroot()}.
#'   Default is \code{1e-6}.
#'
#' @return An object of class \code{"rcpgsd_fs"}, a list containing:
#' \describe{
#'   \item{\code{result}}{Data frame with one row per method and columns
#'     \code{method}, \code{fs_r1} (required proportion of events in region 1
#'     at the final analysis), \code{events_r1} (corresponding number of
#'     events), \code{events_r1_ceiling} (that number rounded up to an
#'     integer), \code{value_ceiling} (value of the criterion at the rounded
#'     number) and \code{at_lower} (\code{TRUE} if the target is already met
#'     at \code{interval[1]}).}
#'   \item{\code{fit}}{Named list with, for each method, the
#'     \code{\link{rcpgsd}} result at the required proportion (or
#'     \code{NULL} if no solution was found).}
#'   \item{\code{design}, \code{events}, \code{fs}, \code{hr}, \code{r},
#'     \code{method}, \code{type}, \code{target}, \code{scale},
#'     \code{hr_method}, \code{gamma_M1}, \code{gamma_M2}}{Settings used,
#'     with \code{fs} the template.}
#' }
#'
#' @references
#' Ikeda K, Bretz F (2010). Sample size and proportion of Japanese patients
#' in multi-regional trials. \emph{Pharmaceutical Statistics}, 9, 207-216.
#'
#' Hayashi N, Itoh Y (2018). A re-examination of Japanese sample size
#' calculation for multi-regional clinical trial evaluating survival
#' endpoint. \emph{Japanese Journal of Biometrics}, 38, 79-92.
#'
#' @seealso \code{\link{rcpgsd}}, \code{\link{rcpgsd_events}},
#'   \code{\link{print.rcpgsd_fs}}
#'
#' @importFrom stats uniroot
#'
#' @examples
#' # Two analyses; proportion of region 1 constant across analyses
#' design <- list(upper = c(2.963, 1.969), alpha = 0.025)
#' req <- rcpgsd_get_fs(design, events = c(200, 400),
#'                      fs = c(0.15, 0.85), hr = 0.7)
#' print(req)
#'
#' # Fixed design, hazard reduction scale, unconditional probability
#' # (Hayashi and Itoh, 2018, Table 1: about 29.8 Japanese events)
#' req_hi <- rcpgsd_get_fs(list(upper = qnorm(0.975)), events = 161,
#'                         fs = c(0.2, 0.8), hr = 0.6, method = "M1",
#'                         type = "unconditional", scale = "hr")
#' req_hi$result
#'
#' @export
rcpgsd_get_fs <- function(design,
                          events,
                          fs,
                          hr,
                          target    = 0.8,
                          r         = 1,
                          method    = c("M1", "M2"),
                          type      = c("overall", "stagewise",
                                        "unconditional"),
                          scale     = c("loghr", "hr"),
                          gamma_M1  = 0.5,
                          gamma_M2  = 0,
                          hr_method = c("exact", "delta"),
                          futility  = c("obey", "ignore"),
                          interval  = c(0.01, 0.99),
                          tol       = 1e-6,
                          n_nodes   = 64L) {

  # ========== Argument checks ==========
  method    <- match.arg(method, c("M1", "M2"), several.ok = TRUE)
  method    <- intersect(c("M1", "M2"), method)
  type      <- match.arg(type)
  scale     <- match.arg(scale)
  hr_method <- match.arg(hr_method)
  futility  <- match.arg(futility)
  if (missing(hr) || is.null(hr)) {
    stop("'hr' must be supplied.", call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L || target <= 0 ||
      target >= 1) {
    stop("'target' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(interval) || length(interval) != 2L ||
      interval[1L] <= 0 || interval[2L] >= 1 ||
      interval[1L] >= interval[2L]) {
    stop("'interval' must be an increasing numeric vector of length 2 ",
         "within (0, 1).", call. = FALSE)
  }

  bnd  <- .rcpgsd_bounds(design)
  K    <- bnd$K
  inp  <- .rcpgsd_check_inputs(events, fs, hr, r, K)
  tmpl <- inp$fs
  gamma_M1 <- .rcpgsd_check_gamma(gamma_M1, K, "gamma_M1", scale)
  gamma_M2 <- .rcpgsd_check_gamma(gamma_M2, K, "gamma_M2", scale)

  fut <- .rcpgsd_futility(bnd, futility)
  info_frac <- inp$events / inp$events[K]
  if (!is.null(bnd$timing) && max(abs(info_frac - bnd$timing)) > 0.01) {
    warning("The information fractions implied by 'events' differ from ",
            "the planned information fractions of 'design'; the boundaries ",
            "of 'design' are used without recalculation.", call. = FALSE)
  }

  # Upper limit of the search: region 1 must stay below 1 at every
  # analysis, and the numbers of events in the other regions, which shrink
  # as region 1 grows, must remain non-decreasing across analyses. With
  # a_k = f_{s,k} d_k / (1 - f_{1,k}) and p_k = f_{1,k} / f_{1,K}, the
  # condition for analyses k and k + 1 is
  # a_{k+1} - a_k >= x (a_{k+1} p_{k+1} - a_k p_k).
  profile <- tmpl[1L, ] / tmpl[1L, K]
  x_feas  <- Inf
  if (K > 1L) {
    for (s in seq_len(nrow(tmpl))[-1L]) {
      a <- tmpl[s, ] * inp$events / (1 - tmpl[1L, ])
      for (k in seq_len(K - 1L)) {
        cc <- a[k + 1L] * profile[k + 1L] - a[k] * profile[k]
        if (cc > 0) {
          x_feas <- min(x_feas, (a[k + 1L] - a[k]) / cc)
        }
      }
    }
  }
  x_max <- min(interval[2L], (1 - 1e-8) / max(profile),
               x_feas * (1 - 1e-8))
  if (interval[1L] >= x_max) {
    stop("'interval' is not compatible with the shape of 'fs' for region 1.",
         call. = FALSE)
  }

  # ========== Criterion as a function of the proportion of region 1 ==========
  criterion <- function(x, m) {
    fs_x <- .rcpgsd_fs_scale(x, tmpl)
    ev_x <- sweep(fs_x, 2L, inp$events, "*")
    if (K > 1L && any(ev_x[, -1L, drop = FALSE] -
                      ev_x[, -K, drop = FALSE] < -1e-8)) {
      stop(sprintf(paste0("For a proportion of %.4f in region 1, the ",
                          "implied regional numbers of events decrease ",
                          "across analyses; modify 'fs' or 'interval'."), x),
           call. = FALSE)
    }
    st  <- .rcpgsd_setup(inp$events, fs_x, inp$hr, inp$r)
    gam <- if (m == "M1") gamma_M1 else gamma_M2
    if (type == "unconditional") {
      return(.rcpgsd_prob_uncond(st, K, m, scale, gam[K], hr_method))
    }
    p_eff <- vapply(seq_len(K), function(k) {
      .rcpgsd_prob_stop(st, bnd$upper, fut$lower, k, "efficacy")
    }, numeric(1L))
    p_joint <- vapply(seq_len(K), function(k) {
      if (p_eff[k] <= 0) {
        return(0)
      }
      .rcpgsd_prob_joint(st, bnd$upper, fut$lower, k, m, scale, gam[k],
                         hr_method, n_nodes)
    }, numeric(1L))
    if (type == "overall") {
      sum(p_joint) / sum(p_eff)
    } else {
      ok <- p_eff > 0
      min(p_joint[ok] / p_eff[ok])
    }
  }

  # ========== Search ==========
  solve_x <- function(m) {
    f  <- function(x) criterion(x, m) - target
    lo <- interval[1L]
    flo <- f(lo)
    if (flo >= 0) {
      return(list(x = lo, at_lower = TRUE))
    }
    repeat {
      hi  <- min(2 * lo, x_max)
      fhi <- f(hi)
      if (fhi >= 0) break
      if (hi >= x_max) {
        warning(sprintf(paste0("The target was not reached for %s within ",
                               "the search interval (upper limit %.4f); NA ",
                               "is returned."), m, x_max),
                call. = FALSE)
        return(list(x = NA_real_, at_lower = FALSE))
      }
      lo  <- hi
      flo <- fhi
    }
    x <- if (fhi == 0) hi else
      stats::uniroot(f, c(lo, hi), f.lower = flo, f.upper = fhi,
                     tol = tol)$root
    list(x = x, at_lower = FALSE)
  }

  result <- vector("list", length(method))
  fit    <- vector("list", length(method))
  names(result) <- method
  names(fit)    <- method
  for (m in method) {
    sol <- solve_x(m)
    x   <- sol$x
    if (is.na(x)) {
      result[[m]] <- data.frame(method = m, fs_r1 = NA_real_,
                                events_r1 = NA_real_,
                                events_r1_ceiling = NA_real_,
                                value_ceiling = NA_real_,
                                at_lower = FALSE, stringsAsFactors = FALSE)
      fit[m] <- list(NULL)
      next
    }
    d_K   <- inp$events[K]
    e_r1  <- x * d_K
    e_ceil <- ceiling(e_r1 - 1e-8)
    x_ceil <- min(e_ceil / d_K, x_max)
    result[[m]] <- data.frame(method            = m,
                              fs_r1             = x,
                              events_r1         = e_r1,
                              events_r1_ceiling = e_ceil,
                              value_ceiling     = criterion(x_ceil, m),
                              at_lower          = sol$at_lower,
                              stringsAsFactors  = FALSE)
    fit[[m]] <- suppressWarnings(
      rcpgsd(design, inp$events, .rcpgsd_fs_scale(x, tmpl), inp$hr,
             r = inp$r, method = m, scale = scale, gamma_M1 = gamma_M1,
             gamma_M2 = gamma_M2, hr_method = hr_method,
             futility = futility, n_nodes = n_nodes)
    )
  }
  result <- do.call(rbind, result)
  rownames(result) <- NULL

  structure(
    list(result    = result,
         fit       = fit,
         design    = list(K        = K,
                          upper    = bnd$upper,
                          lower    = fut$lower,
                          binding  = bnd$binding,
                          futility = fut$label,
                          alpha    = bnd$alpha,
                          source   = bnd$source),
         events    = inp$events,
         fs        = tmpl,
         hr        = inp$hr,
         r         = inp$r,
         method    = method,
         type      = type,
         target    = target,
         scale     = scale,
         hr_method = hr_method,
         gamma_M1  = gamma_M1,
         gamma_M2  = gamma_M2),
    class = "rcpgsd_fs"
  )
}
