# =============================================================================
# data_generate_main.R
#
# Computation stage for the manuscript figures and tables of the rcpgsd
# article (regional consistency probability for multi-regional clinical
# trials with a time-to-event endpoint and a group sequential design).
#
# This script performs all numerical computation (closed-form RCPs, adjusted
# thresholds, expected regional event proportions under delayed regional
# starts, required proportions of events in the region of interest, and the
# patient-level simulation) and writes the results to inst/scripts/data/.
# The companion script table_and_figure_manuscript.R reads these files and
# produces the figures and the LaTeX tables without any computation.
#
# Usage
#   Open rcpgsd.Rproj (the working directory must be the package root), load
#   the package, and source this file:
#     devtools::load_all()
#     source("inst/scripts/data_generate_main.R")
#   Parts of the computation can be skipped with the logical vector RUN below;
#   the data files of the skipped parts are left unchanged.
#
# Required packages: rcpgsd, rpact
#
# Running time (rough guide; the time of each step is written to the log)
#   fig1, fig2, delay : a few minutes each
#   sim               : dominates; 8 scenarios with NSIM = 100000 trials
#                       each. Set NSIM = 10000 for a quick check.
#
# Outputs (relative to the package root)
#   inst/scripts/data/fig1_data.rds    Figure 1 (RCP under H1, no delay)
#   inst/scripts/data/fig2_data.rds    Figure 2 (null RCP, adjusted thresholds)
#   inst/scripts/data/fig3_data.rds    Figure 3 (delayed start in Japan)
#   inst/scripts/data/fig4_data.rds    Figure 4 (required proportion of events)
#   inst/scripts/data/fig5_data.rds    Figure 5 (patient-level simulation)
#   inst/scripts/data/table1_data.rds  Table 1 (application, by analysis)
#   inst/scripts/data/table2_data.rds  Table 2 (application, required events)
#   inst/scripts/data/csv/*.csv        Data frames of the .rds files as CSV
#   inst/scripts/data/data_generate_main.log
#   inst/scripts/data/sessionInfo_data_generate_main.txt
# =============================================================================

# --------------------------------------------------------------------------- #
# Working directory and packages
# --------------------------------------------------------------------------- #
if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "rcpgsd")) {
  stop("Run this script from the root of the rcpgsd package ",
       "(open rcpgsd.Rproj first).", call. = FALSE)
}
# A version loaded by devtools::load_all() is used if present; otherwise the
# installed package is loaded.
if (!requireNamespace("rcpgsd", quietly = TRUE)) {
  stop("The rcpgsd package is neither loaded nor installed.", call. = FALSE)
}
if (!requireNamespace("rpact", quietly = TRUE)) {
  stop("The rpact package is required: install.packages(\"rpact\").",
       call. = FALSE)
}

# Parts of the computation to run. RUN can also be assigned before this
# file is sourced to rerun only some parts, e.g.
#   RUN <- c(fig1 = FALSE, fig2 = FALSE, delay = TRUE, sim = FALSE)
if (!exists("RUN", inherits = FALSE)) {
  RUN <- c(fig1  = TRUE,   # Figure 1
           fig2  = TRUE,   # Figure 2
           delay = TRUE,   # Figures 3 and 4, Tables 1 and 2
           sim   = TRUE)   # Figure 5 (patient-level simulation, slow)
}
if (!is.logical(RUN) ||
    !setequal(names(RUN), c("fig1", "fig2", "delay", "sim"))) {
  stop("RUN must be a logical vector with the elements fig1, fig2, delay ",
       "and sim.", call. = FALSE)
}

# --------------------------------------------------------------------------- #
# Output location and logging
# --------------------------------------------------------------------------- #
data_dir <- file.path("inst", "scripts", "data")
csv_dir  <- file.path(data_dir, "csv")
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_main.log")
cat(sprintf("data_generate_main.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}

# Evaluate an expression, log its running time, and write its warnings to
# the log instead of the console
run_step <- function(label, expr) {
  t0  <- proc.time()[["elapsed"]]
  val <- withCallingHandlers(expr, warning = function(w) {
    log_msg("  [warning] %s: %s", label, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  log_msg("%-44s %8.1f sec", label, proc.time()[["elapsed"]] - t0)
  val
}

# Save a list as .rds and each of its data frame elements as .csv
save_data <- function(obj, name) {
  obj$settings <- SETTINGS
  obj$meta     <- META
  saveRDS(obj, file.path(data_dir, paste0(name, ".rds")))
  for (el in names(obj)) {
    if (is.data.frame(obj[[el]])) {
      utils::write.csv(obj[[el]],
                       file.path(csv_dir, paste0(name, "_", el, ".csv")),
                       row.names = FALSE)
    }
  }
  log_msg("Saved %s", file.path(data_dir, paste0(name, ".rds")))
}

# Short commit hash and modified tracked files of the package repository
git_info <- function() {
  run_git <- function(args) {
    out <- tryCatch(suppressWarnings(system2("git", args, stdout = TRUE,
                                             stderr = FALSE)),
                    error = function(e) NULL)
    if (is.null(out) || !is.null(attr(out, "status"))) NULL else out
  }
  hash  <- run_git(c("rev-parse", "--short", "HEAD"))
  dirty <- run_git(c("status", "--porcelain", "--untracked-files=no"))
  list(commit   = if (length(hash) == 1L) hash else NA_character_,
       modified = if (is.null(dirty)) character(0) else trimws(dirty))
}

GIT  <- git_info()
META <- list(date           = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
             rcpgsd_version = as.character(getNamespaceVersion("rcpgsd")),
             rcpgsd_path    = getNamespaceInfo("rcpgsd", "path"),
             rpact_version  = as.character(utils::packageVersion("rpact")),
             git_commit     = GIT$commit,
             git_modified   = GIT$modified,
             R_version      = R.version.string)

log_msg("rcpgsd %s loaded from %s", META$rcpgsd_version, META$rcpgsd_path)
log_msg("git commit %s; modified tracked files: %s", META$git_commit,
        if (length(META$git_modified) == 0L) "none" else
          paste(META$git_modified, collapse = "; "))

# =============================================================================
# Common settings
# =============================================================================
ALPHA    <- 0.025          # one-sided significance level
BETA     <- 0.1            # type II error rate (power 0.9)
HR       <- 0.7            # true hazard ratio in all regions
R_ALLOC  <- 1              # allocation ratio (experimental : control)
K        <- 3L             # number of analyses
TIMING   <- c(1, 2, 3) / 3 # planned information fractions
GAMMA_M1 <- 0.5            # unadjusted threshold, Method 1
GAMMA_M2 <- 0              # unadjusted threshold, Method 2
TARGET   <- 0.8            # target value for the required proportion
F1_GRID  <- seq(0.05, 0.30, by = 0.01)   # proportion of events in Japan

# Accrual and event-time setting
N_TOTAL   <- 600                     # total number of patients
T_ACCRUAL <- 24                      # accrual period (months)
RATE      <- N_TOTAL / T_ACCRUAL     # overall accrual rate (per month)
N_REGION  <- c(60, 270, 270)         # patients by region (Japan first)
MEDIAN_C  <- 12                      # median survival, control (months)
LAMBDA_C  <- log(2) / MEDIAN_C       # control hazard

# Fixed-design number of events (Schoenfeld)
D_FIX <- (1 + R_ALLOC)^2 / R_ALLOC *
  (qnorm(1 - ALPHA) + qnorm(1 - BETA))^2 / log(HR)^2

# Overall design with rpact: Lan-DeMets spending functions (efficacy only)
# and the numbers of events required for power 1 - BETA (Schoenfeld, with
# the inflation of the group sequential design). With the number of
# patients and the accrual period given, rpact determines the follow-up
# and the expected analysis times.
make_design <- function(type) {
  rpact::getDesignGroupSequential(kMax = K, alpha = ALPHA, beta = BETA,
                                  sided = 1, typeOfDesign = type,
                                  informationRates = TIMING)
}
sample_size <- function(design = NULL) {
  args <- list(lambda2 = LAMBDA_C, hazardRatio = HR,
               allocationRatioPlanned = R_ALLOC,
               accrualTime = c(0, T_ACCRUAL),
               maxNumberOfSubjects = N_TOTAL,
               typeOfComputation = "Schoenfeld")
  if (is.null(design)) {
    do.call(rpact::getSampleSizeSurvival,
            c(list(alpha = ALPHA, beta = BETA, sided = 1), args))
  } else {
    do.call(rpact::getSampleSizeSurvival, c(list(design = design), args))
  }
}
DESIGNS     <- list(OBF    = make_design("asOF"),
                    Pocock = make_design("asP"))
SAMPLE_SIZE <- lapply(DESIGNS, sample_size)

# The fixed-design number of events of rpact must equal the formula
SS_FIX <- sample_size()
if (abs(SS_FIX$maxNumberOfEvents[1L] - D_FIX) > 1e-6) {
  stop("The fixed-design number of events of rpact differs from the ",
       "Schoenfeld formula.", call. = FALSE)
}

# Planned numbers of events: final number rounded up, interim numbers
# rounded to the planned information fractions
plan_events <- function(ss) {
  d_K <- ceiling(ss$maxNumberOfEvents[1L])
  c(round(TIMING[-K] * d_K), d_K)
}
EVENTS <- lapply(SAMPLE_SIZE, plan_events)

# Fixed design with the fixed-design number of events
DES_FIX <- list(upper = qnorm(1 - ALPHA), alpha = ALPHA)
EV_FIX  <- ceiling(D_FIX)

# Regional proportions: Japan (region 1) and two other regions of equal size
fs_const <- function(f1) c(f1, (1 - f1) / 2, (1 - f1) / 2)
# Delays in Japan. Even values give integer numbers of patients without
# rounding, so the expected analysis times equal those without delay.
TAU_GRID  <- seq(0, 14, by = 2)
TAU_LIMIT <- (1 - N_REGION[1L] / N_TOTAL) * T_ACCRUAL   # feasibility limit
# Step of the delays used to find the longest delay at which the stagewise
# target is attainable. Multiples of 0.08 months give integer numbers of
# patients in every accrual interval (12.5 patients per month in each of
# the other regions before the start in Japan).
TAU_STEP  <- 0.08
APP_TAU   <- 6                       # delay in the application (months)
TYPES     <- c("overall", "stagewise", "unconditional")

# Patient-level simulation (Figure 5)
P1_SIM  <- c(0.05, 0.10, 0.15, 0.20) # proportion of Japanese patients
TAU_SIM <- c(0, APP_TAU)             # delays in Japan (months)
NSIM      <- 100000                  # simulated trials per scenario
SEED      <- 2026
ESTIMATOR <- "cox"                   # HR estimator, overall and by region

stopifnot(APP_TAU %in% TAU_GRID,
          sum(N_REGION) == N_TOTAL,
          all(round(P1_SIM * N_TOTAL) %% 2 == 0))

design_info <- function(name) {
  des <- DESIGNS[[name]]
  ss  <- SAMPLE_SIZE[[name]]
  data.frame(boundary   = name,
             stage      = seq_len(K),
             timing     = des$informationRates,
             n_I        = as.numeric(ss$eventsPerStage[, 1L]),
             time_rpact = as.numeric(ss$analysisTime[, 1L]),
             events     = EVENTS[[name]],
             info_frac  = EVENTS[[name]] / EVENTS[[name]][K],
             upper      = des$criticalValues,
             stringsAsFactors = FALSE)
}
DESIGN_INFO <- do.call(rbind, lapply(names(DESIGNS), design_info))

# Independent check of rcpgsd against rpact for the design of the
# application (OBF, no delay): expected analysis times (rpact uses the
# unrounded interim numbers of events) and power
pw_rpact <- rpact::getPowerSurvival(DESIGNS$OBF, lambda2 = LAMBDA_C,
                                    hazardRatio = HR, directionUpper = FALSE,
                                    allocationRatioPlanned = R_ALLOC,
                                    accrualTime = c(0, T_ACCRUAL),
                                    maxNumberOfSubjects = N_TOTAL,
                                    maxNumberOfEvents = EVENTS$OBF[K])
acc_chk <- rcpgsd::rcpgsd_accrual(n_region = N_REGION,
                                  accrual_time = c(0, T_ACCRUAL),
                                  rate_overall = RATE)
ev_chk  <- rcpgsd::rcpgsd_events(events = EVENTS$OBF, hr = HR,
                                 lambda_c = LAMBDA_C,
                                 accrual_time = acc_chk$accrual_time,
                                 n_accrual = acc_chk$n_accrual, r = R_ALLOC)
pw_chk  <- rcpgsd::rcpgsd(DESIGNS$OBF, events = EVENTS$OBF,
                          fs = fs_const(N_REGION[1L] / N_TOTAL), hr = HR,
                          r = R_ALLOC, method = "M1")$overall$power
DESIGN_CHECK <- data.frame(
  quantity = c(sprintf("events at analysis %d", seq_len(K)),
               sprintf("time of analysis %d", seq_len(K)), "power"),
  rpact    = c(as.numeric(pw_rpact$eventsPerStage[, 1L]),
               as.numeric(pw_rpact$analysisTime[, 1L]),
               pw_rpact$overallReject[1L]),
  rcpgsd   = c(EVENTS$OBF, ev_chk$analysis$time, pw_chk),
  stringsAsFactors = FALSE)
DESIGN_CHECK$difference <- DESIGN_CHECK$rcpgsd - DESIGN_CHECK$rpact
# The final analysis (same number of events) must agree closely; the power
# differs only through the rounding of the interim numbers of events
chk_time  <- abs(DESIGN_CHECK$difference[DESIGN_CHECK$quantity ==
                                           sprintf("time of analysis %d", K)])
chk_power <- abs(DESIGN_CHECK$difference[DESIGN_CHECK$quantity == "power"])
if (chk_time > 1e-4 || chk_power > 0.005) {
  stop(sprintf(paste0("rcpgsd and rpact disagree (final analysis time ",
                      "difference %.2e, power difference %.4f)."),
               chk_time, chk_power), call. = FALSE)
}

SETTINGS <- list(ALPHA = ALPHA, BETA = BETA, HR = HR, R_ALLOC = R_ALLOC,
                 K = K, TIMING = TIMING, GAMMA_M1 = GAMMA_M1,
                 GAMMA_M2 = GAMMA_M2, TARGET = TARGET, F1_GRID = F1_GRID,
                 D_FIX = D_FIX, EV_FIX = EV_FIX, EVENTS = EVENTS,
                 DESIGN_INFO = DESIGN_INFO, DESIGN_CHECK = DESIGN_CHECK,
                 N_TOTAL = N_TOTAL, T_ACCRUAL = T_ACCRUAL, RATE = RATE,
                 N_REGION = N_REGION, MEDIAN_C = MEDIAN_C,
                 LAMBDA_C = LAMBDA_C, TAU_GRID = TAU_GRID,
                 TAU_LIMIT = TAU_LIMIT, TAU_STEP = TAU_STEP,
                 APP_TAU = APP_TAU, TYPES = TYPES,
                 P1_SIM = P1_SIM, TAU_SIM = TAU_SIM, NSIM = NSIM,
                 SEED = SEED, ESTIMATOR = ESTIMATOR)

log_msg("Settings: alpha = %.3f, power = %.2f, HR = %.2f, K = %d, D_FIX = %.2f",
        ALPHA, 1 - BETA, HR, K, D_FIX)
for (b in names(EVENTS)) {
  log_msg("  %-6s events = %s, upper = %s", b,
          paste(EVENTS[[b]], collapse = "/"),
          paste(sprintf("%.4f", DESIGNS[[b]]$criticalValues), collapse = "/"))
}
log_msg("  Check against rpact (OBF, no delay; rcpgsd - rpact):")
for (i in seq_len(nrow(DESIGN_CHECK))) {
  log_msg("    %-22s rpact %10.4f  rcpgsd %10.4f  difference %9.2e",
          DESIGN_CHECK$quantity[i], DESIGN_CHECK$rpact[i],
          DESIGN_CHECK$rcpgsd[i], DESIGN_CHECK$difference[i])
}
log_msg("  NSIM = %d, SEED = %d, ESTIMATOR = %s, RUN = %s", NSIM, SEED,
        ESTIMATOR, paste(names(RUN)[RUN], collapse = ", "))

# --------------------------------------------------------------------------- #
# Compute-side helpers
# --------------------------------------------------------------------------- #
rcp_h1 <- function(design, events, fs, gamma_M1 = GAMMA_M1,
                   gamma_M2 = GAMMA_M2) {
  rcpgsd::rcpgsd(design, events = events, fs = fs, hr = HR, r = R_ALLOC,
                 gamma_M1 = gamma_M1, gamma_M2 = gamma_M2)
}

# Accrual by region (Japan starting at month tau) under the overall accrual
# plan, and the expected analysis times and regional event proportions
make_delay <- function(n_region, tau) {
  acc <- rcpgsd::rcpgsd_accrual(
    n_region     = n_region,
    start        = c(tau, rep(0, length(n_region) - 1L)),
    accrual_time = c(0, T_ACCRUAL),
    rate_overall = RATE
  )
  ev <- rcpgsd::rcpgsd_events(events = EVENTS$OBF, hr = HR,
                              lambda_c = LAMBDA_C,
                              accrual_time = acc$accrual_time,
                              n_accrual = acc$n_accrual, r = R_ALLOC)
  list(acc = acc, ev = ev,
       max_round = max(abs(acc$n_accrual - acc$n_accrual_exact)))
}

# Search interval for the required proportion of events in Japan at the
# final analysis. Under the overall accrual plan, the expected numbers of
# events in Japan are proportional to n_1 while the analysis times do not
# change, so a proportion x corresponds to n_1 = x n_1^0 / f_{1,K}^0, where
# n_1^0 and f_{1,K}^0 are the planned values. Condition (4) of the article,
# n_1 <= (1 - tau / T) N, then bounds x by f_{1,K}^0 (1 - tau / T) N / n_1^0.
fs_interval <- function(fs, tau) {
  c(0.01, min(0.99, fs[1L, ncol(fs)] * (1 - tau / T_ACCRUAL) * N_TOTAL /
                N_REGION[1L]))
}

# Whether the stagewise target is attainable for a delay tau in Japan
stagewise_ok <- function(tau, m) {
  d <- make_delay(N_REGION, tau)
  if (d$max_round > 1e-8) {
    stop(sprintf("The accrual by region is not integer for a delay of %g.",
                 tau), call. = FALSE)
  }
  fsr <- suppressWarnings(rcpgsd::rcpgsd_get_fs(
    DESIGNS$OBF, events = EVENTS$OBF, fs = d$ev$fs, hr = HR,
    target = TARGET, r = R_ALLOC, method = m, type = "stagewise",
    gamma_M1 = GAMMA_M1, gamma_M2 = GAMMA_M2,
    interval = fs_interval(d$ev$fs, tau)))
  !is.na(fsr$result$fs_r1)
}

# Long data frame with one row per method and analysis
stage_df <- function(res, ...) {
  data.frame(..., res$stage, stringsAsFactors = FALSE)
}

# Marginal RCP of Teng et al. (2018) at each analysis: the fixed-design
# probability Pr(C_k | Z_k >= u_k), which conditions on crossing the
# boundary at analysis k but not on the trial having continued to analysis
# k. It is computed as a single-stage design with critical value u_k and
# the data of analysis k.
marginal_rcp <- function(design, events, fs) {
  fs <- matrix(fs, nrow = if (is.matrix(fs)) nrow(fs) else length(fs),
               ncol = length(events))
  do.call(rbind, lapply(seq_along(events), function(k) {
    res <- rcp_h1(list(upper = design$criticalValues[k], alpha = ALPHA),
                  events[k], fs[, k])
    data.frame(method = res$stage$method, stage = k,
               upper = res$stage$upper, p_cross = res$stage$p_eff,
               rcp_marginal = res$stage$rcp, stringsAsFactors = FALSE)
  }))
}

# =============================================================================
# Figure 1 data: RCP under H1 without delay (closed form), by boundary and f1
# =============================================================================
if (RUN[["fig1"]]) {
  fig1 <- run_step("Figure 1 (RCP under H1)", {
    st <- list()
    ov <- list()
    mg <- list()
    for (b in names(DESIGNS)) {
      for (f1 in F1_GRID) {
        res <- rcp_h1(DESIGNS[[b]], EVENTS[[b]], fs_const(f1))
        st[[length(st) + 1L]] <- stage_df(res, boundary = b, f1 = f1)
        ov[[length(ov) + 1L]] <- data.frame(boundary = b, f1 = f1,
                                            res$overall)
        mg[[length(mg) + 1L]] <- data.frame(
          boundary = b, f1 = f1,
          marginal_rcp(DESIGNS[[b]], EVENTS[[b]], fs_const(f1)))
      }
    }
    fx <- lapply(F1_GRID, function(f1) {
      res <- rcp_h1(DES_FIX, EV_FIX, fs_const(f1))
      data.frame(f1 = f1, method = res$stage$method,
                 events = res$stage$events, power = res$stage$p_eff,
                 rcp = res$stage$rcp, rcp_uncond = res$stage$rcp_uncond,
                 stringsAsFactors = FALSE)
    })
    list(stage    = do.call(rbind, st),
         overall  = do.call(rbind, ov),
         fixed    = do.call(rbind, fx),
         marginal = do.call(rbind, mg))
  })
  save_data(fig1, "fig1_data")
}

# =============================================================================
# Figure 2 data: null RCP and adjusted thresholds (stagewise), by boundary
# and f1
# =============================================================================
if (RUN[["fig2"]]) {
  fig2 <- run_step("Figure 2 (null RCP, adjusted thresholds)", {
    st  <- list()
    ref <- list()
    ov  <- list()
    for (b in names(DESIGNS)) {
      for (f1 in F1_GRID) {
        g <- rcpgsd::rcpgsd_get_gamma(DESIGNS[[b]], events = EVENTS[[b]],
                                      fs = fs_const(f1), r = R_ALLOC,
                                      gamma_M1 = GAMMA_M1,
                                      gamma_M2 = GAMMA_M2,
                                      type = "stagewise")
        st[[length(st) + 1L]] <- data.frame(boundary = b, f1 = f1,
                                            g$result)
        ref[[length(ref) + 1L]] <- data.frame(boundary = b, f1 = f1,
                                              method = names(g$ref),
                                              ref = unname(g$ref),
                                              stringsAsFactors = FALSE)
        ov[[length(ov) + 1L]] <- data.frame(boundary = b, f1 = f1,
                                            g$overall)
      }
    }
    list(stage   = do.call(rbind, st),
         ref     = do.call(rbind, ref),
         overall = do.call(rbind, ov))
  })
  save_data(fig2, "fig2_data")
}

# =============================================================================
# Figures 3 and 4, Tables 1 and 2: delayed start in Japan (OBF boundaries)
# =============================================================================
if (RUN[["delay"]]) {
  delay <- run_step("Delay: accrual and expected events", {
    out <- lapply(TAU_GRID, function(tau) make_delay(N_REGION, tau))
    names(out) <- paste0("tau_", TAU_GRID)
    out
  })
  max_round <- max(vapply(delay, `[[`, numeric(1L), "max_round"))
  if (max_round > 1e-8) {
    stop("The accrual by region is not integer for some delays; choose ",
         "TAU_GRID so that no rounding is needed.", call. = FALSE)
  }

  # Feasibility limit of the delay: TAU_LIMIT is feasible, a longer delay is
  # not
  feas_at_limit <- !inherits(try(rcpgsd::rcpgsd_accrual(
    n_region = N_REGION, start = c(TAU_LIMIT, 0, 0),
    accrual_time = c(0, T_ACCRUAL), rate_overall = RATE), silent = TRUE),
    "try-error")
  feas_beyond <- !inherits(try(rcpgsd::rcpgsd_accrual(
    n_region = N_REGION, start = c(TAU_LIMIT + 0.1, 0, 0),
    accrual_time = c(0, T_ACCRUAL), rate_overall = RATE), silent = TRUE),
    "try-error")
  log_msg("Feasibility: tau = %.2f feasible = %s; tau = %.2f feasible = %s",
          TAU_LIMIT, feas_at_limit, TAU_LIMIT + 0.1, feas_beyond)

  # ---- Figure 3: expected times, regional proportions and RCP by delay ----
  fig3 <- run_step("Figure 3 (RCP by delay)", {
    an <- list()
    st <- list()
    ov <- list()
    for (i in seq_along(TAU_GRID)) {
      tau <- TAU_GRID[i]
      d   <- delay[[i]]
      an[[i]] <- data.frame(tau = tau, d$ev$analysis,
                            events_r1 = d$ev$events_region[1L, ],
                            enrolled_r1 = d$ev$enrolled_region[1L, ],
                            f1 = d$ev$fs[1L, ])
      res <- rcp_h1(DESIGNS$OBF, EVENTS$OBF, d$ev$fs)
      st[[i]] <- stage_df(res, tau = tau)
      ov[[i]] <- data.frame(tau = tau, res$overall)
    }
    an <- do.call(rbind, an)
    t0 <- an$time[an$tau == 0]
    an$time_diff <- an$time - t0[an$stage]
    list(analysis = an,
         stage    = do.call(rbind, st),
         overall  = do.call(rbind, ov),
         feasibility = data.frame(tau_limit = TAU_LIMIT,
                                  feasible_at_limit = feas_at_limit,
                                  feasible_beyond = feas_beyond))
  })
  log_msg("  max |analysis time - time without delay| = %.2e months",
          max(abs(fig3$analysis$time_diff)))
  save_data(fig3, "fig3_data")

  # ---- Figure 4: required proportion of events in Japan by delay ----
  fig4 <- run_step("Figure 4 (required proportion by delay)", {
    rq <- list()
    for (i in seq_along(TAU_GRID)) {
      for (ty in TYPES) {
        fsr <- rcpgsd::rcpgsd_get_fs(DESIGNS$OBF, events = EVENTS$OBF,
                                     fs = delay[[i]]$ev$fs, hr = HR,
                                     target = TARGET, r = R_ALLOC,
                                     type = ty, gamma_M1 = GAMMA_M1,
                                     gamma_M2 = GAMMA_M2,
                                     interval = fs_interval(delay[[i]]$ev$fs,
                                                            TAU_GRID[i]))
        rq[[length(rq) + 1L]] <- data.frame(tau = TAU_GRID[i], type = ty,
                                            fsr$result,
                                            stringsAsFactors = FALSE)
      }
    }
    # Fixed design with the fixed-design number of events (reference)
    fsr_fix <- rcpgsd::rcpgsd_get_fs(DES_FIX, events = EV_FIX,
                                     fs = fs_const(N_REGION[1L] / N_TOTAL),
                                     hr = HR, target = TARGET, r = R_ALLOC,
                                     type = "unconditional",
                                     gamma_M1 = GAMMA_M1,
                                     gamma_M2 = GAMMA_M2)
    list(required = do.call(rbind, rq),
         fixed    = data.frame(events = EV_FIX, fsr_fix$result))
  })

  # Longest delay at which the stagewise target is attainable, by bisection
  # over the delays that are multiples of TAU_STEP. The limit is reported
  # as the midpoint of the final bracket.
  fig4$critical <- run_step("Figure 4 (longest attainable delay)", {
    do.call(rbind, lapply(c("M1", "M2"), function(m) {
      lo <- 0L
      hi <- as.integer(round(max(TAU_GRID) / TAU_STEP))
      if (!stagewise_ok(lo * TAU_STEP, m) || stagewise_ok(hi * TAU_STEP, m)) {
        stop(sprintf(paste0("The stagewise target of %s must be attainable ",
                            "without delay and not attainable at %g months."),
                     m, max(TAU_GRID)), call. = FALSE)
      }
      while (hi - lo > 1L) {
        mid <- (lo + hi) %/% 2L
        if (stagewise_ok(mid * TAU_STEP, m)) lo <- mid else hi <- mid
      }
      data.frame(method             = m,
                 tau_attainable     = lo * TAU_STEP,
                 tau_not_attainable = hi * TAU_STEP,
                 tau_limit          = (lo + hi) / 2 * TAU_STEP,
                 stringsAsFactors   = FALSE)
    }))
  })
  # The delays of TAU_GRID must agree with the limit
  for (m in c("M1", "M2")) {
    rq  <- fig4$required[fig4$required$type == "stagewise" &
                           fig4$required$method == m, ]
    lim <- fig4$critical$tau_limit[fig4$critical$method == m]
    if (!identical(!is.na(rq$fs_r1), rq$tau < lim)) {
      stop(sprintf(paste0("The stagewise results of %s on TAU_GRID disagree ",
                          "with the longest attainable delay."), m),
           call. = FALSE)
    }
    log_msg("  %s: stagewise target attainable at %.2f, not at %.2f months",
            m, fig4$critical$tau_attainable[fig4$critical$method == m],
            fig4$critical$tau_not_attainable[fig4$critical$method == m])
  }
  save_data(fig4, "fig4_data")

  # ---- Table 1: application (delay APP_TAU), unadjusted and adjusted ----
  table1 <- run_step("Table 1 (application)", {
    d   <- delay[[paste0("tau_", APP_TAU)]]
    res <- rcp_h1(DESIGNS$OBF, EVENTS$OBF, d$ev$fs)
    gam <- rcpgsd::rcpgsd_get_gamma(DESIGNS$OBF, events = EVENTS$OBF,
                                    fs = d$ev$fs, r = R_ALLOC,
                                    gamma_M1 = GAMMA_M1,
                                    gamma_M2 = GAMMA_M2,
                                    type = "stagewise")
    res_adj <- rcp_h1(DESIGNS$OBF, EVENTS$OBF, d$ev$fs,
                      gamma_M1 = gam$result$gamma_M1_adj,
                      gamma_M2 = gam$result$gamma_M2_adj)
    list(analysis     = data.frame(d$ev$analysis,
                                   events_r1 = d$ev$events_region[1L, ],
                                   f1 = d$ev$fs[1L, ]),
         stage_unadj  = res$stage,
         overall_unadj = res$overall,
         stage_adj    = res_adj$stage,
         marginal     = marginal_rcp(DESIGNS$OBF, EVENTS$OBF, d$ev$fs),
         overall_adj  = res_adj$overall,
         gamma        = gam$result,
         gamma_overall = gam$overall,
         ref          = data.frame(method = names(gam$ref),
                                   ref = unname(gam$ref),
                                   stringsAsFactors = FALSE),
         accrual      = data.frame(region = seq_along(N_REGION),
                                   n_region = N_REGION,
                                   start = d$acc$start),
         n_accrual    = d$acc$n_accrual,
         accrual_time = d$acc$accrual_time)
  })
  save_data(table1, "table1_data")

  # ---- Table 2: required proportion in the application (no delay and
  #      delay APP_TAU), from the Figure 4 data ----
  table2 <- list(required = fig4$required[fig4$required$tau %in%
                                            c(0, APP_TAU), ],
                 fixed    = fig4$fixed)
  save_data(table2, "table2_data")
}

# =============================================================================
# Figure 5 data: patient-level simulation versus closed form (OBF)
# =============================================================================
if (RUN[["sim"]]) {
  fig5 <- run_step("Figure 5 (patient-level simulation)", {
    st <- list()
    ov <- list()
    an <- list()
    for (tau in TAU_SIM) {
      for (p1 in P1_SIM) {
        n1    <- round(p1 * N_TOTAL)
        n_reg <- c(n1, rep((N_TOTAL - n1) / 2, 2L))
        d     <- make_delay(n_reg, tau)
        sim   <- run_step(sprintf("  simulation tau = %g, p1 = %.2f",
                                  tau, p1),
                          rcpgsd::rcpgsd_sim(
                            DESIGNS$OBF, events = EVENTS$OBF, hr = HR,
                            lambda_c = LAMBDA_C,
                            accrual_time = d$acc$accrual_time,
                            n_accrual = d$acc$n_accrual, r = R_ALLOC,
                            gamma_M1 = GAMMA_M1, gamma_M2 = GAMMA_M2,
                            estimator = ESTIMATOR, nsim = NSIM,
                            seed = SEED))
        cf_plan <- rcp_h1(DESIGNS$OBF, EVENTS$OBF, d$ev$fs)
        cf_mean <- rcp_h1(DESIGNS$OBF, EVENTS$OBF, sim$fs_mean)
        key <- list(tau = tau, p1 = p1, n1 = n1)
        st[[length(st) + 1L]] <- rbind(
          data.frame(key, source = "simulation",
                     sim$stage, stringsAsFactors = FALSE),
          data.frame(key, source = "formula_expected",
                     cf_plan$stage, se_rcp = NA_real_,
                     stringsAsFactors = FALSE),
          data.frame(key, source = "formula_mean",
                     cf_mean$stage, se_rcp = NA_real_,
                     stringsAsFactors = FALSE))
        ov[[length(ov) + 1L]] <- rbind(
          data.frame(key, source = "simulation", sim$overall,
                     stringsAsFactors = FALSE),
          data.frame(key, source = "formula_expected", cf_plan$overall,
                     stringsAsFactors = FALSE),
          data.frame(key, source = "formula_mean", cf_mean$overall,
                     stringsAsFactors = FALSE))
        an[[length(an) + 1L]] <- data.frame(
          key, sim$analysis,
          time_expected = d$ev$analysis$time,
          events_r1_expected = d$ev$events_region[1L, ],
          f1_expected = d$ev$fs[1L, ],
          f1_mean = sim$fs_mean[1L, ])
      }
    }
    list(stage    = do.call(rbind, st),
         overall  = do.call(rbind, ov),
         analysis = do.call(rbind, an))
  })
  save_data(fig5, "fig5_data")
}

# =============================================================================
# Record the execution environment for reproducibility
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_main.txt")
)

log_msg("All requested data written to '%s/'.", data_dir)

# A RUN assigned for this run must not carry over to the next one
rm(RUN)
