# =============================================================================
# data_generate_supplement.R
#
# Computation stage for the figures of the supplementary material of the
# rcpgsd article. The results are written to inst/sup_info/data/ and are
# rendered by table_and_figure_supplement.R without any computation.
#
#   S1  RCP under H1 with the adjusted and the unadjusted thresholds
#   S2  Futility boundaries: binding, non-binding obeyed, non-binding ignored
#   S3  Hazard reduction scale (Method 1): exact and delta-method evaluation
#   S4  Hazard ratio in Japan different from the other regions
#
# The settings (significance level, power, hazard ratio, analysis timing,
# thresholds, grid of proportions) are read from the main data, so that the
# supplementary results use exactly the settings of the manuscript. The
# adjusted thresholds of S1 are those of Figure 2.
#
# Usage
#   Open rcpgsd.Rproj (the working directory must be the package root), load
#   the package, run data_generate_main.R first, and source this file:
#     devtools::load_all()
#     source("inst/sup_info/data_generate_supplement.R")
#   Parts of the computation can be skipped with the logical vector RUN below.
#
# Required packages: rcpgsd, rpact
#
# Running time (rough guide; the time of each step is written to the log)
#   S1, S2, S4 : about one minute each
#   S3         : a few minutes (numerical integration)
#
# Inputs (relative to the package root)
#   inst/scripts/data/fig1_data.rds    settings of the manuscript
#   inst/scripts/data/fig2_data.rds    adjusted thresholds (S1)
#
# Outputs (relative to the package root)
#   inst/sup_info/data/siS1_data.rds ... siS4_data.rds
#   inst/sup_info/data/csv/*.csv       Data frames of the .rds files as CSV
#   inst/sup_info/data/data_generate_supplement.log
#   inst/sup_info/data/sessionInfo_data_generate_supplement.txt
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
if (!requireNamespace("rcpgsd", quietly = TRUE)) {
  stop("The rcpgsd package is neither loaded nor installed.", call. = FALSE)
}
if (!requireNamespace("rpact", quietly = TRUE)) {
  stop("The rpact package is required: install.packages(\"rpact\").",
       call. = FALSE)
}

# Parts of the computation to run
RUN <- c(S1 = TRUE, S2 = TRUE, S3 = TRUE, S4 = TRUE)

# --------------------------------------------------------------------------- #
# Output location and logging
# --------------------------------------------------------------------------- #
main_dir <- file.path("inst", "scripts", "data")
data_dir <- file.path("inst", "sup_info", "data")
csv_dir  <- file.path(data_dir, "csv")
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_supplement.log")
cat(sprintf("data_generate_supplement.R log -- %s\n", format(Sys.time())),
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
# Settings of the manuscript
# =============================================================================
read_main <- function(name) {
  path <- file.path(main_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run data_generate_main.R first.", path),
         call. = FALSE)
  }
  readRDS(path)
}
MAIN <- read_main("fig1_data")$settings

ALPHA    <- MAIN$ALPHA
BETA     <- MAIN$BETA
HR       <- MAIN$HR
R_ALLOC  <- MAIN$R_ALLOC
K        <- MAIN$K
TIMING   <- MAIN$TIMING
GAMMA_M1 <- MAIN$GAMMA_M1
GAMMA_M2 <- MAIN$GAMMA_M2
F1_GRID  <- MAIN$F1_GRID
D_FIX    <- MAIN$D_FIX
LAMBDA_C <- MAIN$LAMBDA_C
N_TOTAL  <- MAIN$N_TOTAL
T_ACCRUAL <- MAIN$T_ACCRUAL

# Futility boundaries (S2): Hwang-Shih-DeCani beta spending
GAMMA_B <- -2
# Hazard ratio in Japan (S4); proportion of events in Japan
HR1_GRID <- seq(0.6, 1.0, by = 0.02)
F1_S4    <- 0.1

fs_const <- function(f1) c(f1, (1 - f1) / 2, (1 - f1) / 2)

# Designs and numbers of events with rpact, as in data_generate_main.R
make_design <- function(type, binding = NULL) {
  args <- list(kMax = K, alpha = ALPHA, beta = BETA, sided = 1,
               typeOfDesign = type, informationRates = TIMING)
  if (!is.null(binding)) {
    args <- c(args, list(typeBetaSpending = "bsHSD", gammaB = GAMMA_B,
                         bindingFutility = binding))
  }
  do.call(rpact::getDesignGroupSequential, args)
}
plan_events <- function(design) {
  ss <- rpact::getSampleSizeSurvival(design, lambda2 = LAMBDA_C,
                                     hazardRatio = HR,
                                     allocationRatioPlanned = R_ALLOC,
                                     accrualTime = c(0, T_ACCRUAL),
                                     maxNumberOfSubjects = N_TOTAL,
                                     typeOfComputation = "Schoenfeld")
  d_K <- ceiling(ss$maxNumberOfEvents[1L])
  structure(c(round(TIMING[-K] * d_K), d_K),
            n_I = as.numeric(ss$eventsPerStage[, 1L]))
}

# Efficacy-only designs of the manuscript
DESIGNS <- list(OBF    = make_design("asOF"),
                Pocock = make_design("asP"))
EVENTS  <- lapply(DESIGNS, plan_events)

# The designs must reproduce those of the manuscript
for (b in names(DESIGNS)) {
  di <- MAIN$DESIGN_INFO[MAIN$DESIGN_INFO$boundary == b, ]
  if (!identical(as.numeric(EVENTS[[b]]), as.numeric(MAIN$EVENTS[[b]])) ||
      max(abs(DESIGNS[[b]]$criticalValues - di$upper)) > 1e-8) {
    stop(sprintf("The %s design differs from that of the manuscript.", b),
         call. = FALSE)
  }
}

# Designs with futility boundaries (OBF efficacy boundaries, HSD beta
# spending), binding and non-binding
DES_FUT <- list(binding    = make_design("asOF", binding = TRUE),
                nonbinding = make_design("asOF", binding = FALSE))
EV_FUT  <- lapply(DES_FUT, plan_events)

fut_info <- function(name) {
  des <- DES_FUT[[name]]
  data.frame(design = name, stage = seq_len(K),
             timing = des$informationRates,
             n_I = attr(EV_FUT[[name]], "n_I"),
             events = as.numeric(EV_FUT[[name]]),
             upper = des$criticalValues,
             lower = c(des$futilityBounds, NA_real_),
             stringsAsFactors = FALSE)
}

SETTINGS <- list(MAIN = MAIN, GAMMA_B = GAMMA_B, HR1_GRID = HR1_GRID,
                 F1_S4 = F1_S4, EV_FUT = EV_FUT,
                 FUT_INFO = do.call(rbind, lapply(names(DES_FUT), fut_info)))

log_msg("Settings from %s: alpha = %.3f, power = %.2f, HR = %.2f, K = %d",
        file.path(main_dir, "fig1_data.rds"), ALPHA, 1 - BETA, HR, K)
for (b in names(EV_FUT)) {
  log_msg("  futility (%s): events = %s, upper = %s, lower = %s", b,
          paste(EV_FUT[[b]], collapse = "/"),
          paste(sprintf("%.4f", DES_FUT[[b]]$criticalValues), collapse = "/"),
          paste(sprintf("%.4f", DES_FUT[[b]]$futilityBounds), collapse = "/"))
}
log_msg("  RUN = %s", paste(names(RUN)[RUN], collapse = ", "))

stage_df <- function(res, ...) {
  data.frame(..., res$stage, stringsAsFactors = FALSE)
}
overall_df <- function(res, ...) {
  data.frame(..., res$overall, stringsAsFactors = FALSE)
}

# =============================================================================
# S1: RCP under H1 with the adjusted thresholds of Figure 2
# =============================================================================
if (RUN[["S1"]]) {
  g2 <- read_main("fig2_data")$stage
  siS1 <- run_step("S1 (adjusted thresholds under H1)", {
    st <- list()
    ov <- list()
    for (b in names(DESIGNS)) {
      for (f1 in F1_GRID) {
        g <- g2[g2$boundary == b & abs(g2$f1 - f1) < 1e-9, ]
        g <- g[order(g$stage), ]
        stopifnot(nrow(g) == K)
        res <- rcpgsd::rcpgsd(DESIGNS[[b]], events = EVENTS[[b]],
                              fs = fs_const(f1), hr = HR, r = R_ALLOC,
                              gamma_M1 = g$gamma_M1_adj,
                              gamma_M2 = g$gamma_M2_adj)
        st[[length(st) + 1L]] <- stage_df(res, boundary = b, f1 = f1)
        ov[[length(ov) + 1L]] <- overall_df(res, boundary = b, f1 = f1)
      }
    }
    list(stage = do.call(rbind, st), overall = do.call(rbind, ov))
  })
  save_data(siS1, "siS1_data")
}

# =============================================================================
# S2: futility boundaries (OBF efficacy boundaries, HSD beta spending)
# =============================================================================
if (RUN[["S2"]]) {
  cases <- list(
    list(label = "binding",            design = "binding",
         futility = "obey"),
    list(label = "nonbinding_obeyed",  design = "nonbinding",
         futility = "obey"),
    list(label = "nonbinding_ignored", design = "nonbinding",
         futility = "ignore")
  )
  siS2 <- run_step("S2 (futility boundaries)", {
    st <- list()
    ov <- list()
    for (cs in cases) {
      for (f1 in F1_GRID) {
        res <- rcpgsd::rcpgsd(DES_FUT[[cs$design]],
                              events = EV_FUT[[cs$design]],
                              fs = fs_const(f1), hr = HR, r = R_ALLOC,
                              gamma_M1 = GAMMA_M1, gamma_M2 = GAMMA_M2,
                              futility = cs$futility)
        st[[length(st) + 1L]] <- stage_df(res, case = cs$label, f1 = f1)
        ov[[length(ov) + 1L]] <- overall_df(res, case = cs$label, f1 = f1)
      }
    }
    list(stage = do.call(rbind, st), overall = do.call(rbind, ov))
  })
  save_data(siS2, "siS2_data")
}

# =============================================================================
# S3: hazard reduction scale (Method 1), exact and delta-method evaluation
# =============================================================================
if (RUN[["S3"]]) {
  siS3 <- run_step("S3 (hazard reduction scale)", {
    st <- list()
    ov <- list()
    for (b in names(DESIGNS)) {
      for (hm in c("exact", "delta")) {
        for (f1 in F1_GRID) {
          res <- rcpgsd::rcpgsd(DESIGNS[[b]], events = EVENTS[[b]],
                                fs = fs_const(f1), hr = HR, r = R_ALLOC,
                                method = "M1", scale = "hr",
                                gamma_M1 = GAMMA_M1, hr_method = hm)
          st[[length(st) + 1L]] <- stage_df(res, boundary = b,
                                            hr_method = hm, f1 = f1)
          ov[[length(ov) + 1L]] <- overall_df(res, boundary = b,
                                              hr_method = hm, f1 = f1)
        }
      }
    }
    list(stage = do.call(rbind, st), overall = do.call(rbind, ov))
  })
  save_data(siS3, "siS3_data")
}

# =============================================================================
# S4: hazard ratio in Japan different from the other regions
# =============================================================================
if (RUN[["S4"]]) {
  siS4 <- run_step("S4 (hazard ratio in Japan)", {
    st <- list()
    ov <- list()
    for (b in names(DESIGNS)) {
      for (hr1 in HR1_GRID) {
        res <- rcpgsd::rcpgsd(DESIGNS[[b]], events = EVENTS[[b]],
                              fs = fs_const(F1_S4), hr = c(hr1, HR, HR),
                              r = R_ALLOC, gamma_M1 = GAMMA_M1,
                              gamma_M2 = GAMMA_M2)
        st[[length(st) + 1L]] <- stage_df(res, boundary = b, hr1 = hr1)
        ov[[length(ov) + 1L]] <- overall_df(res, boundary = b, hr1 = hr1)
      }
    }
    list(stage = do.call(rbind, st), overall = do.call(rbind, ov))
  })
  save_data(siS4, "siS4_data")
}

# =============================================================================
# Record the execution environment for reproducibility
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_supplement.txt")
)

log_msg("All requested data written to '%s/'.", data_dir)
