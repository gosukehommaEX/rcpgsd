# =============================================================================
# table_and_figure_supplement.R
#
# Rendering stage for the supplementary material of the rcpgsd article. This
# script reads the data written by data_generate_supplement.R
# (inst/sup_info/data/*.rds) and, for Figures S1 and S5 and Table S1, the
# main data (inst/scripts/data/fig1_data.rds and fig5_data.rds). It produces the
# figures, the LaTeX table and the ledger of the numbers cited in the
# supplementary material, without any computation of probabilities.
#
#   Figure S1  RCP under H1 with the adjusted and the unadjusted thresholds
#   Figure S2  Futility boundaries
#   Figure S3  Hazard reduction scale, exact and delta-method evaluation
#   Figure S4  Hazard ratio in Japan different from the other regions
#   Figure S5  RCP given efficacy stop versus the marginal RCP of
#              Teng et al. (2018)
#   Table S1   Patient-level simulation versus closed form (details)
#
# Usage
#   Open rcpgsd.Rproj (the working directory must be the package root) and
#   source this file:
#     source("inst/sup_info/table_and_figure_supplement.R")
#
# Required packages: ggplot2
#
# Outputs (relative to the package root)
#   inst/sup_info/results/FigureS1.pdf, .eps ... FigureS5.pdf, .eps
#   inst/sup_info/results/tableS1_simulation.tex
#   inst/sup_info/results/supplement_numbers.csv   Ledger of cited numbers
# The table uses the rules of the booktabs package (\usepackage{booktabs}).
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "rcpgsd")) {
  stop("Run this script from the root of the rcpgsd package ",
       "(open rcpgsd.Rproj first).", call. = FALSE)
}

library(ggplot2)

# --------------------------------------------------------------------------- #
# Read the data
# --------------------------------------------------------------------------- #
main_dir    <- file.path("inst", "scripts", "data")
data_dir    <- file.path("inst", "sup_info", "data")
results_dir <- file.path("inst", "sup_info", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

read_data <- function(dir, name, script) {
  path <- file.path(dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run %s first.", path, script),
         call. = FALSE)
  }
  readRDS(path)
}

siS1 <- read_data(data_dir, "siS1_data", "data_generate_supplement.R")
siS2 <- read_data(data_dir, "siS2_data", "data_generate_supplement.R")
siS3 <- read_data(data_dir, "siS3_data", "data_generate_supplement.R")
siS4 <- read_data(data_dir, "siS4_data", "data_generate_supplement.R")
fig1 <- read_data(main_dir, "fig1_data", "data_generate_main.R")
fig5 <- read_data(main_dir, "fig5_data", "data_generate_main.R")

SET  <- siS1$settings
MAIN <- SET$MAIN
K    <- MAIN$K

# The supplementary data must use the settings of the main data
same_settings <- function(a, b) {
  isTRUE(all.equal(a[c("ALPHA", "BETA", "HR", "K", "TIMING", "GAMMA_M1",
                       "GAMMA_M2", "F1_GRID", "EVENTS")],
                   b[c("ALPHA", "BETA", "HR", "K", "TIMING", "GAMMA_M1",
                       "GAMMA_M2", "F1_GRID", "EVENTS")]))
}
if (!same_settings(MAIN, fig1$settings)) {
  stop("The supplementary data were generated with settings that differ ",
       "from the current main data; re-run data_generate_supplement.R.",
       call. = FALSE)
}

# --------------------------------------------------------------------------- #
# Plotting constants and helpers (as in table_and_figure_manuscript.R)
# --------------------------------------------------------------------------- #
BASE_SIZE <- 24
LW        <- 1.2
PS        <- 3
LKW       <- unit(2.5, "cm")
DPI       <- 800

OKABE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#000000")
SERIES_COL <- c("Analysis 1" = OKABE[6],
                "Analysis 2" = OKABE[5],
                "Analysis 3" = OKABE[3],
                "Overall"    = "black")
TYPE_LTY   <- c("Given efficacy stop" = "solid",
                "Unconditional"       = "dashed")

BOUNDARY_LAB <- c(OBF = "O'Brien-Fleming", Pocock = "Pocock")
METHOD_LAB   <- c(M1 = "Method 1", M2 = "Method 2")
X_F1  <- "Proportion of events in Japan"
Y_RCP <- "Regional consistency probability"

theme_paper <- function() {
  theme_bw(base_size = BASE_SIZE) +
    theme(strip.background = element_rect(fill = "white", colour = "black"),
          strip.text       = element_text(face = "plain"),
          panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.box       = "vertical",
          legend.key.width = LKW,
          legend.spacing.y = unit(0.1, "cm"),
          plot.margin      = margin(8, 14, 8, 8))
}

use_cairo <- isTRUE(capabilities("cairo"))
save_figure <- function(plot, name, width, height) {
  f_pdf <- file.path(results_dir, paste0(name, ".pdf"))
  f_eps <- file.path(results_dir, paste0(name, ".eps"))
  if (use_cairo) {
    ggsave(f_pdf, plot, device = grDevices::cairo_pdf, width = width,
           height = height, units = "in", dpi = DPI,
           fallback_resolution = DPI)
    ggsave(f_eps, plot, device = grDevices::cairo_ps, width = width,
           height = height, units = "in", dpi = DPI,
           fallback_resolution = DPI, onefile = FALSE)
  } else {
    ggsave(f_pdf, plot, device = "pdf", width = width, height = height,
           units = "in", dpi = DPI)
    ggsave(f_eps, plot, device = "eps", width = width, height = height,
           units = "in", dpi = DPI)
  }
  message("Figure saved: ", f_pdf, ", ", f_eps)
}

write_tex <- function(lines, filename) {
  writeLines(lines, con = filename)
  message("LaTeX table saved: ", filename)
}

fmt <- function(x, digits = 3, na = "---") {
  ifelse(is.na(x), na, sprintf(paste0("%.", digits, "f"), x))
}

# Long data frame of the RCP given efficacy stop, the unconditional
# probability and the overall RCP
rcp_long <- function(stage, overall, keys) {
  rbind(
    data.frame(stage[, keys, drop = FALSE], method = stage$method,
               series = paste("Analysis", stage$stage),
               type = "Given efficacy stop", value = stage$rcp),
    data.frame(stage[, keys, drop = FALSE], method = stage$method,
               series = paste("Analysis", stage$stage),
               type = "Unconditional", value = stage$rcp_uncond),
    data.frame(overall[, keys, drop = FALSE], method = overall$method,
               series = "Overall", type = "Given efficacy stop",
               value = overall$rcp)
  )
}

# =============================================================================
# Consistency checks (breakdowns and totals)
# =============================================================================
check_sums <- function(stage, overall, keys, label) {
  grp <- interaction(stage[, keys], drop = TRUE)
  p_eff   <- tapply(stage$p_eff, grp, sum)
  p_joint <- tapply(stage$p_joint, grp, sum)
  og <- as.character(interaction(overall[, keys], drop = TRUE))
  err <- max(abs(p_eff[og] - overall$power),
             abs(p_joint[og] - overall$p_joint),
             abs(p_joint[og] / p_eff[og] - overall$rcp))
  if (err > 1e-10) {
    stop(sprintf("%s: stagewise probabilities do not add up (error %.2e).",
                 label, err), call. = FALSE)
  }
  invisible(err)
}
check_sums(siS1$stage, siS1$overall, c("boundary", "f1", "method"), "S1")
check_sums(siS2$stage, siS2$overall, c("case", "f1", "method"), "S2")
check_sums(siS3$stage, siS3$overall, c("boundary", "hr_method", "f1",
                                       "method"), "S3")
check_sums(siS4$stage, siS4$overall, c("boundary", "hr1", "method"), "S4")

# With the hazard ratio of the manuscript in all regions, S4 reproduces
# Figure 1 at the same proportion of events
s4_ref <- siS4$stage[abs(siS4$stage$hr1 - MAIN$HR) < 1e-9, ]
f1_ref <- fig1$stage[abs(fig1$stage$f1 - SET$F1_S4) < 1e-9, ]
key4 <- function(d) paste(d$boundary, d$method, d$stage)
stopifnot(nrow(s4_ref) == nrow(f1_ref),
          max(abs(s4_ref$rcp - f1_ref$rcp[match(key4(s4_ref),
                                                key4(f1_ref))])) < 1e-8)

# =============================================================================
# Figure S1: RCP given efficacy stop under H1, adjusted and unadjusted
# thresholds
# =============================================================================
THR_LTY <- c("Adjusted threshold" = "solid",
             "Unadjusted threshold" = "dashed")
s1a <- siS1$stage
s1u <- fig1$stage
figS1_df <- rbind(
  data.frame(boundary = s1a$boundary, method = s1a$method, f1 = s1a$f1,
             series = paste("Analysis", s1a$stage),
             threshold = "Adjusted threshold", value = s1a$rcp),
  data.frame(boundary = s1u$boundary, method = s1u$method, f1 = s1u$f1,
             series = paste("Analysis", s1u$stage),
             threshold = "Unadjusted threshold", value = s1u$rcp)
)
figS1_df$boundary  <- factor(BOUNDARY_LAB[figS1_df$boundary],
                             levels = BOUNDARY_LAB)
figS1_df$method    <- factor(METHOD_LAB[figS1_df$method],
                             levels = METHOD_LAB)
figS1_df$series    <- factor(figS1_df$series, levels = names(SERIES_COL))
figS1_df$threshold <- factor(figS1_df$threshold, levels = names(THR_LTY))

FigureS1 <- ggplot(figS1_df,
                   aes(x = f1, y = value, colour = series,
                       linetype = threshold,
                       group = interaction(series, threshold))) +
  geom_line(linewidth = LW) +
  facet_grid(boundary ~ method) +
  scale_x_continuous(name = X_F1, breaks = seq(0.05, 0.30, by = 0.05)) +
  scale_y_continuous(name = Y_RCP) +
  scale_colour_manual(name = NULL, values = SERIES_COL,
                      breaks = c("Analysis 1", "Analysis 2",
                                 "Analysis 3")) +
  scale_linetype_manual(name = NULL, values = THR_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_paper()

print(FigureS1)
save_figure(FigureS1, "FigureS1", width = 16, height = 13)

# =============================================================================
# Figure S2: futility boundaries
# =============================================================================
CASE_LAB <- c(binding            = "Binding",
              nonbinding_obeyed  = "Non-binding, obeyed",
              nonbinding_ignored = "Non-binding, ignored")
figS2_df <- rcp_long(siS2$stage, siS2$overall, c("case", "f1"))
figS2_df$case   <- factor(CASE_LAB[figS2_df$case], levels = CASE_LAB)
figS2_df$method <- factor(METHOD_LAB[figS2_df$method], levels = METHOD_LAB)
figS2_df$series <- factor(figS2_df$series, levels = names(SERIES_COL))
figS2_df$type   <- factor(figS2_df$type, levels = names(TYPE_LTY))

FigureS2 <- ggplot(figS2_df,
                   aes(x = f1, y = value, colour = series, linetype = type,
                       group = interaction(series, type))) +
  geom_line(linewidth = LW) +
  facet_grid(method ~ case) +
  scale_x_continuous(name = X_F1, breaks = c(0.1, 0.2, 0.3)) +
  scale_y_continuous(name = Y_RCP) +
  scale_colour_manual(name = NULL, values = SERIES_COL) +
  scale_linetype_manual(name = NULL, values = TYPE_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_paper()

print(FigureS2)
save_figure(FigureS2, "FigureS2", width = 18, height = 13)

# =============================================================================
# Figure S3: hazard reduction scale (Method 1), exact and delta method
# =============================================================================
HRM_LTY <- c("Exact" = "solid", "Delta method" = "dashed")
s3 <- siS3$stage
figS3_df <- rbind(
  data.frame(boundary = s3$boundary, hr_method = s3$hr_method, f1 = s3$f1,
             series = paste("Analysis", s3$stage),
             quantity = "Given efficacy stop", value = s3$rcp),
  data.frame(boundary = s3$boundary, hr_method = s3$hr_method, f1 = s3$f1,
             series = paste("Analysis", s3$stage),
             quantity = "Unconditional", value = s3$rcp_uncond)
)
figS3_df$boundary  <- factor(BOUNDARY_LAB[figS3_df$boundary],
                             levels = BOUNDARY_LAB)
figS3_df$hr_method <- factor(c(exact = "Exact",
                               delta = "Delta method")[figS3_df$hr_method],
                             levels = names(HRM_LTY))
figS3_df$quantity  <- factor(figS3_df$quantity, levels = names(TYPE_LTY))
figS3_df$series    <- factor(figS3_df$series, levels = names(SERIES_COL))

FigureS3 <- ggplot(figS3_df,
                   aes(x = f1, y = value, colour = series,
                       linetype = hr_method,
                       group = interaction(series, hr_method))) +
  geom_line(linewidth = LW) +
  facet_grid(quantity ~ boundary) +
  scale_x_continuous(name = X_F1, breaks = seq(0.05, 0.30, by = 0.05)) +
  scale_y_continuous(name = Y_RCP) +
  scale_colour_manual(name = NULL, values = SERIES_COL,
                      breaks = c("Analysis 1", "Analysis 2",
                                 "Analysis 3")) +
  scale_linetype_manual(name = NULL, values = HRM_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_paper()

print(FigureS3)
save_figure(FigureS3, "FigureS3", width = 16, height = 13)

# =============================================================================
# Figure S4: hazard ratio in Japan different from the other regions
# =============================================================================
figS4_df <- rcp_long(siS4$stage, siS4$overall, c("boundary", "hr1"))
figS4_df$boundary <- factor(BOUNDARY_LAB[figS4_df$boundary],
                            levels = BOUNDARY_LAB)
figS4_df$method   <- factor(METHOD_LAB[figS4_df$method], levels = METHOD_LAB)
figS4_df$series   <- factor(figS4_df$series, levels = names(SERIES_COL))
figS4_df$type     <- factor(figS4_df$type, levels = names(TYPE_LTY))

FigureS4 <- ggplot(figS4_df,
                   aes(x = hr1, y = value, colour = series, linetype = type,
                       group = interaction(series, type))) +
  geom_vline(xintercept = MAIN$HR, colour = "grey55", linetype = "dotted",
             linewidth = LW) +
  geom_line(linewidth = LW) +
  facet_grid(boundary ~ method) +
  scale_x_continuous(name = "Hazard ratio in Japan",
                     breaks = seq(0.6, 1.0, by = 0.1)) +
  scale_y_continuous(name = Y_RCP) +
  scale_colour_manual(name = NULL, values = SERIES_COL) +
  scale_linetype_manual(name = NULL, values = TYPE_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_paper()

print(FigureS4)
save_figure(FigureS4, "FigureS4", width = 16, height = 13)

# =============================================================================
# Table S1: patient-level simulation versus closed form (OBF boundaries)
# =============================================================================
s5  <- fig5$stage
an5 <- fig5$analysis
SET5 <- fig5$settings
sim5 <- s5[s5$source == "simulation", ]
cf5  <- s5[s5$source == "formula_expected", ]
key5 <- function(d) paste(d$tau, d$p1, d$method, d$stage)

tab5_rows <- character(0)
for (tau in SET5$TAU_SIM) {
  if (length(tab5_rows) > 0L) tab5_rows <- c(tab5_rows, "\\midrule")
  for (p1 in SET5$P1_SIM) {
    for (k in seq_len(K)) {
      a <- an5[an5$tau == tau & an5$p1 == p1 & an5$stage == k, ]
      stopifnot(nrow(a) == 1L)
      cells <- c(if (p1 == SET5$P1_SIM[1L] && k == 1L) sprintf("%g", tau)
                 else "",
                 if (k == 1L) fmt(p1, 2) else "",
                 sprintf("%d", k),
                 fmt(a$events_r1_expected, 1),
                 fmt(a$p_no_est_r1, 3),
                 fmt(a$p_inf_r1, 3))
      for (m in c("M1", "M2")) {
        id <- paste(tau, p1, m, k)
        s  <- sim5[key5(sim5) == id, ]
        f  <- cf5[key5(cf5) == id, ]
        stopifnot(nrow(s) == 1L, nrow(f) == 1L)
        cells <- c(cells, fmt(s$rcp, 3), fmt(s$se_rcp, 3), fmt(f$rcp, 3))
      }
      tab5_rows <- c(tab5_rows, paste(paste(cells, collapse = " & "),
                                      "\\\\"))
    }
  }
}

tab5_caption <- paste0(
  "Patient-level simulation and closed form of the regional consistency ",
  "probability given efficacy stopping (O'Brien-Fleming-type boundaries, ",
  "$K = ", K, "$ analyses with $", paste(SET5$EVENTS$OBF, collapse = "/"),
  "$ events, hazard ratio $", sprintf("%.2f", SET5$HR), "$, $",
  SET5$N_TOTAL, "$ patients enrolled over ", SET5$T_ACCRUAL, " months, ",
  "median survival time of ", SET5$MEDIAN_C, " months in the control ",
  "group). ",
  "$\\tau$, start of enrollment in Japan (months); ",
  "$p_{1}$, proportion of Japanese patients; ",
  "$k$, analysis; ",
  "$d_{1,k}$, expected number of events in Japan; ",
  "no est., proportion of simulated trials without a Japanese estimate; ",
  "infinite, proportion with an infinite Japanese estimate; ",
  "Sim., simulated probability ($",
  formatC(SET5$NSIM, format = "d", big.mark = "{,}"),
  "$ trials per scenario) with its Monte Carlo standard error (SE); ",
  "Formula, closed form with the expected regional proportions of events. ",
  "The hazard ratios in the overall population and in each region are ",
  "estimated by ",
  c(cox = "Cox regression", peto = "the method of Peto",
    pike = "the method of Pike")[[SET5$ESTIMATOR]],
  ".")

tab5_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  paste0("\\caption{", tab5_caption, "}"),
  "\\label{tab:simulation}",
  "\\begin{tabular}{rrr rrr rrr rrr}",
  "\\toprule",
  paste(" & & & & & &", "\\multicolumn{3}{c}{Method 1} &",
        "\\multicolumn{3}{c}{Method 2} \\\\"),
  "\\cmidrule(lr){7-9}\\cmidrule(lr){10-12}",
  paste("$\\tau$ & $p_{1}$ & $k$ & $d_{1,k}$ & No est. & Infinite &",
        "Sim. & SE & Formula & Sim. & SE & Formula \\\\"),
  "\\midrule",
  tab5_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
write_tex(tab5_lines, file.path(results_dir, "tableS1_simulation.tex"))

# =============================================================================
# Figure S5: RCP given efficacy stop (proposed) versus the marginal RCP of
# Teng et al. (2018), Pr(C_k | Z_k >= u_k), which ignores that the trial
# has not stopped before analysis k
# =============================================================================
MRG_LTY <- c("Given efficacy stop at the analysis" = "solid",
             "Given crossing the boundary (Teng et al., 2018)" = "dashed")
mg <- fig1$marginal
key_mg <- function(d) paste(d$boundary, d$f1, d$method, d$stage)
s1m <- s1u[match(key_mg(mg), key_mg(s1u)), ]
# At the first analysis both probabilities are the same
stopifnot(max(abs(mg$rcp_marginal[mg$stage == 1L] -
                    s1m$rcp[mg$stage == 1L])) < 1e-6)
figS5_df <- rbind(
  data.frame(boundary = mg$boundary, method = mg$method, f1 = mg$f1,
             series = paste("Analysis", mg$stage),
             type = names(MRG_LTY)[1L], value = s1m$rcp),
  data.frame(boundary = mg$boundary, method = mg$method, f1 = mg$f1,
             series = paste("Analysis", mg$stage),
             type = names(MRG_LTY)[2L], value = mg$rcp_marginal)
)
figS5_df$boundary <- factor(BOUNDARY_LAB[figS5_df$boundary],
                            levels = BOUNDARY_LAB)
figS5_df$method   <- factor(METHOD_LAB[figS5_df$method], levels = METHOD_LAB)
figS5_df$series   <- factor(figS5_df$series, levels = names(SERIES_COL))
figS5_df$type     <- factor(figS5_df$type, levels = names(MRG_LTY))

FigureS5 <- ggplot(figS5_df,
                   aes(x = f1, y = value, colour = series, linetype = type,
                       group = interaction(series, type))) +
  geom_line(linewidth = LW) +
  facet_grid(boundary ~ method) +
  scale_x_continuous(name = X_F1, breaks = seq(0.05, 0.30, by = 0.05)) +
  scale_y_continuous(name = Y_RCP) +
  scale_colour_manual(name = NULL, values = SERIES_COL,
                      breaks = c("Analysis 1", "Analysis 2",
                                 "Analysis 3")) +
  scale_linetype_manual(name = NULL, values = MRG_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 2)) +
  theme_paper()

print(FigureS5)
save_figure(FigureS5, "FigureS5", width = 16, height = 13)

# =============================================================================
# Ledger of the numbers cited in the supplementary material:
# results/supplement_numbers.csv (same columns as manuscript_numbers.csv)
# =============================================================================
ledger <- list()
add_num <- function(key, value, digits, description, source) {
  if (key %in% names(ledger)) {
    stop(sprintf("Duplicated ledger key '%s'.", key), call. = FALSE)
  }
  if (length(value) != 1L || !is.finite(value)) {
    stop(sprintf("Ledger value '%s' must be a single finite number.", key),
         call. = FALSE)
  }
  ledger[[key]] <<- data.frame(key = key, value = value, digits = digits,
                               text = sprintf(paste0("%.", digits, "f"),
                                              value),
                               description = description, source = source,
                               stringsAsFactors = FALSE)
}

# S1: largest reduction of the RCP by the adjustment
for (b in names(BOUNDARY_LAB)) {
  for (m in c("M1", "M2")) {
    a <- s1a[s1a$boundary == b & s1a$method == m, ]
    u <- s1u[s1u$boundary == b & s1u$method == m, ]
    u <- u[match(paste(a$f1, a$stage), paste(u$f1, u$stage)), ]
    add_num(sprintf("siS1_%s_%s_maxdiff", b, m), max(u$rcp - a$rcp), 3,
            sprintf(paste0("Largest reduction of the RCP given efficacy ",
                           "stop by the adjustment, %s, %s"), b, m),
            "siS1_data$stage$rcp, fig1_data$stage$rcp")
  }
}

# S2: power and largest difference from the binding case
o2 <- siS2$overall
for (cs in names(CASE_LAB)) {
  add_num(sprintf("siS2_power_%s", cs), o2$power[o2$case == cs][1L], 3,
          sprintf("Power, %s", CASE_LAB[[cs]]), "siS2_data$overall$power")
}
st2 <- siS2$stage
for (m in c("M1", "M2")) {
  ref <- st2[st2$case == "binding" & st2$method == m, ]
  for (cs in c("nonbinding_obeyed", "nonbinding_ignored")) {
    d <- st2[st2$case == cs & st2$method == m, ]
    d <- d[match(paste(ref$f1, ref$stage), paste(d$f1, d$stage)), ]
    add_num(sprintf("siS2_%s_%s_maxdiff", cs, m), max(abs(d$rcp - ref$rcp)),
            3, sprintf(paste0("Largest |RCP difference| from the binding ",
                              "case, %s, %s"), CASE_LAB[[cs]], m),
            "siS2_data$stage$rcp")
  }
}

# S3: largest difference between the exact and the delta-method evaluation
for (b in names(BOUNDARY_LAB)) {
  ex <- s3[s3$boundary == b & s3$hr_method == "exact", ]
  de <- s3[s3$boundary == b & s3$hr_method == "delta", ]
  de <- de[match(paste(ex$f1, ex$stage), paste(de$f1, de$stage)), ]
  add_num(sprintf("siS3_%s_maxdiff_cond", b), max(abs(ex$rcp - de$rcp)), 3,
          sprintf("Largest |exact - delta| of the RCP given efficacy stop, %s",
                  b), "siS3_data$stage$rcp")
  add_num(sprintf("siS3_%s_maxdiff_uncond", b),
          max(abs(ex$rcp_uncond - de$rcp_uncond)), 3,
          sprintf("Largest |exact - delta| of the unconditional RCP, %s", b),
          "siS3_data$stage$rcp_uncond")
}

# S5: largest overestimation by the marginal RCP of Teng et al. (2018)
for (b in names(BOUNDARY_LAB)) {
  for (m in c("M1", "M2")) {
    for (k in 2:K) {
      sel <- mg$boundary == b & mg$method == m & mg$stage == k
      add_num(sprintf("siS5_%s_%s_k%d_maxdiff", b, m, k),
              max(mg$rcp_marginal[sel] - s1m$rcp[sel]), 3,
              sprintf(paste0("Largest (marginal - given efficacy stop) RCP ",
                             "over f1, analysis %d, %s, %s"), k, b, m),
              "fig1_data$marginal, fig1_data$stage")
    }
  }
}

# S4: no treatment effect in Japan (hazard ratio 1)
o4 <- siS4$overall
for (b in names(BOUNDARY_LAB)) {
  for (m in c("M1", "M2")) {
    sel <- o4$boundary == b & o4$method == m & abs(o4$hr1 - 1) < 1e-9
    add_num(sprintf("siS4_%s_%s_overall_hr1", b, m), o4$rcp[sel], 3,
            sprintf("Overall RCP given success, HR in Japan 1, %s, %s", b,
                    m), "siS4_data$overall$rcp")
  }
  add_num(sprintf("siS4_%s_power_hr1", b),
          o4$power[o4$boundary == b & abs(o4$hr1 - 1) < 1e-9][1L], 3,
          sprintf("Power, HR in Japan 1, %s", b), "siS4_data$overall$power")
}

ledger_df <- do.call(rbind, ledger)
rownames(ledger_df) <- NULL
utils::write.csv(ledger_df, file.path(results_dir, "supplement_numbers.csv"),
                 row.names = FALSE)
message("Ledger saved: ", file.path(results_dir, "supplement_numbers.csv"),
        " (", nrow(ledger_df), " numbers)")
