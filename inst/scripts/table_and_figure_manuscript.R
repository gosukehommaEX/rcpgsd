# =============================================================================
# table_and_figure_manuscript.R
#
# Rendering stage for the manuscript figures and tables of the rcpgsd
# article. This script reads the data written by data_generate_main.R
# (inst/scripts/data/*.rds) and produces the figures, the LaTeX tables and
# the ledger of the numbers cited in the text. It performs no numerical
# computation of probabilities, so figures and tables can be adjusted
# without re-running data_generate_main.R.
#
# Usage
#   Open rcpgsd.Rproj (the working directory must be the package root) and
#   source this file:
#     source("inst/scripts/table_and_figure_manuscript.R")
#
# Prerequisite
#   inst/scripts/data/fig1_data.rds ... fig5_data.rds, table1_data.rds and
#   table2_data.rds produced by data_generate_main.R.
#
# Required packages: ggplot2
#
# Outputs (relative to the package root)
#   inst/scripts/results/Figure1.pdf, Figure1.eps   RCP under H1 (no delay)
#   inst/scripts/results/Figure2.pdf, Figure2.eps   Null RCP and thresholds
#   inst/scripts/results/Figure3.pdf, Figure3.eps   Delayed start in Japan
#   inst/scripts/results/Figure4.pdf, Figure4.eps   Required proportion
#   inst/scripts/results/Figure5.pdf, Figure5.eps   Patient-level simulation
#   inst/scripts/results/table1_application.tex     Table 1
#   inst/scripts/results/table2_required.tex        Table 2
#   inst/scripts/results/manuscript_numbers.csv     Ledger of cited numbers
#
# In the manuscript, each table is included with \input{table1_application}
# and each figure with \includegraphics{Figure1}; the captions of the
# figures are written in the manuscript, those of the tables here. The
# tables use the rules of the booktabs package (\usepackage{booktabs}).
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "rcpgsd")) {
  stop("Run this script from the root of the rcpgsd package ",
       "(open rcpgsd.Rproj first).", call. = FALSE)
}

library(ggplot2)

# --------------------------------------------------------------------------- #
# Read the data produced by data_generate_main.R
# --------------------------------------------------------------------------- #
data_dir    <- file.path("inst", "scripts", "data")
results_dir <- file.path("inst", "scripts", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

read_data <- function(name) {
  path <- file.path(data_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop(sprintf("'%s' not found; run data_generate_main.R first.", path),
         call. = FALSE)
  }
  readRDS(path)
}

fig1_data   <- read_data("fig1_data")
fig2_data   <- read_data("fig2_data")
fig3_data   <- read_data("fig3_data")
fig4_data   <- read_data("fig4_data")
fig5_data   <- read_data("fig5_data")
table1_data <- read_data("table1_data")
table2_data <- read_data("table2_data")

SET <- fig1_data$settings

# --------------------------------------------------------------------------- #
# Plotting constants
# --------------------------------------------------------------------------- #
BASE_SIZE <- 24                  # base font size
LW        <- 1.2                 # line width
PS        <- 3                   # point size
LKW       <- unit(2.5, "cm")     # legend key width
DPI       <- 800                 # resolution of the EPS output

# Okabe-Ito colour palette (colour-blind friendly)
OKABE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#000000")

SERIES_COL <- c("Analysis 1"   = OKABE[6],
                "Analysis 2"   = OKABE[5],
                "Analysis 3"   = OKABE[3],
                "Overall"      = "black",
                "Fixed design" = "grey55")
TYPE_LTY   <- c("Given efficacy stop" = "solid",
                "Unconditional"       = "dashed")
CRIT_COL   <- c("Overall"       = OKABE[5],
                "Stagewise"     = OKABE[6],
                "Unconditional" = OKABE[3])

BOUNDARY_LAB <- c(OBF = "O'Brien-Fleming", Pocock = "Pocock")
METHOD_LAB   <- c(M1 = "Method 1", M2 = "Method 2")
CRIT_LAB     <- c(overall = "Overall", stagewise = "Stagewise",
                  unconditional = "Unconditional")

X_F1  <- "Proportion of events in Japan"
X_TAU <- "Start of enrollment in Japan (months)"
Y_RCP <- "Regional consistency probability"

# --------------------------------------------------------------------------- #
# Rendering helpers
# --------------------------------------------------------------------------- #
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

# Save a figure in PDF and EPS formats. The cairo devices embed the fonts;
# no transparency is used, so the EPS file is fully vector.
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

# Fixed-point formatting; NA becomes "---"
fmt <- function(x, digits = 3, na = "---") {
  ifelse(is.na(x), na, sprintf(paste0("%.", digits, "f"), x))
}

# Value of a column for one method, ordered by analysis
by_stage <- function(df, method, column) {
  sub <- df[df$method == method, , drop = FALSE]
  sub[order(sub$stage), column]
}

# =============================================================================
# Consistency checks of the data (breakdowns and totals)
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

check_sums(fig1_data$stage, fig1_data$overall,
           c("boundary", "f1", "method"), "Figure 1")
check_sums(fig3_data$stage, fig3_data$overall, c("tau", "method"),
           "Figure 3")
check_sums(table1_data$stage_unadj, table1_data$overall_unadj, "method",
           "Table 1 (unadjusted)")
check_sums(table1_data$stage_adj, table1_data$overall_adj, "method",
           "Table 1 (adjusted)")

# Adjusted null RCP equals the reference wherever the threshold was adjusted
check_adjusted <- function(res, ref, label) {
  base <- c(M1 = SET$GAMMA_M1, M2 = SET$GAMMA_M2)
  for (m in c("M1", "M2")) {
    g   <- res[[paste0("gamma_", m, "_adj")]]
    adj <- res[[paste0("rcp0_", m, "_adj")]]
    sel <- !is.na(g) & g > base[[m]] + 1e-9
    if (any(abs(adj[sel] - ref[[m]][sel]) > 1e-5)) {
      stop(sprintf("%s: adjusted null RCP differs from the reference (%s).",
                   label, m), call. = FALSE)
    }
  }
}
ref2 <- fig2_data$ref
ref_of <- function(m) {
  r <- ref2[ref2$method == m, ]
  r$ref[match(paste(fig2_data$stage$boundary, fig2_data$stage$f1),
              paste(r$boundary, r$f1))]
}
check_adjusted(fig2_data$stage, list(M1 = ref_of("M1"), M2 = ref_of("M2")),
               "Figure 2")
ref1 <- setNames(table1_data$ref$ref, table1_data$ref$method)
check_adjusted(table1_data$gamma,
               list(M1 = rep(ref1[["M1"]], SET$K),
                    M2 = rep(ref1[["M2"]], SET$K)), "Table 1")

# The reference value does not depend on the boundaries
ref_wide <- reshape(ref2[, c("boundary", "f1", "method", "ref")],
                    idvar = c("f1", "method"), timevar = "boundary",
                    direction = "wide")
stopifnot(max(abs(ref_wide$ref.OBF - ref_wide$ref.Pocock)) < 1e-8)

# A delayed start leaves the expected analysis times unchanged
stopifnot(max(abs(fig3_data$analysis$time_diff)) < 1e-6)

# All analyses were reached in the simulation (no drop-out)
stopifnot(all(fig5_data$analysis$p_not_reached == 0))

# =============================================================================
# Figure 1: RCP under H1 without delay, by boundary and proportion of events
# =============================================================================
s1 <- fig1_data$stage
o1 <- fig1_data$overall
x1 <- fig1_data$fixed

fig1_df <- rbind(
  data.frame(boundary = s1$boundary, method = s1$method, f1 = s1$f1,
             series = paste("Analysis", s1$stage),
             type = "Given efficacy stop", value = s1$rcp),
  data.frame(boundary = s1$boundary, method = s1$method, f1 = s1$f1,
             series = paste("Analysis", s1$stage),
             type = "Unconditional", value = s1$rcp_uncond),
  data.frame(boundary = o1$boundary, method = o1$method, f1 = o1$f1,
             series = "Overall", type = "Given efficacy stop",
             value = o1$rcp),
  do.call(rbind, lapply(names(BOUNDARY_LAB), function(b) {
    data.frame(boundary = b, method = x1$method, f1 = x1$f1,
               series = "Fixed design", type = "Given efficacy stop",
               value = x1$rcp)
  }))
)
fig1_df$boundary <- factor(BOUNDARY_LAB[fig1_df$boundary],
                           levels = BOUNDARY_LAB)
fig1_df$method   <- factor(METHOD_LAB[fig1_df$method], levels = METHOD_LAB)
fig1_df$series   <- factor(fig1_df$series, levels = names(SERIES_COL))
fig1_df$type     <- factor(fig1_df$type, levels = names(TYPE_LTY))

Figure1 <- ggplot(fig1_df,
                  aes(x = f1, y = value, colour = series, linetype = type,
                      group = interaction(series, type))) +
  geom_line(linewidth = LW) +
  facet_grid(boundary ~ method) +
  scale_x_continuous(name = X_F1, breaks = seq(0.05, 0.30, by = 0.05)) +
  scale_y_continuous(name = Y_RCP, breaks = seq(0.5, 1, by = 0.1)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_colour_manual(name = NULL, values = SERIES_COL) +
  scale_linetype_manual(name = NULL, values = TYPE_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_paper()

print(Figure1)
save_figure(Figure1, "Figure1", width = 16, height = 13)

# =============================================================================
# Figure 2: null RCP with the unadjusted thresholds and adjusted thresholds
# =============================================================================
s2 <- fig2_data$stage
QUANTITY_LAB <- c(null  = "Null RCP (unadjusted threshold)",
                  gamma = "Adjusted threshold")

fig2_df <- do.call(rbind, lapply(c("M1", "M2"), function(m) {
  rbind(
    data.frame(boundary = s2$boundary, method = m, f1 = s2$f1,
               series = paste("Analysis", s2$stage), quantity = "null",
               value = s2[[paste0("rcp0_", m, "_unadj")]]),
    data.frame(boundary = s2$boundary, method = m, f1 = s2$f1,
               series = paste("Analysis", s2$stage), quantity = "gamma",
               value = s2[[paste0("gamma_", m, "_adj")]])
  )
}))
fig2_df <- rbind(fig2_df,
                 data.frame(boundary = ref2$boundary, method = ref2$method,
                            f1 = ref2$f1, series = "Fixed design",
                            quantity = "null", value = ref2$ref))
fig2_df$boundary <- factor(BOUNDARY_LAB[fig2_df$boundary],
                           levels = BOUNDARY_LAB)
fig2_df$method   <- factor(METHOD_LAB[fig2_df$method], levels = METHOD_LAB)
fig2_df$quantity <- factor(QUANTITY_LAB[fig2_df$quantity],
                           levels = QUANTITY_LAB)
fig2_df$series   <- factor(fig2_df$series, levels = names(SERIES_COL))
fig2_df <- fig2_df[!is.na(fig2_df$value), ]

Figure2 <- ggplot(fig2_df,
                  aes(x = f1, y = value, colour = series, group = series)) +
  geom_line(linewidth = LW) +
  facet_grid(quantity ~ method + boundary, scales = "free_y") +
  scale_x_continuous(name = X_F1, breaks = c(0.1, 0.2, 0.3)) +
  scale_y_continuous(name = NULL) +
  scale_colour_manual(name = NULL, values = SERIES_COL,
                      breaks = c("Analysis 1", "Analysis 2", "Analysis 3",
                                 "Fixed design"),
                      labels = c("Analysis 1", "Analysis 2", "Analysis 3",
                                 "Fixed-design reference")) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_paper()

print(Figure2)
save_figure(Figure2, "Figure2", width = 18, height = 11)

# =============================================================================
# Figure 3: delayed start in Japan (OBF boundaries, application setting)
# =============================================================================
a3 <- fig3_data$analysis
s3 <- fig3_data$stage
o3 <- fig3_data$overall
PANEL3_LAB <- c(f  = "Proportion of events in Japan",
                M1 = "RCP, Method 1",
                M2 = "RCP, Method 2")

fig3_df <- rbind(
  data.frame(panel = "f", tau = a3$tau, series = paste("Analysis", a3$stage),
             type = "Given efficacy stop", value = a3$f1),
  data.frame(panel = s3$method, tau = s3$tau,
             series = paste("Analysis", s3$stage),
             type = "Given efficacy stop", value = s3$rcp),
  data.frame(panel = s3$method, tau = s3$tau,
             series = paste("Analysis", s3$stage),
             type = "Unconditional", value = s3$rcp_uncond),
  data.frame(panel = o3$method, tau = o3$tau, series = "Overall",
             type = "Given efficacy stop", value = o3$rcp)
)
fig3_df$panel  <- factor(PANEL3_LAB[fig3_df$panel], levels = PANEL3_LAB)
fig3_df$series <- factor(fig3_df$series, levels = names(SERIES_COL))
fig3_df$type   <- factor(fig3_df$type, levels = names(TYPE_LTY))

# Proportion of Japanese patients (reference line in the first panel)
fig3_ref <- data.frame(panel = factor(PANEL3_LAB[["f"]], levels = PANEL3_LAB),
                       value = SET$N_REGION[1L] / SET$N_TOTAL)

Figure3 <- ggplot(fig3_df,
                  aes(x = tau, y = value, colour = series, linetype = type,
                      group = interaction(series, type))) +
  geom_hline(data = fig3_ref, aes(yintercept = value), colour = "grey55",
             linetype = "dotted", linewidth = LW, inherit.aes = FALSE) +
  geom_line(linewidth = LW) +
  geom_point(size = PS, shape = 16) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_x_continuous(name = X_TAU, breaks = SET$TAU_GRID) +
  scale_y_continuous(name = NULL) +
  scale_colour_manual(name = NULL, values = SERIES_COL,
                      breaks = c("Analysis 1", "Analysis 2", "Analysis 3",
                                 "Overall")) +
  scale_linetype_manual(name = NULL, values = TYPE_LTY) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1)) +
  theme_paper()

print(Figure3)
save_figure(Figure3, "Figure3", width = 20, height = 9)

# =============================================================================
# Figure 4: required proportion of events in Japan at the final analysis
# =============================================================================
r4 <- fig4_data$required
fig4_df <- data.frame(method = r4$method, tau = r4$tau,
                      criterion = CRIT_LAB[r4$type], value = r4$fs_r1)
fig4_df$method    <- factor(METHOD_LAB[fig4_df$method], levels = METHOD_LAB)
fig4_df$criterion <- factor(fig4_df$criterion, levels = CRIT_LAB)

# Delays for which the target cannot be attained (crosses on the x-axis)
fig4_na <- fig4_df[is.na(fig4_df$value), ]
fig4_na$value <- 0

# Fixed design with the fixed-design number of events (dotted line)
fig4_fix <- data.frame(method = factor(METHOD_LAB[fig4_data$fixed$method],
                                       levels = METHOD_LAB),
                       value = fig4_data$fixed$fs_r1)

Figure4 <- ggplot(fig4_df[!is.na(fig4_df$value), ],
                  aes(x = tau, y = value, colour = criterion,
                      group = criterion)) +
  geom_hline(data = fig4_fix, aes(yintercept = value), colour = "grey55",
             linetype = "dotted", linewidth = LW, inherit.aes = FALSE) +
  geom_line(linewidth = LW) +
  geom_point(size = PS, shape = 16) +
  geom_point(data = fig4_na, shape = 4, size = 2 * PS, stroke = 2,
             show.legend = FALSE) +
  facet_wrap(~ method, nrow = 1, scales = "free_y") +
  expand_limits(y = 0) +
  scale_x_continuous(name = X_TAU, breaks = SET$TAU_GRID) +
  scale_y_continuous(name = "Required proportion of events in Japan") +
  scale_colour_manual(name = "Criterion", values = CRIT_COL) +
  theme_paper()

print(Figure4)
save_figure(Figure4, "Figure4", width = 16, height = 9)

# =============================================================================
# Figure 5: patient-level simulation versus closed form (OBF boundaries)
# =============================================================================
s5 <- fig5_data$stage
DELAY_LAB <- c("0" = "No delay",
               setNames(sprintf("Enrollment in Japan from month %g",
                                SET$APP_TAU), SET$APP_TAU))

fig5_sim <- s5[s5$source == "simulation", ]
fig5_cf  <- s5[s5$source == "formula_expected", ]
prep5 <- function(d) {
  data.frame(method = factor(METHOD_LAB[d$method], levels = METHOD_LAB),
             delay = factor(DELAY_LAB[as.character(d$tau)],
                            levels = DELAY_LAB),
             p1 = d$p1,
             series = factor(paste("Analysis", d$stage),
                             levels = names(SERIES_COL)),
             rcp = d$rcp,
             lower = d$rcp - qnorm(0.975) * d$se_rcp,
             upper = d$rcp + qnorm(0.975) * d$se_rcp)
}
fig5_sim <- prep5(fig5_sim)
fig5_cf  <- prep5(fig5_cf)

Figure5 <- ggplot(mapping = aes(x = p1, y = rcp, colour = series)) +
  geom_line(data = fig5_cf, aes(group = series), linewidth = LW) +
  geom_errorbar(data = fig5_sim, aes(ymin = lower, ymax = upper),
                width = 0.006, linewidth = LW) +
  geom_point(data = fig5_sim, size = PS, shape = 16) +
  facet_grid(method ~ delay) +
  scale_x_continuous(name = "Proportion of Japanese patients",
                     breaks = SET$P1_SIM) +
  scale_y_continuous(name = Y_RCP) +
  scale_colour_manual(name = NULL, values = SERIES_COL,
                      breaks = c("Analysis 1", "Analysis 2", "Analysis 3")) +
  guides(colour = guide_legend(nrow = 1)) +
  theme_paper()

print(Figure5)
save_figure(Figure5, "Figure5", width = 16, height = 13)

# =============================================================================
# Table 1: application (delayed start in Japan), by analysis
# =============================================================================
t1  <- table1_data
an1 <- t1$analysis
su1 <- t1$stage_unadj
sa1 <- t1$stage_adj
ou1 <- t1$overall_unadj
oa1 <- t1$overall_adj
K   <- SET$K

tab1_body <- vapply(seq_len(K), function(k) {
  cells <- c(sprintf("%d", k),
             fmt(an1$time[k], 1),
             sprintf("%d", an1$events[k]),
             fmt(an1$events_r1[k], 1),
             fmt(by_stage(su1, "M1", "upper")[k], 3),
             fmt(by_stage(su1, "M1", "p_eff")[k], 3))
  for (m in c("M1", "M2")) {
    cells <- c(cells,
               fmt(by_stage(su1, m, "rcp")[k], 3),
               fmt(by_stage(su1, m, "rcp_uncond")[k], 3),
               fmt(by_stage(sa1, m, "gamma")[k], 3),
               fmt(by_stage(sa1, m, "rcp")[k], 3))
  }
  paste(paste(cells, collapse = " & "), "\\\\")
}, character(1L))

tab1_overall <- paste(paste(c(
  "Overall", "", "", "", "",
  fmt(ou1$power[ou1$method == "M1"], 3),
  fmt(ou1$rcp[ou1$method == "M1"], 3), "", "",
  fmt(oa1$rcp[oa1$method == "M1"], 3),
  fmt(ou1$rcp[ou1$method == "M2"], 3), "", "",
  fmt(oa1$rcp[oa1$method == "M2"], 3)), collapse = " & "), "\\\\")

acc1 <- t1$accrual
tab1_caption <- paste0(
  "Application to a hypothetical multiregional clinical trial with a ",
  "delayed start of enrollment in Japan. ",
  "Group sequential design with $K = ", K, "$ analyses at information ",
  "fractions $1/3$, $2/3$ and $1$, O'Brien-Fleming-type efficacy ",
  "boundaries (Lan-DeMets spending function), one-sided ",
  "$\\alpha = ", SET$ALPHA, "$ and power $", sprintf("%.2f", 1 - SET$BETA),
  "$ for a hazard ratio of $", sprintf("%.2f", SET$HR), "$ in all regions. ",
  "A total of $", SET$N_TOTAL, "$ patients are enrolled uniformly over ",
  SET$T_ACCRUAL, " months (", SET$RATE, " patients per month), ",
  acc1$n_region[1L], " of them in Japan, where enrollment starts at month ",
  acc1$start[1L], "; the two other regions (",
  if (acc1$n_region[2L] == acc1$n_region[3L]) {
    paste(acc1$n_region[2L], "patients each")
  } else {
    paste(acc1$n_region[2L], "and", acc1$n_region[3L], "patients")
  },
  ") start at month 0 and absorb the delay, ",
  "so that the overall accrual is as planned. ",
  "The median survival time in the control group is ", SET$MEDIAN_C,
  " months. ",
  "$t_{k}$, expected calendar time of analysis $k$ (months); ",
  "$d_{k}$, number of events; ",
  "$d_{1,k}$, expected number of events in Japan; ",
  "$u_{k}$, efficacy boundary; ",
  "$\\Pr(E_{k})$, probability of stopping for efficacy at analysis $k$ ",
  "(overall row, power); ",
  "$\\phi_{k}$, regional consistency probability given efficacy stopping ",
  "at analysis $k$ with $\\gamma_{\\rm M1} = ", SET$GAMMA_M1,
  "$ and $\\gamma_{\\rm M2} = ", SET$GAMMA_M2,
  "$ (overall row, given success); ",
  "$\\psi_{k}$, unconditional probability of consistency at analysis $k$; ",
  "$\\gamma_{k}$, adjusted threshold; ",
  "$\\phi^{\\rm adj}_{k}$, regional consistency probability given efficacy ",
  "stopping with the adjusted threshold. ",
  "The reference null regional consistency probabilities of the fixed ",
  "design are ", fmt(ref1[["M1"]], 3), " (Method 1) and ",
  fmt(ref1[["M2"]], 3), " (Method 2).")

tab1_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  paste0("\\caption{", tab1_caption, "}"),
  "\\label{tab:application}",
  "\\begin{tabular}{rrrrrr rrrr rrrr}",
  "\\toprule",
  paste(" & & & & & &",
        "\\multicolumn{4}{c}{Method 1} &",
        "\\multicolumn{4}{c}{Method 2} \\\\"),
  "\\cmidrule(lr){7-10}\\cmidrule(lr){11-14}",
  paste("$k$ & $t_{k}$ & $d_{k}$ & $d_{1,k}$ & $u_{k}$ & $\\Pr(E_{k})$ &",
        "$\\phi_{k}$ & $\\psi_{k}$ & $\\gamma_{k}$ & $\\phi^{\\rm adj}_{k}$ &",
        "$\\phi_{k}$ & $\\psi_{k}$ & $\\gamma_{k}$ & $\\phi^{\\rm adj}_{k}$",
        "\\\\"),
  "\\midrule",
  tab1_body,
  "\\midrule",
  tab1_overall,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
write_tex(tab1_lines, file.path(results_dir, "table1_application.tex"))

# =============================================================================
# Table 2: required proportion and number of events in Japan (application)
# =============================================================================
r2   <- table2_data$required
fix2 <- table2_data$fixed
TAU2 <- sort(unique(r2$tau))
stopifnot(length(TAU2) == 2L)

req_cells <- function(type, m) {
  unlist(lapply(TAU2, function(tau) {
    row <- r2[r2$type == type & r2$method == m & r2$tau == tau, ]
    stopifnot(nrow(row) == 1L)
    c(fmt(row$fs_r1, 3), fmt(row$events_r1_ceiling, 0))
  }))
}
tab2_body <- vapply(names(CRIT_LAB), function(ty) {
  paste(paste(c(CRIT_LAB[[ty]], req_cells(ty, "M1"), req_cells(ty, "M2")),
              collapse = " & "), "\\\\")
}, character(1L))
fix_cells <- function(m) {
  row <- fix2[fix2$method == m, ]
  rep(c(fmt(row$fs_r1, 3), fmt(row$events_r1_ceiling, 0)), length(TAU2))
}
tab2_fixed <- paste(paste(c("Fixed design", fix_cells("M1"), fix_cells("M2")),
                          collapse = " & "), "\\\\")
any_na2 <- anyNA(r2$fs_r1)

tab2_caption <- paste0(
  "Required proportion $f_{1}$ and number $d_{1}$ of events in Japan at ",
  "the final analysis for a target value of ", SET$TARGET, " in the ",
  "application of Table~\\ref{tab:application}, without delay and with ",
  "enrollment in Japan starting at month ", TAU2[2L], ". ",
  "The proportions of events in Japan at the interim analyses are scaled ",
  "together with the final proportion, keeping the pattern implied by the ",
  "accrual. ",
  "Overall, regional consistency probability given success; ",
  "stagewise, smallest regional consistency probability given efficacy ",
  "stopping over the analyses; ",
  "unconditional, unconditional probability of consistency at the final ",
  "analysis. ",
  "The fixed design has $", fix2$events[1L], "$ events and uses the ",
  "unconditional criterion; it does not depend on the delay. ",
  "$d_{1}$ is rounded up to an integer.",
  if (any_na2) " ---, target not attainable." else "")

tau_head <- vapply(TAU2, function(tau) {
  if (tau == 0) "No delay" else sprintf("Month %g", tau)
}, character(1L))
tab2_lines <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\footnotesize",
  paste0("\\caption{", tab2_caption, "}"),
  "\\label{tab:required}",
  "\\begin{tabular}{l rrrr rrrr}",
  "\\toprule",
  "& \\multicolumn{4}{c}{Method 1} & \\multicolumn{4}{c}{Method 2} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  paste0("Start in Japan & ",
         paste(rep(sprintf("\\multicolumn{2}{c}{%s}", tau_head), 2L),
               collapse = " & "), " \\\\"),
  paste0("\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
         "\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}"),
  paste("Criterion &",
        paste(rep("$f_{1}$ & $d_{1}$", 2L * length(TAU2)),
              collapse = " & "), "\\\\"),
  "\\midrule",
  tab2_body,
  tab2_fixed,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
write_tex(tab2_lines, file.path(results_dir, "table2_required.tex"))

# =============================================================================
# Ledger of the numbers cited in the text: results/manuscript_numbers.csv
#   key          identifier used in number_locations.csv
#   value        unrounded value
#   digits       number of decimals in the text
#   text         value formatted with 'digits' (to be transcribed)
#   description  what the number is
#   source       data file and element
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
F1_REF <- 0.1   # proportion of events in Japan used for cited values

# ---- Design ----
add_num("d_fix", SET$D_FIX, 1,
        "Fixed-design number of events (Schoenfeld)", "settings$D_FIX")
for (b in names(SET$EVENTS)) {
  di <- SET$DESIGN_INFO[SET$DESIGN_INFO$boundary == b, ]
  for (k in seq_len(K)) {
    add_num(sprintf("events_%s_k%d", b, k), di$events[k], 0,
            sprintf("Number of events at analysis %d (%s)", k, b),
            "settings$DESIGN_INFO$events")
    add_num(sprintf("upper_%s_k%d", b, k), di$upper[k], 3,
            sprintf("Efficacy boundary at analysis %d (%s)", k, b),
            "settings$DESIGN_INFO$upper")
  }
}

# ---- Check of the overall design against rpact (OBF, no delay) ----
dc <- SET$DESIGN_CHECK
dc_val <- function(q, col) dc[[col]][dc$quantity == q]
add_num("check_power_rpact", dc_val("power", "rpact"), 3,
        "Power computed by rpact (OBF)", "settings$DESIGN_CHECK")
add_num("check_power_rcpgsd", dc_val("power", "rcpgsd"), 3,
        "Power computed by rcpgsd (OBF)", "settings$DESIGN_CHECK")
for (k in seq_len(K)) {
  q <- sprintf("time of analysis %d", k)
  add_num(sprintf("check_time_k%d_rpact", k), dc_val(q, "rpact"), 2,
          sprintf("Expected time of analysis %d by rpact (months)", k),
          "settings$DESIGN_CHECK")
  add_num(sprintf("check_time_k%d_absdiff", k), abs(dc_val(q, "difference")),
          2, sprintf("|rcpgsd - rpact| of the time of analysis %d (months)",
                     k), "settings$DESIGN_CHECK")
}

# ---- Figure 1 (f1 = F1_REF) ----
for (b in names(BOUNDARY_LAB)) {
  for (m in c("M1", "M2")) {
    sub <- s1[s1$boundary == b & s1$method == m &
                abs(s1$f1 - F1_REF) < 1e-9, ]
    sub <- sub[order(sub$stage), ]
    for (k in seq_len(K)) {
      add_num(sprintf("fig1_%s_%s_cond_k%d", b, m, k), sub$rcp[k], 3,
              sprintf("RCP given efficacy stop, analysis %d, %s, %s, f1 = %g",
                      k, b, m, F1_REF), "fig1_data$stage$rcp")
      add_num(sprintf("fig1_%s_%s_uncond_k%d", b, m, k), sub$rcp_uncond[k],
              3, sprintf("Unconditional RCP, analysis %d, %s, %s, f1 = %g",
                         k, b, m, F1_REF), "fig1_data$stage$rcp_uncond")
    }
    add_num(sprintf("fig1_%s_%s_overall", b, m),
            o1$rcp[o1$boundary == b & o1$method == m &
                     abs(o1$f1 - F1_REF) < 1e-9], 3,
            sprintf("Overall RCP given success, %s, %s, f1 = %g", b, m,
                    F1_REF), "fig1_data$overall$rcp")
  }
}
for (m in c("M1", "M2")) {
  sel <- x1$method == m & abs(x1$f1 - F1_REF) < 1e-9
  add_num(sprintf("fig1_fixed_%s_cond", m), x1$rcp[sel], 3,
          sprintf("Fixed design, RCP given success, %s, f1 = %g", m, F1_REF),
          "fig1_data$fixed$rcp")
  add_num(sprintf("fig1_fixed_%s_uncond", m), x1$rcp_uncond[sel], 3,
          sprintf("Fixed design, unconditional RCP, %s, f1 = %g", m, F1_REF),
          "fig1_data$fixed$rcp_uncond")
}
add_num("fig1_power_OBF", o1$power[o1$boundary == "OBF"][1L], 3,
        "Power, OBF", "fig1_data$overall$power")
add_num("fig1_power_Pocock", o1$power[o1$boundary == "Pocock"][1L], 3,
        "Power, Pocock", "fig1_data$overall$power")

# ---- Marginal RCP of Teng et al. (2018), Pr(C_k | Z_k >= u_k) ----
mg1 <- fig1_data$marginal
for (b in names(BOUNDARY_LAB)) {
  for (m in c("M1", "M2")) {
    for (k in seq_len(K)) {
      add_num(sprintf("fig1_%s_%s_marginal_k%d", b, m, k),
              mg1$rcp_marginal[mg1$boundary == b & mg1$method == m &
                                 mg1$stage == k &
                                 abs(mg1$f1 - F1_REF) < 1e-9], 3,
              sprintf(paste0("Marginal RCP of Teng et al. (2018), ",
                             "analysis %d, %s, %s, f1 = %g"), k, b, m,
                      F1_REF), "fig1_data$marginal$rcp_marginal")
    }
  }
}
mgt <- table1_data$marginal
for (m in c("M1", "M2")) {
  for (k in seq_len(K)) {
    add_num(sprintf("tab1_%s_marginal_k%d", m, k),
            mgt$rcp_marginal[mgt$method == m & mgt$stage == k], 3,
            sprintf(paste0("Marginal RCP of Teng et al. (2018) in the ",
                           "application, analysis %d, %s"), k, m),
            "table1_data$marginal$rcp_marginal")
  }
}

# ---- Figure 2 (f1 = F1_REF, and ranges over f1) ----
for (m in c("M1", "M2")) {
  add_num(sprintf("fig2_ref_%s", m),
          ref2$ref[ref2$method == m & ref2$boundary == "OBF" &
                     abs(ref2$f1 - F1_REF) < 1e-9], 3,
          sprintf("Fixed-design reference null RCP, %s, f1 = %g", m, F1_REF),
          "fig2_data$ref$ref")
  for (b in names(BOUNDARY_LAB)) {
    sub <- s2[s2$boundary == b & abs(s2$f1 - F1_REF) < 1e-9, ]
    sub <- sub[order(sub$stage), ]
    for (k in seq_len(K)) {
      add_num(sprintf("fig2_%s_%s_null_k%d", b, m, k),
              sub[[paste0("rcp0_", m, "_unadj")]][k], 3,
              sprintf("Null RCP (unadjusted), analysis %d, %s, %s, f1 = %g",
                      k, b, m, F1_REF),
              sprintf("fig2_data$stage$rcp0_%s_unadj", m))
      add_num(sprintf("fig2_%s_%s_gamma_k%d", b, m, k),
              sub[[paste0("gamma_", m, "_adj")]][k], 3,
              sprintf("Adjusted threshold, analysis %d, %s, %s, f1 = %g",
                      k, b, m, F1_REF),
              sprintf("fig2_data$stage$gamma_%s_adj", m))
    }
    for (k in seq_len(K)) {
      g <- s2[[paste0("gamma_", m, "_adj")]][s2$boundary == b &
                                                s2$stage == k]
      add_num(sprintf("fig2_%s_%s_gamma_k%d_range", b, m, k),
              max(g) - min(g), 3,
              sprintf("Range of the adjusted threshold over f1, analysis %d, %s, %s",
                      k, b, m), sprintf("fig2_data$stage$gamma_%s_adj", m))
    }
  }
}

# ---- Figure 3 (delay) ----
add_num("fig3_tau_limit", fig3_data$feasibility$tau_limit, 1,
        "Largest delay that the other regions can absorb (months)",
        "fig3_data$feasibility$tau_limit")
for (k in seq_len(K)) {
  add_num(sprintf("fig3_time_k%d", k), a3$time[a3$tau == 0 & a3$stage == k],
          1, sprintf("Expected time of analysis %d (months)", k),
          "fig3_data$analysis$time")
}
TAU_CITE <- c(SET$APP_TAU, max(SET$TAU_GRID))
for (tau in TAU_CITE) {
  for (k in seq_len(K)) {
    add_num(sprintf("fig3_tau%g_f1_k%d", tau, k),
            a3$f1[a3$tau == tau & a3$stage == k], 3,
            sprintf("Proportion of events in Japan, analysis %d, delay %g",
                    k, tau), "fig3_data$analysis$f1")
    add_num(sprintf("fig3_tau%g_events_r1_k%d", tau, k),
            a3$events_r1[a3$tau == tau & a3$stage == k], 1,
            sprintf("Expected events in Japan, analysis %d, delay %g", k,
                    tau), "fig3_data$analysis$events_r1")
    for (m in c("M1", "M2")) {
      add_num(sprintf("fig3_tau%g_%s_cond_k%d", tau, m, k),
              s3$rcp[s3$tau == tau & s3$method == m & s3$stage == k], 3,
              sprintf("RCP given efficacy stop, analysis %d, %s, delay %g",
                      k, m, tau), "fig3_data$stage$rcp")
    }
  }
}
for (tau in c(0, TAU_CITE)) {
  for (m in c("M1", "M2")) {
    add_num(sprintf("fig3_tau%g_%s_overall", tau, m),
            o3$rcp[o3$tau == tau & o3$method == m], 3,
            sprintf("Overall RCP given success, %s, delay %g", m, tau),
            "fig3_data$overall$rcp")
  }
}

# ---- Figure 4 (required proportion) ----
for (tau in c(0, max(SET$TAU_GRID))) {
  for (ty in c("overall", "unconditional")) {
    for (m in c("M1", "M2")) {
      row <- r4[r4$tau == tau & r4$type == ty & r4$method == m, ]
      add_num(sprintf("fig4_tau%g_%s_%s_fs", tau, ty, m), row$fs_r1, 3,
              sprintf("Required proportion of events in Japan, %s, %s, delay %g",
                      ty, m, tau), "fig4_data$required$fs_r1")
      add_num(sprintf("fig4_tau%g_%s_%s_events", tau, ty, m),
              row$events_r1_ceiling, 0,
              sprintf("Required events in Japan, %s, %s, delay %g", ty, m,
                      tau), "fig4_data$required$events_r1_ceiling")
    }
  }
}
for (m in c("M1", "M2")) {
  na_tau <- r4$tau[r4$type == "stagewise" & r4$method == m &
                     is.na(r4$fs_r1)]
  if (length(na_tau) > 0L) {
    add_num(sprintf("fig4_stagewise_%s_first_na", m), min(na_tau), 0,
            sprintf("Smallest delay with the stagewise target not attainable, %s",
                    m), "fig4_data$required$fs_r1")
  }
  add_num(sprintf("fig4_fixed_%s_fs", m),
          fig4_data$fixed$fs_r1[fig4_data$fixed$method == m], 3,
          sprintf("Required proportion of events in Japan, fixed design, %s",
                  m), "fig4_data$fixed$fs_r1")
}

# ---- Figure 5 (patient-level simulation) ----
an5 <- fig5_data$analysis
sim5 <- s5[s5$source == "simulation", ]
cf5  <- s5[s5$source == "formula_expected", ]
key5 <- function(d) paste(d$tau, d$p1, d$method, d$stage)
cf5  <- cf5[match(key5(sim5), key5(cf5)), ]
ev5  <- an5$events_r1_expected[match(paste(sim5$tau, sim5$p1, sim5$stage),
                                     paste(an5$tau, an5$p1, an5$stage))]
diff5 <- sim5$rcp - cf5$rcp
EV_MIN <- 5   # cells with at least this many expected events in Japan
for (m in c("M1", "M2")) {
  sel <- sim5$method == m & ev5 >= EV_MIN
  add_num(sprintf("fig5_%s_maxdiff_ev%d", m, EV_MIN), max(abs(diff5[sel])),
          3, sprintf(paste0("Largest |simulation - closed form| of the RCP ",
                            "given efficacy stop, %s, cells with at least %d ",
                            "expected events in Japan"), m, EV_MIN),
          "fig5_data$stage$rcp")
}
worst <- sim5$tau == SET$APP_TAU & sim5$p1 == min(SET$P1_SIM) &
  sim5$stage == 1L
for (m in c("M1", "M2")) {
  add_num(sprintf("fig5_worst_%s_sim", m), sim5$rcp[worst &
                                                      sim5$method == m], 3,
          sprintf("Simulated RCP, delay %g, p1 = %g, analysis 1, %s",
                  SET$APP_TAU, min(SET$P1_SIM), m), "fig5_data$stage$rcp")
  add_num(sprintf("fig5_worst_%s_formula", m), cf5$rcp[worst &
                                                         sim5$method == m],
          3, sprintf("Closed-form RCP, delay %g, p1 = %g, analysis 1, %s",
                     SET$APP_TAU, min(SET$P1_SIM), m), "fig5_data$stage$rcp")
}
w_an <- an5[an5$tau == SET$APP_TAU & an5$p1 == min(SET$P1_SIM) &
              an5$stage == 1L, ]
add_num("fig5_worst_events_r1", w_an$events_r1_expected, 1,
        "Expected events in Japan at analysis 1 in the worst scenario",
        "fig5_data$analysis$events_r1_expected")
add_num("fig5_worst_p_inf", w_an$p_inf_r1, 2,
        "Proportion of infinite Japanese estimates at analysis 1 (worst)",
        "fig5_data$analysis$p_inf_r1")
add_num("fig5_worst_p_no_est", w_an$p_no_est_r1, 3,
        "Proportion of trials without a Japanese estimate at analysis 1 (worst)",
        "fig5_data$analysis$p_no_est_r1")
add_num("fig5_nsim", SET$NSIM, 0, "Number of simulated trials per scenario",
        "settings$NSIM")
pw5 <- fig5_data$overall
pw5 <- pw5$power[pw5$source == "simulation" & pw5$method == "M1"]
add_num("fig5_power_min", min(pw5), 3, "Smallest simulated power",
        "fig5_data$overall$power")
add_num("fig5_power_max", max(pw5), 3, "Largest simulated power",
        "fig5_data$overall$power")

# ---- Table 1 (application) ----
for (m in c("M1", "M2")) {
  add_num(sprintf("tab1_ref_%s", m), ref1[[m]], 3,
          sprintf("Reference null RCP in the application, %s", m),
          "table1_data$ref$ref")
  add_num(sprintf("tab1_overall_unadj_%s", m), ou1$rcp[ou1$method == m], 3,
          sprintf("Overall RCP given success, unadjusted, %s", m),
          "table1_data$overall_unadj$rcp")
  add_num(sprintf("tab1_overall_adj_%s", m), oa1$rcp[oa1$method == m], 3,
          sprintf("Overall RCP given success, adjusted, %s", m),
          "table1_data$overall_adj$rcp")
}

# ---- Table 1 by analysis and Table 2 (values cited in the application) ----
for (k in seq_len(K)) {
  add_num(sprintf("tab1_peff_k%d", k), by_stage(su1, "M1", "p_eff")[k], 3,
          sprintf("Probability of stopping for efficacy at analysis %d (application)",
                  k), "table1_data$stage_unadj$p_eff")
  for (m in c("M1", "M2")) {
    add_num(sprintf("tab1_%s_cond_k%d", m, k), by_stage(su1, m, "rcp")[k], 3,
            sprintf("RCP given efficacy stop, analysis %d, %s (application)",
                    k, m), "table1_data$stage_unadj$rcp")
    add_num(sprintf("tab1_%s_uncond_k%d", m, k),
            by_stage(su1, m, "rcp_uncond")[k], 3,
            sprintf("Unconditional RCP, analysis %d, %s (application)", k, m),
            "table1_data$stage_unadj$rcp_uncond")
    add_num(sprintf("tab1_%s_gamma_k%d", m, k), by_stage(sa1, m, "gamma")[k],
            3, sprintf("Adjusted threshold, analysis %d, %s (application)", k,
                       m), "table1_data$stage_adj$gamma")
    add_num(sprintf("tab1_%s_adj_k%d", m, k), by_stage(sa1, m, "rcp")[k], 3,
            sprintf(paste0("RCP given efficacy stop with the adjusted ",
                           "threshold, analysis %d, %s (application)"), k, m),
            "table1_data$stage_adj$rcp")
  }
}
add_num("tab1_power", ou1$power[ou1$method == "M1"], 3,
        "Power (application)", "table1_data$overall_unadj$power")
for (tau in TAU2) {
  for (ty in names(CRIT_LAB)) {
    for (m in c("M1", "M2")) {
      row <- r2[r2$tau == tau & r2$type == ty & r2$method == m, ]
      if (is.na(row$fs_r1)) next
      add_num(sprintf("tab2_tau%g_%s_%s_fs", tau, ty, m), row$fs_r1, 3,
              sprintf("Required proportion of events in Japan, %s, %s, delay %g (application)",
                      ty, m, tau), "table2_data$required$fs_r1")
      add_num(sprintf("tab2_tau%g_%s_%s_events", tau, ty, m),
              row$events_r1_ceiling, 0,
              sprintf("Required events in Japan, %s, %s, delay %g (application)",
                      ty, m, tau), "table2_data$required$events_r1_ceiling")
    }
  }
}
for (m in c("M1", "M2")) {
  add_num(sprintf("fig4_fixed_%s_events", m),
          fix2$events_r1_ceiling[fix2$method == m], 0,
          sprintf("Required events in Japan, fixed design, %s", m),
          "table2_data$fixed$events_r1_ceiling")
}

ledger_df <- do.call(rbind, ledger)
rownames(ledger_df) <- NULL
utils::write.csv(ledger_df, file.path(results_dir, "manuscript_numbers.csv"),
                 row.names = FALSE)
message("Ledger saved: ", file.path(results_dir, "manuscript_numbers.csv"),
        " (", nrow(ledger_df), " numbers)")
