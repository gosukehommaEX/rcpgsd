# Scripts for the figures and tables of the article

This folder and `inst/sup_info/` contain the code that produces every figure,
table and number of the article on regional consistency probabilities (RCPs)
under group sequential designs. Computation and rendering are separated, so
that figures and tables can be changed without repeating the computation.

All scripts are run from the package root (open `rcpgsd.Rproj`).

## Folder structure

```
inst/scripts/
  data_generate_main.R           computation  -> data/*.rds, data/csv/*.csv,
                                                 data/data_generate_main.log
  table_and_figure_manuscript.R  rendering    -> results/Figure1-5.pdf/.eps,
                                                 results/table1_application.tex,
                                                 results/table2_required.tex,
                                                 results/manuscript_numbers.csv
  number_check.R                 check of the numbers in the manuscript
                                              -> number_check_output.txt
  number_locations.csv           positions of the numbers in the manuscript
                                 (maintained by hand)
  copy_results.R                 copy of figures and tables to the manuscript
                                 folder with MD5 verification
inst/sup_info/
  data_generate_supplement.R     computation  -> data/siS1-siS4_data.rds, ...
  table_and_figure_supplement.R  rendering    -> results/FigureS1-S5.pdf/.eps,
                                                 results/tableS1_simulation.tex,
                                                 results/supplement_numbers.csv
```

## Order of execution

```r
devtools::load_all()
source("inst/scripts/data_generate_main.R")           # about 40 minutes
source("inst/scripts/table_and_figure_manuscript.R")  # less than 1 minute
source("inst/sup_info/data_generate_supplement.R")    # about 3 minutes
source("inst/sup_info/table_and_figure_supplement.R") # less than 1 minute
source("inst/scripts/copy_results.R")
source("inst/scripts/number_check.R")
```

The patient-level simulation (Figure 5) takes most of the time. Parts of the
computation can be skipped with the vector `RUN` at the top of each
`data_generate_*.R` script. The supplementary scripts read the settings and
the adjusted thresholds from the main data, so `data_generate_main.R` must be
run first.

## Figures and tables

| Output | Content | Data (`data/`) | Computation |
|---|---|---|---|
| Figure 1 | RCP under H1 by analysis (given efficacy stop and unconditional), overall RCP and fixed design, versus the proportion of events in Japan | `fig1_data.rds` | `rcpgsd()` |
| Figure 2 | Null RCP with the unadjusted thresholds and adjusted thresholds | `fig2_data.rds` | `rcpgsd_get_gamma()` |
| Figure 3 | Delayed start of enrolment in Japan: proportion of events in Japan and RCP | `fig3_data.rds` | `rcpgsd_accrual()`, `rcpgsd_events()`, `rcpgsd()` |
| Figure 4 | Required proportion of events in Japan by delay and criterion | `fig4_data.rds` | `rcpgsd_get_fs()` |
| Figure 5 | Patient-level simulation versus closed form | `fig5_data.rds` | `rcpgsd_sim()`, `rcpgsd()` |
| Table 1 | Application with a 6-month delay in Japan, unadjusted and adjusted thresholds | `table1_data.rds` | as Figures 2 and 3 |
| Table 2 | Required proportion and number of events in Japan in the application | `table2_data.rds` | subset of Figure 4 |
| Figure S1 | RCP under H1 with the adjusted and the unadjusted thresholds | `siS1_data.rds`, `fig1_data.rds` | `rcpgsd()` |
| Figure S2 | Binding and non-binding futility boundaries | `siS2_data.rds` | `rcpgsd()` |
| Figure S3 | Hazard reduction scale, exact and delta-method evaluation | `siS3_data.rds` | `rcpgsd()` |
| Figure S4 | Hazard ratio in Japan different from the other regions | `siS4_data.rds` | `rcpgsd()` |
| Figure S5 | RCP given efficacy stop versus the marginal RCP of Teng et al. (2018), Pr(C_k \| Z_k >= u_k) | `fig1_data.rds` | `rcpgsd()` |
| Table S1 | Details of the patient-level simulation | `fig5_data.rds` | as Figure 5 |

The figures are saved as PDF (fonts embedded) and EPS (vector, no
transparency). The tables are complete `table` environments with caption and
label; they use the booktabs rules (`\usepackage{booktabs}`) and are included
with `\input{table1_application}`. Every data frame of the `.rds` files is
also written to `data/csv/` for inspection.

## Settings

The overall design is computed with rpact: Lan-DeMets O'Brien-Fleming-type
or Pocock-type efficacy boundaries (`getDesignGroupSequential()`) and the
numbers of events for a power of 0.90 (`getSampleSizeSurvival()`, Schoenfeld
approximation, consistent with the information used for the RCPs). The final
number of events is rounded up and the interim numbers are rounded to the
planned information fractions. `data_generate_main.R` checks that the
expected analysis times and the power computed by rcpgsd agree with rpact
and stops otherwise; the comparison is written to the log and to
`settings$DESIGN_CHECK` of each `.rds` file. All settings are stored in
`settings` of each `.rds` file, and the package version, git commit and
session information in `meta` and in the log files.

## Checking the numbers in the text

1. `table_and_figure_manuscript.R` writes the ledger
   `results/manuscript_numbers.csv` (key, unrounded value, number of
   decimals, text to be transcribed, description, source).
2. Transcribe the `text` column into the manuscript.
3. Add one row per number to `number_locations.csv`: the ledger key and the
   words around the number, with `{}` in place of the number, e.g.
   `fig1_OBF_M1_cond_k1,the first interim analysis was {} for Method 1,,`.
   Use `scale = 100` for percentages and key `-` for numbers that are
   deliberately not in the ledger.
4. Run `number_check.R`. It reports, for each location, the number in the
   text, the recomputed value and the expected text, lists decimal numbers of
   the text that are not registered, and stops with an error if any check
   fails. The manuscript files are `../../Manuscript/main.tex` and
   `../../Manuscript/supporting_information.tex` by default
   (`MANUSCRIPT_FILES` can be assigned before sourcing).

The ledger of the supplementary material is
`inst/sup_info/results/supplement_numbers.csv`.
