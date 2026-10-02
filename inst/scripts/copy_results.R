# =============================================================================
# copy_results.R
#
# Copies the figures and LaTeX tables of the manuscript and of the
# supplementary material to the manuscript folder and checks that each copy
# is identical to its source (MD5 checksum). Use this script instead of
# copying the files by hand.
#
# Usage
#   Open rcpgsd.Rproj (the working directory must be the package root) and
#   source this file after table_and_figure_manuscript.R and
#   table_and_figure_supplement.R:
#     source("inst/scripts/copy_results.R")
#   To copy to another folder (e.g., a new submission), assign DEST_DIR
#   before sourcing this file. A folder of a submitted version should not be
#   overwritten; create a new folder instead.
#
# Inputs (relative to the package root)
#   inst/scripts/results/   Figure1-5 (.pdf, .eps), table1, table2 (.tex)
#   inst/sup_info/results/  FigureS1-S5 (.pdf, .eps), tableS1 (.tex)
#
# Output
#   The files in DEST_DIR (all in one folder, so that the manuscript uses
#   \includegraphics{Figure1} and \input{table1_application}) and the log
#   DEST_DIR/copy_results.log with the MD5 checksums
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "rcpgsd")) {
  stop("Run this script from the root of the rcpgsd package ",
       "(open rcpgsd.Rproj first).", call. = FALSE)
}

if (!exists("DEST_DIR")) {
  DEST_DIR <- file.path("..", "..", "Manuscript")
}
if (!dir.exists(DEST_DIR)) {
  stop(sprintf("The destination folder '%s' does not exist.", DEST_DIR),
       call. = FALSE)
}

SOURCE_DIRS <- file.path("inst", c("scripts", "sup_info"), "results")
src <- unlist(lapply(SOURCE_DIRS, list.files,
                     pattern = "\\.(pdf|eps|tex)$", full.names = TRUE))
if (length(src) == 0L) {
  stop("No figures or tables found; run the rendering scripts first.",
       call. = FALSE)
}
if (anyDuplicated(basename(src))) {
  stop("Files with the same name exist in the result folders: ",
       paste(unique(basename(src)[duplicated(basename(src))]),
             collapse = ", "), call. = FALSE)
}

dest <- file.path(DEST_DIR, basename(src))
ok <- file.copy(src, dest, overwrite = TRUE, copy.date = TRUE)
if (!all(ok)) {
  stop("Copy failed: ", paste(basename(src)[!ok], collapse = ", "),
       call. = FALSE)
}

md5_src  <- unname(tools::md5sum(src))
md5_dest <- unname(tools::md5sum(dest))
res <- data.frame(file   = basename(src),
                  source = dirname(src),
                  md5    = md5_src,
                  status = ifelse(md5_src == md5_dest, "OK", "MISMATCH"),
                  stringsAsFactors = FALSE)

log_path <- file.path(DEST_DIR, "copy_results.log")
out <- c(sprintf("copy_results.R -- %s", format(Sys.time())),
         sprintf("Destination: %s", normalizePath(DEST_DIR)),
         "",
         utils::capture.output(print(res, row.names = FALSE, right = FALSE)))
writeLines(out, log_path)
writeLines(out)

if (any(res$status != "OK")) {
  stop("MD5 mismatch for: ",
       paste(res$file[res$status != "OK"], collapse = ", "), call. = FALSE)
}
message(sprintf("%d files copied to %s and verified by MD5.", nrow(res),
                DEST_DIR))
