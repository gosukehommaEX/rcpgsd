# =============================================================================
# number_check.R
#
# Checks the numbers transcribed into the manuscript against the values
# computed in R. The reference values are the ledger
# inst/scripts/results/manuscript_numbers.csv written by
# table_and_figure_manuscript.R from the .rds files of data_generate_main.R.
# The positions of the numbers in the manuscript are given by
# inst/scripts/number_locations.csv, so the .tex files need no markers.
#
# Usage
#   Open rcpgsd.Rproj (the working directory must be the package root) and
#   source this file:
#     source("inst/scripts/number_check.R")
#   The manuscript files are MANUSCRIPT_FILES below. To check other files
#   (e.g., the supplementary material or a response letter), assign
#   MANUSCRIPT_FILES before sourcing this file.
#
# number_locations.csv (one row per number in the manuscript)
#   key      ledger key, or "-" for a number that is deliberately not in the
#            ledger (e.g., a value quoted from another article); such rows
#            are only located, not compared
#   context  the words around the number as they appear in the .tex file,
#            with {} in place of the number, e.g.
#              "efficacy stopping at the first analysis was {} for Method 1"
#            Runs of spaces match any white space. The context must lie on
#            one line of the .tex file (one paragraph per line) and must
#            match exactly once in MANUSCRIPT_FILES.
#   scale    optional factor applied to the ledger value before the
#            comparison, e.g. 100 when 0.36 is written as 36\%; empty = 1
#   note     optional free text (not used)
#
# Checks and outcome
#   1. The ledger is newer than the .rds files and the rendering script
#      (otherwise table_and_figure_manuscript.R must be re-run).
#   2. Each context is found exactly once.
#   3. The number in the text equals the ledger value rounded to the number
#      of decimals of the ledger (both as a string and numerically, within
#      half a unit of the last decimal).
#   4. Decimal numbers and percentages in the text that are not covered by
#      number_locations.csv are listed as unregistered numbers (for review;
#      not a failure). Comments, display mathematics, \cite, \ref, \label,
#      \input, \includegraphics and lengths are not scanned.
#   The result is printed and written to inst/scripts/number_check_output.txt
#   with the session information. The script stops with an error at the end
#   if any check fails.
# =============================================================================

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "rcpgsd")) {
  stop("Run this script from the root of the rcpgsd package ",
       "(open rcpgsd.Rproj first).", call. = FALSE)
}

if (!exists("MANUSCRIPT_FILES")) {
  MANUSCRIPT_FILES <- file.path("..", "..", "Manuscript", "main.tex")
}

script_dir  <- file.path("inst", "scripts")
data_dir    <- file.path(script_dir, "data")
ledger_path <- file.path(script_dir, "results", "manuscript_numbers.csv")
loc_path    <- file.path(script_dir, "number_locations.csv")
out_path    <- file.path(script_dir, "number_check_output.txt")
render_path <- file.path(script_dir, "table_and_figure_manuscript.R")
if (file.exists(out_path)) invisible(file.remove(out_path))

emit <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = out_path, append = TRUE)
}
emit_table <- function(df) {
  if (nrow(df) == 0L) {
    emit("  (none)")
    return(invisible(NULL))
  }
  out <- utils::capture.output(print(df, row.names = FALSE, right = FALSE))
  for (l in out) emit("  ", l)
}

n_fail <- 0L

emit("==================================================================")
emit(" Number check: manuscript versus ledger")
emit(" ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
emit("==================================================================")

# =============================================================================
# 1. Ledger up to date
# =============================================================================
if (!file.exists(ledger_path)) {
  stop(sprintf("'%s' not found; run table_and_figure_manuscript.R first.",
               ledger_path), call. = FALSE)
}
ledger <- utils::read.csv(ledger_path, stringsAsFactors = FALSE,
                          colClasses = c(key = "character",
                                         value = "numeric",
                                         digits = "integer",
                                         text = "character",
                                         description = "character",
                                         source = "character"))
t_ledger <- file.mtime(ledger_path)
inputs   <- c(list.files(data_dir, pattern = "\\.rds$", full.names = TRUE),
              render_path)
newer    <- inputs[file.mtime(inputs) > t_ledger]
emit("")
emit("[1] Ledger: ", ledger_path, " (", nrow(ledger), " numbers, ",
     format(t_ledger, "%Y-%m-%d %H:%M:%S"), ")")
if (length(newer) > 0L) {
  n_fail <- n_fail + 1L
  emit("  FAIL  files newer than the ledger (re-run ",
       "table_and_figure_manuscript.R):")
  for (f in newer) emit("        ", f)
} else {
  emit("  PASS  the ledger is newer than the .rds files and the ",
       "rendering script")
}

# =============================================================================
# 2. Manuscript files and locations
# =============================================================================
locs <- utils::read.csv(loc_path, stringsAsFactors = FALSE,
                        colClasses = "character")
if (!all(c("key", "context") %in% names(locs))) {
  stop("number_locations.csv must have the columns 'key' and 'context'.",
       call. = FALSE)
}
if (!"scale" %in% names(locs)) locs$scale <- ""
locs$scale[is.na(locs$scale) | trimws(locs$scale) == ""] <- "1"

emit("")
emit("[2] Manuscript files and locations (", nrow(locs), " rows in ",
     loc_path, ")")
tex_ok <- file.exists(MANUSCRIPT_FILES)
for (i in seq_along(MANUSCRIPT_FILES)) {
  emit("  ", if (tex_ok[i]) "found     " else "NOT FOUND ",
       MANUSCRIPT_FILES[i])
}
if (!all(tex_ok)) {
  if (nrow(locs) > 0L) {
    n_fail <- n_fail + 1L
    emit("  FAIL  manuscript file(s) not found")
  } else {
    emit("  INFO  no manuscript and no locations yet; only the ledger was ",
         "checked")
  }
}

bad_key <- locs$key != "-" & !(locs$key %in% ledger$key)
if (any(bad_key)) {
  n_fail <- n_fail + 1L
  emit("  FAIL  keys not in the ledger: ",
       paste(unique(locs$key[bad_key]), collapse = ", "))
}
n_ph <- vapply(gregexpr("{}", locs$context, fixed = TRUE),
               function(m) sum(m > 0), integer(1L))
if (any(n_ph != 1L)) {
  n_fail <- n_fail + 1L
  emit("  FAIL  contexts without exactly one {} placeholder (rows ",
       paste(which(n_ph != 1L), collapse = ", "), ")")
}

# Read the .tex files: comments removed, and a masked copy in which the
# regions that are not scanned for unregistered numbers are blanked out
DISPLAY_ENVS <- "(equation|eqnarray|align|gather|multline|displaymath)\\*?"
DISPLAY_RE   <- paste0("\\\\begin\\{", DISPLAY_ENVS, "\\}|\\\\end\\{",
                       DISPLAY_ENVS, "\\}|(?<!\\\\)\\\\\\[|(?<!\\\\)\\\\\\]")
SKIP_RE <- paste0(
  "\\\\(cite[a-z]*|ref|eqref|autoref|label|input|include|includegraphics|",
  "url|href|bibliography|bibliographystyle|hspace|vspace|setlength|",
  "addtolength|documentclass|usepackage)\\*?(\\[[^]]*\\])*\\{[^}]*\\}",
  "|[0-9]*\\.?[0-9]+(pt|cm|mm|in|em|ex|bp)\\b",
  "|[0-9]*\\.?[0-9]+\\s*\\\\(textwidth|linewidth|columnwidth|textheight)")

blank <- function(chars, from, len) {
  if (len > 0L) chars[from:(from + len - 1L)] <- " "
  chars
}

read_tex <- function(path) {
  raw   <- readLines(path, encoding = "UTF-8", warn = FALSE)
  text  <- gsub("(?<!\\\\)%.*$", "", raw, perl = TRUE)
  mask  <- text
  in_math <- FALSE
  for (i in seq_along(text)) {
    chars <- strsplit(text[i], "")[[1L]]
    if (length(chars) == 0L) next
    m <- gregexpr(DISPLAY_RE, text[i], perl = TRUE)[[1L]]
    pos <- 1L
    if (m[1L] > 0L) {
      for (j in seq_along(m)) {
        tok <- substr(text[i], m[j], m[j] + attr(m, "match.length")[j] - 1L)
        opening <- grepl("^\\\\begin|^\\\\\\[", tok)
        if (opening && !in_math) {
          in_math <- TRUE
          pos <- m[j]
        } else if (!opening && in_math) {
          chars <- blank(chars, pos, m[j] + attr(m, "match.length")[j] - pos)
          in_math <- FALSE
        }
      }
    }
    if (in_math) chars <- blank(chars, pos, length(chars) - pos + 1L)
    s <- gregexpr(SKIP_RE, text[i], perl = TRUE)[[1L]]
    if (s[1L] > 0L) {
      for (j in seq_along(s)) {
        chars <- blank(chars, s[j], attr(s, "match.length")[j])
      }
    }
    mask[i] <- paste(chars, collapse = "")
  }
  list(path = path, text = text, mask = mask)
}
tex <- lapply(MANUSCRIPT_FILES[tex_ok], read_tex)

# =============================================================================
# 3. Compare the numbers at the locations with the ledger
# =============================================================================
NUM_RE <- "(-?[0-9]+(?:\\.[0-9]+)?)"
escape_re <- function(x) {
  x <- gsub("([][{}()|.^$*+?\\\\])", "\\\\\\1", x, perl = TRUE)
  gsub("\\s+", "\\\\s+", x, perl = TRUE)
}
context_re <- function(context) {
  p <- regexpr("{}", context, fixed = TRUE)
  before <- substr(context, 1L, p - 1L)
  after  <- substr(context, p + 2L, nchar(context))
  paste0(escape_re(before), "(?<![0-9.])", NUM_RE, "(?![0-9])",
         escape_re(after))
}

covered <- list()   # spans of located numbers (file, line, start, end)
res <- vector("list", nrow(locs))
for (i in seq_len(nrow(locs))) {
  key <- locs$key[i]
  row <- data.frame(key = key, location = "", in_text = "",
                    recomputed = "", expected = "", status = "FAIL",
                    note = "", stringsAsFactors = FALSE)
  if (n_ph[i] != 1L || length(tex) == 0L) {
    row$note <- if (n_ph[i] != 1L) "no single {} in context" else
      "no manuscript"
    res[[i]] <- row
    next
  }
  re   <- context_re(locs$context[i])
  hits <- list()
  for (f in seq_along(tex)) {
    for (l in seq_along(tex[[f]]$text)) {
      m <- gregexpr(re, tex[[f]]$text[l], perl = TRUE)[[1L]]
      if (m[1L] < 0L) next
      cs <- attr(m, "capture.start")[, 1L]
      cl <- attr(m, "capture.length")[, 1L]
      for (j in seq_along(cs)) {
        hits[[length(hits) + 1L]] <- list(
          file = f, line = l, start = cs[j], end = cs[j] + cl[j] - 1L,
          num = substr(tex[[f]]$text[l], cs[j], cs[j] + cl[j] - 1L))
      }
    }
  }
  if (length(hits) != 1L) {
    row$note <- if (length(hits) == 0L) "context not found" else
      sprintf("context found %d times", length(hits))
    if (length(hits) > 1L) {
      row$location <- paste(vapply(hits, function(h) {
        sprintf("%s:%d", basename(tex[[h$file]]$path), h$line)
      }, character(1L)), collapse = ", ")
    }
    res[[i]] <- row
    next
  }
  h <- hits[[1L]]
  covered[[length(covered) + 1L]] <- h
  row$location <- sprintf("%s:%d", basename(tex[[h$file]]$path), h$line)
  row$in_text  <- h$num
  if (key == "-") {
    row$status <- "SKIP"
    row$note   <- "not in the ledger by design"
    res[[i]] <- row
    next
  }
  if (!key %in% ledger$key) {
    row$note <- "key not in the ledger"
    res[[i]] <- row
    next
  }
  lg     <- ledger[ledger$key == key, ]
  scale  <- as.numeric(locs$scale[i])
  d_text <- max(0L, lg$digits - as.integer(round(log10(scale))))
  target <- lg$value * scale
  row$recomputed <- format(signif(target, 8))
  row$expected   <- sprintf("%.*f", d_text, target)
  ok_text <- identical(h$num, row$expected)
  ok_num  <- abs(as.numeric(h$num) - target) <= 0.5 * 10^(-d_text) + 1e-9
  if (ok_text && ok_num) {
    row$status <- "PASS"
  } else {
    row$note <- if (ok_num) "rounding or number of decimals differs" else
      "value differs"
  }
  res[[i]] <- row
}
res <- if (length(res) > 0L) do.call(rbind, res) else
  data.frame(key = character(0), location = character(0),
             in_text = character(0), recomputed = character(0),
             expected = character(0), status = character(0),
             note = character(0))
n_fail <- n_fail + sum(res$status == "FAIL")

emit("")
emit("[3] Numbers at the registered locations")
emit_table(res)
emit("  PASS ", sum(res$status == "PASS"), ", FAIL ",
     sum(res$status == "FAIL"), ", SKIP ", sum(res$status == "SKIP"))

# =============================================================================
# 4. Unregistered numbers in the text (for review)
# =============================================================================
SCAN_RE <- "(?<![0-9.A-Za-z_])([0-9]*\\.[0-9]+|[0-9]+(?=\\\\%))"
unreg <- list()
for (f in seq_along(tex)) {
  for (l in seq_along(tex[[f]]$mask)) {
    m <- gregexpr(SCAN_RE, tex[[f]]$mask[l], perl = TRUE)[[1L]]
    if (m[1L] < 0L) next
    for (j in seq_along(m)) {
      st <- m[j]
      en <- m[j] + attr(m, "match.length")[j] - 1L
      is_cov <- any(vapply(covered, function(h) {
        h$file == f && h$line == l && st <= h$end && en >= h$start
      }, logical(1L)))
      if (is_cov) next
      line <- tex[[f]]$text[l]
      unreg[[length(unreg) + 1L]] <- data.frame(
        location = sprintf("%s:%d", basename(tex[[f]]$path), l),
        number   = substr(line, st, en),
        snippet  = substr(line, max(1L, st - 35L),
                          min(nchar(line), en + 35L)),
        stringsAsFactors = FALSE)
    }
  }
}
unreg <- if (length(unreg) > 0L) do.call(rbind, unreg) else
  data.frame(location = character(0), number = character(0),
             snippet = character(0))

emit("")
emit("[4] Unregistered numbers (decimals and percentages not covered by ",
     "number_locations.csv)")
emit_table(unreg)

unused <- setdiff(ledger$key, locs$key)
emit("")
emit("[5] Ledger keys not cited in the manuscript: ", length(unused), " of ",
     nrow(ledger), " (for information)")

# =============================================================================
# Summary and session information
# =============================================================================
emit("")
emit("==================================================================")
emit(" Result: ", if (n_fail == 0L) "ALL CHECKS PASSED" else
  sprintf("%d FAILURE(S)", n_fail), "; ", nrow(unreg),
  " unregistered number(s) to review")
emit("==================================================================")
emit("")
for (l in utils::capture.output(sessionInfo())) emit(l)

if (n_fail > 0L) {
  stop(sprintf("number_check: %d failure(s); see %s.", n_fail, out_path),
       call. = FALSE)
}
