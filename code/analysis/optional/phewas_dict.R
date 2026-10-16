# Code dictionary of the NHANES III adult, exam and lab files, used by phewas_clean.R to recognise missing codes.
#   Fill codes (88, 888, ...) depend on field width and can be valid values (e.g. 88 for ALP or MCV), so each variable's
#   codes, labels and frequencies are read from the codebook PDFs and missing codes are identified by label.
#   Codebook frequencies are checked against counts in DuckDB; variables that do not match are excluded.
# Input: data/raw/nhanes3/{adult,exam,lab}.sas, data/raw/nhanes3/docs/{adult,exam,lab}-acc.pdf (not shipped; needs pdftotext)
# Output: results/analysis/phewas/ewas_codes.csv (one row per variable and code), ewas_vars.csv, ewas_dict_check.csv
# Run from the repository root: Rscript code/analysis/optional/phewas_dict.R
suppressMessages({library(DBI); library(duckdb)})
raw <- here::here("data", "raw", "nhanes3"); out <- here::here("results", "analysis", "phewas")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
FILES <- c(adult = "nhanes3_adult", exam = "nhanes3_exam", lab = "nhanes3_lab")
NROW <- c(adult = 20050, exam = 31311, lab = 29314)

# ---------- 1. SAS layout: name, start/end position, width, character type, label ----------
read_sas <- function(f) {
  L <- readLines(file.path(raw, paste0(f, ".sas")), warn = FALSE)
  i0 <- grep("^\\s*INPUT\\s*$", L); i1 <- i0 + which(grepl("^\\s*;\\s*$", L[(i0 + 1):length(L)]))[1]
  m <- regmatches(L[i0:i1], regexec("^\\s*([A-Z][A-Z0-9_]*)\\s+(\\$\\s*)?([0-9]+)(-([0-9]+))?\\s*$", L[i0:i1]))
  m <- do.call(rbind, lapply(m[lengths(m) > 0], function(x) data.frame(variable = x[2], char = nzchar(x[3]), start = as.integer(x[4]),
    end = ifelse(nzchar(x[6]), as.integer(x[6]), as.integer(x[4])))))
  j0 <- grep("^\\s*LABEL\\s*$", L); lb <- regmatches(L[j0:length(L)], regexec('^\\s*([A-Z][A-Z0-9_]*)\\s*=\\s*"(.*)"\\s*$', L[j0:length(L)]))
  lb <- do.call(rbind, lapply(lb[lengths(lb) > 0], function(x) data.frame(variable = x[2], sas_label = x[3])))
  m <- merge(m, lb, by = "variable", all.x = TRUE); m$width <- m$end - m$start + 1; m$file <- f
  m[order(m$start), ]
}

# ---------- 2. PDF to text ----------
pdf_txt <- function(f) {
  tf <- tempfile(fileext = ".txt")
  stopifnot(system2("pdftotext", c("-layout", shQuote(file.path(raw, "docs", paste0(f, "-acc.pdf"))), shQuote(tf))) == 0)
  gsub("\f", "", readLines(tf, warn = FALSE))
}

# ---------- 3. Index: variable -> section / subsection ----------
# index lines look like "Doctor ever told you had: arthritis .........   HAC1A              1465"; section headings are
#   all caps without leader dots, in column 0 (sections) or indented (subsections)
read_index <- function(L, sas) {
  hdr <- grepl("^\\s*[A-Z][A-Z0-9 ,/()&:;.'+-]{3,}$", L) & !grepl("\\.{3,}|NHANES III|^\\s*(NOTES|CODE|TOPIC)\\b|\\.\\s*$|[a-z]", L)
  cur_top <- NA; cur_sub <- NA; res <- list()
  vre <- "^(.*?)\\s*\\.{2,}\\s*([A-Z][A-Z0-9]+)\\s+([0-9]+)(-([0-9]+))?\\s*$"
  for (k in seq_along(L)) {
    s <- L[k]
    if (hdr[k]) { if (grepl("^[A-Z]", s)) { cur_top <- trimws(s); cur_sub <- NA } else cur_sub <- trimws(s); next }
    m <- regmatches(s, regexec(vre, s))[[1]]
    if (length(m) && m[3] %in% sas$variable) {
      st <- as.integer(m[4]); v <- m[3]
      if (st == sas$start[sas$variable == v]) res[[v]] <- data.frame(variable = v, section = cur_top, subsection = cur_sub, index_label = trimws(m[2]))
    }
  }
  do.call(rbind, res)
}

# ---------- 4. Body: code lines of each variable ----------
# variable block: position line ("    1423-1426    Lead (ug/dL)"), name line ("  PBP   2342   00.7 Below level of detection"),
#   then code lines ("   24476   0001-71.8"). A name line counts only after a position line matching the SAS layout.
read_body <- function(L, sas) {
  pos <- setNames(sas$start, sas$variable); en <- setNames(sas$end, sas$variable)
  nre <- "^\\s{0,8}([A-Z][A-Z0-9]{1,9})(\\s+(.*))?$"; cre <- "^\\s+([0-9]+)\\s+(\\S.*?)\\s*$"
  prev_nonempty <- function(k) { j <- k - 1; while (j > 0 && !nzchar(trimws(L[j]))) j <- j - 1; j }
  is_pos <- function(s, v) { m <- regmatches(s, regexec("^\\s*([0-9]+)(-([0-9]+))?(\\s{2,}.*)?$", s))[[1]]
    length(m) > 0 && as.integer(m[2]) == pos[[v]] && (if (nzchar(m[4])) as.integer(m[4]) == en[[v]] else pos[[v]] == en[[v]]) }
  starts <- integer(); vars <- character()
  for (k in seq_along(L)) {
    m <- regmatches(L[k], regexec(nre, L[k]))[[1]]
    if (!length(m) || !(m[2] %in% names(pos))) next
    j <- prev_nonempty(k)
    if (j > 0 && is_pos(L[j], m[2])) { starts <- c(starts, k); vars <- c(vars, m[2]) }
  }
  out <- list()
  for (i in seq_along(starts)) {
    k <- starts[i]; v <- vars[i]; stop_at <- if (i < length(starts)) starts[i + 1] - 1 else length(L)
    # the next variable's position line is not part of this block
    if (i < length(starts)) stop_at <- prev_nonempty(starts[i + 1]) - 1
    m <- regmatches(L[k], regexec(nre, L[k]))[[1]]; rest <- trimws(m[4])
    lines <- if (grepl("^[0-9]+\\s+\\S", rest)) paste("  ", rest) else character()
    blk <- L[(k + 1):max(k + 1, stop_at)]
    nt <- which(trimws(blk) == "NOTES"); if (length(nt)) blk <- blk[seq_len(nt[1] - 1)]
    lines <- c(lines, blk[grepl(cre, blk)])
    if (!length(lines)) { out[[length(out) + 1]] <- data.frame(variable = v, count = NA, code_raw = NA, line = k); next }
    cm <- regmatches(lines, regexec(cre, lines))
    out[[length(out) + 1]] <- data.frame(variable = v, count = as.numeric(sapply(cm, `[`, 2)), code_raw = sapply(cm, `[`, 3), line = k)
  }
  do.call(rbind, out)
}

# split code lines into: range / single value + label / Blank
SPECIAL_RE <- "blank but applicable|don'?t know|refused|not ascertained|cannot be determined|could not be determined|varies|varied|fewer than|less than one|missing"
parse_code <- function(x) {
  x <- trimws(x); blank <- grepl("^Blank\\s*$", x, ignore.case = TRUE)
  rng <- regmatches(x, regexec("^(-?[0-9]*\\.?[0-9]+)\\s*-\\s*(-?[0-9]*\\.?[0-9]+)(\\s+(.*))?$", x))
  one <- regmatches(x, regexec("^(-?[0-9]*\\.?[0-9]+)(\\s+(.*))?$", x))
  data.frame(
    is_blank = blank,
    is_range = lengths(rng) > 0,
    lo = ifelse(lengths(rng) > 0, suppressWarnings(as.numeric(sapply(rng, `[`, 2))), suppressWarnings(as.numeric(sapply(one, `[`, 2)))),
    hi = ifelse(lengths(rng) > 0, suppressWarnings(as.numeric(sapply(rng, `[`, 3))), suppressWarnings(as.numeric(sapply(one, `[`, 2)))),
    code_label = ifelse(lengths(rng) > 0, trimws(sapply(rng, `[`, 5)), ifelse(lengths(one) > 0, trimws(sapply(one, `[`, 4)), ifelse(blank, "Blank", x))),
    parsed = blank | lengths(rng) > 0 | lengths(one) > 0)
}

# ---------- 5. Per-file summary ----------
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
CODES <- NULL; VARS <- NULL; CHECK <- NULL
for (f in names(FILES)) {
  sas <- read_sas(f); L <- pdf_txt(f)
  idx <- read_index(L, sas); body <- read_body(L, sas)
  pc <- parse_code(body$code_raw); body <- cbind(body, pc)
  body$code_label[is.na(body$code_label)] <- ""
  body$is_special <- grepl(SPECIAL_RE, body$code_label, ignore.case = TRUE)
  body$file <- f
  # if a variable appears more than once in the body (index mismatch), keep the block with the most code lines
  nb <- aggregate(count ~ variable + line, body, length); nb <- nb[order(nb$variable, -nb$count), ]; keep <- nb[!duplicated(nb$variable), c("variable", "line")]
  body <- merge(body, keep, by = c("variable", "line"))
  cat(sprintf("\n== %s: SAS %d variables; index matched %d; code blocks in body %d; code lines %d (unparsed %d)\n", f, nrow(sas), nrow(idx), length(unique(body$variable)),
              nrow(body), sum(!body$parsed, na.rm = TRUE)))
  # unparsed lines are label continuations (e.g. "1 million population or more"), not codes -> dropped; character variables are not in the EWAS
  if (any(!body$parsed, na.rm = TRUE)) print(head(body[!body$parsed & !is.na(body$parsed), c("variable", "count", "code_raw")], 8))
  body <- body[is.na(body$parsed) | body$parsed, ]

  # --- actual DuckDB columns and types
  ty <- dbGetQuery(con, sprintf("select column_name as vname, column_type from (describe %s)", FILES[[f]])); names(ty)[1] <- "variable"
  v <- merge(sas, idx, by = "variable", all.x = TRUE); v <- merge(v, ty, by = "variable", all.x = TRUE)

  # --- check single-code frequencies, Blank vs NULL counts and total frequency against DuckDB
  num <- v$variable[v$column_type %in% c("DOUBLE", "INTEGER", "TINYINT", "BIGINT", "SMALLINT") & v$variable %in% body$variable]
  chk <- lapply(num, function(x) {
    b <- body[body$variable == x, ]
    q <- dbGetQuery(con, sprintf('select "%s" val, count(*) n from %s group by 1', x, FILES[[f]]))
    n_null <- sum(q$n[is.na(q$val)]); qv <- q[!is.na(q$val), ]
    single <- b[!b$is_range & !b$is_blank & b$parsed & !is.na(b$count), ]
    ok_single <- all(vapply(seq_len(nrow(single)), function(r) sum(qv$n[abs(qv$val - single$lo[r]) < 1e-9]) == single$count[r], TRUE))
    rng <- b[b$is_range & !is.na(b$count), ]
    in_rng <- vapply(seq_len(nrow(rng)), function(r) sum(qv$n[qv$val >= rng$lo[r] - 1e-9 & qv$val <= rng$hi[r] + 1e-9]), 0)
    # range counts may overlap single-value codes inside them (e.g. 0 = never) -> allow >=
    ok_rng <- all(in_rng >= rng$count)
    covered <- vapply(qv$val, function(z) any(abs(z - single$lo) < 1e-9) || any(z >= rng$lo - 1e-9 & z <= rng$hi + 1e-9), TRUE)
    blank_n <- sum(b$count[b$is_blank], na.rm = TRUE)
    data.frame(file = f, variable = x, n_codes = nrow(b), sum_counts = sum(b$count, na.rm = TRUE), n_rows = NROW[[f]],
      ok_total = sum(b$count, na.rm = TRUE) == NROW[[f]], ok_single = ok_single, ok_range = ok_rng,
      ok_blank = blank_n == n_null, n_uncovered_values = sum(qv$n[!covered]),
      uncovered_examples = paste(head(qv$val[!covered], 5), collapse = " "))
  })
  chk <- do.call(rbind, chk)
  chk$ok_all <- with(chk, ok_total & ok_single & ok_range & ok_blank & n_uncovered_values == 0)
  cat(sprintf("   numeric columns %d: total frequency = rows %d; all single-value frequencies match %d; ranges %d; Blank = NULL %d; no uncovered values %d; all passed %d\n", nrow(chk),
              sum(chk$ok_total), sum(chk$ok_single), sum(chk$ok_range), sum(chk$ok_blank), sum(chk$n_uncovered_values == 0), sum(chk$ok_all)))
  CODES <- rbind(CODES, body[, c("file", "variable", "count", "code_raw", "is_blank", "is_range", "lo", "hi", "code_label", "is_special", "parsed")])
  VARS <- rbind(VARS, v); CHECK <- rbind(CHECK, chk)
}
cat("\n== Labels classified as missing codes (all files) ==\n"); print(sort(table(tolower(CODES$code_label[CODES$is_special])), decreasing = TRUE))
cat("\n== Labels of repeated-digit single-value codes not classified as missing (for manual check, first 60) ==\n")
tok <- sub("\\s.*$", "", trimws(CODES$code_raw))
rep_digit <- with(CODES, !is_range & !is_blank & parsed & !is_special & grepl("^(6+|7+|8+|9+)$", sub("^0+", "", tok)))
print(head(sort(table(CODES$code_label[rep_digit]), decreasing = TRUE), 60))
write.csv(CODES, file.path(out, "ewas_codes.csv"), row.names = FALSE)
write.csv(VARS, file.path(out, "ewas_vars.csv"), row.names = FALSE)
write.csv(CHECK, file.path(out, "ewas_dict_check.csv"), row.names = FALSE)
cat("\n-> ewas_codes.csv, ewas_vars.csv, ewas_dict_check.csv\n")
