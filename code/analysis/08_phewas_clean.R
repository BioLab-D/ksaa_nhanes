suppressMessages({library(DBI); library(duckdb)})
res <- here::here("results")
out <- file.path(res, "analysis")
ew <- file.path(out, "phewas")
d <- readRDS(file.path(out, "analysis_data.rds"))
V <- read.csv(file.path(ew, "phewas_vars.csv"))
C <- read.csv(file.path(ew, "phewas_codes.csv"))
K <- read.csv(file.path(ew, "phewas_dict_check.csv"))
V$section[is.na(V$section)] <- ""
V$subsection[is.na(V$subsection)] <- ""
C$code_label[is.na(C$code_label)] <- ""
V$label <- ifelse(is.na(V$sas_label) | !nzchar(V$sas_label), V$index_label, V$sas_label)

# ---------- 1. Candidates and exclusions ----------
V <- V[V$column_type %in% c("DOUBLE", "INTEGER", "TINYINT", "SMALLINT", "BIGINT"), ]
V <- V[!(V$file != "adult" & V$variable %in% V$variable[V$file == "adult"]), ]          # keep the shared header from adult only
V <- V[!(V$file == "lab" & V$variable %in% V$variable[V$file == "exam"]), ]
V$sec <- paste(V$section, "|", V$subsection)
V$excl <- ""
ex <- function(cond, why) V$excl[V$excl == "" & cond] <<- why
ex(grepl("^WT", V$variable) | grepl("REPLICATE|SAMPLING WEIGHTS", V$sec), "survey weight")
ex(grepl("^SDP", V$variable) | grepl("SURVEY DESIGN", V$sec), "survey design")
ex(V$variable %in% c("SEQN", "DMPFSEQ", "DMPSTAT", "HSSEX", "DMARETHN", "DMARACER", "DMAETHNR", "HSAGEIR", "HSAGEU",
                     "HSAITMOR", "DMPCNTYR", "DMPFIPSR", "DMPCREGN", "HSDOIMO", "HYAITMO", "MXPAXTMR", "HXPAXTMR",
                     "MXPSESSR", "HXPSESSR", "MXPTIMO", "HXPTIMO", "MXPTIDW", "HXPTIDW"),
   "identifier, covariate or exam timing")
ex(grepl("NAME/SSN", V$sec), "identifier")
ex(grepl("DENTAL CARIES|PERIODONTAL ASSESSMENTS|RESTORATIONS AND TOOTH CONDITIONS|TYMPANOGRAM", V$sec), "tooth-, site- or point-level item")
ex(grepl("CONDITIONS AFFECTING|INTRODUCTORY INFORMATION|SCREENING QUESTIONS|PHLEBOTOMY SCREENING|\\| GENERAL INFORMATION$|RESPONDENT", V$sec) &
   !(V$variable %in% c("DMPPIR", "HSFSIZER", "HSHSIZER", "DMPMETRO")), "administrative / pre-analytic section")
# variables with only an SI version (e.g. ALP) are kept
ex(grepl("SI$", V$variable) & sub("SI$", "", V$variable) %in% V$variable, "SI-unit duplicate")
ex(grepl("check item|interviewer|examiner|technician|recorder|comment|\\bflag\\b|reliab|language of|version|proxy|who (answered|completed)|reason|why (not|did not)|\\bstatus\\b|session|date of|time of day|day of (the )?week|today's date|\\bmonth of\\b|\\byear of (exam|interview)|serial|machine|calibration|\\bcode\\b",
         V$label, ignore.case = TRUE), "administrative field (label)")
bad <- K[with(K, !ok_all & (n_uncovered_values > 0.01 * n_rows | abs(sum_counts - n_rows) > 0.01 * n_rows)), ]
ex(paste(V$file, V$variable) %in% paste(bad$file, bad$variable), "codebook/data mismatch > 1%")
# manual exclusions after reviewing each candidate label
mk <- function(vars, why) setNames(rep(why, length(vars)), vars)
EXCL_MANUAL <- c(
  mk(c("HAA1", "HAA2", "HAN8", "FPP1228", "DEPEXFLG", "DEPLQUAD", "DEPUQUAD", "DEPCDRQ1", "DEPSMRQ1"), "survey/exam process item"),
  mk(c("HAZA2", "HAZA3", "HAZA4", "HAZA7A", "HAZA12", "PEP6A", "PEP6B", "PEP6F"), "BP procedure (cuff, arm, inflation level, pre-measurement check)"),
  mk(c("HAZNOK1R", "HAZNOK5R", "PEPNOK1R", "PEPNOK5R"), "number of BP readings averaged"),
  mk(c("HAZA11A", "HAZA11B"), "unedited duplicate of HAZA11AR/HAZA11BR (years 5-6 only)"),
  mk(c("HAZA5"), "30-second count duplicate of HAZA5R"),
  mk(c("BMPHTIN", "BMPWTLBS", "RUPUNIT"), "same measurement in other units"),
  mk(c("BDPD0", "BDPK", "BDPSCAN", "BDPSIDE"), "DXA calibration / QC / side scanned"),
  mk(c("SPPMANEU", "SPPRELIA", "SPPTEMP", "SPPTRIAL"), "spirometry QC / device"),
  mk(c("GHPMETH", "I1P2PFLG", "G1PTIM1", "G1PTIM2"), "assay method / kit / draw timing"),
  mk(c("GUPR1ID", "GUPR1CNF", "GUPR1QLT", "GUPR1DX1", "GUPTDX1R"),
      "reader ID / confidence / quality; preliminary reads (final adjudicated GUPFDX1R kept)"),
  mk(c("DEPSTDX1"), "lesion diagnosis of first lesion only (conditional descriptor)"),
  mk(grep("^HFH", V$variable, value = TRUE), "family reference person, not the participant (age, sex, ID or duplicate of own item)"))
ex(V$variable %in% names(EXCL_MANUAL), "manual")
V$excl_note <- ifelse(V$variable %in% names(EXCL_MANUAL), EXCL_MANUAL[V$variable], "")
# excluded by the rules but restored: urban/rural, orientation (date, weekday, address), observed mobility, marital status,
#   diet changed for illness, edentulism, physician-rated health, insulin use, reason for no food
RESCUE <- c("DMPMETRO", "HAA4", "HAA5", "HAA6D", "HAM14", grep("^HAM15", V$variable, value = TRUE), grep("^HAY11", V$variable, value = TRUE),
            "HFA12", "HFF6A", "HFF6B", "HFF6C", "HFF6D", "DEPEDENT", "PEP13A", "PHPINSU")
V$excl_note[V$variable %in% RESCUE & V$excl != ""] <- paste("rescued from rule:", V$excl[V$variable %in% RESCUE & V$excl != ""])
V$excl[V$variable %in% RESCUE] <- ""

# ---------- 2. Cohort data ----------
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
dbWriteTable(con, "coh", data.frame(SEQN = d$SEQN), temporary = TRUE)
RAW <- data.frame(SEQN = d$SEQN)
for (f in c("adult", "exam", "lab")) {
  vs <- V$variable[V$file == f & V$excl == ""]
  x <- dbGetQuery(con, sprintf('select SEQN, %s from nhanes3_%s join coh using (SEQN)', paste0('"', vs, '"', collapse = ", "), f))
  RAW <- merge(RAW, x, by = "SEQN", all.x = TRUE)
}
dbDisconnect(con, shutdown = TRUE)
RAW <- RAW[match(d$SEQN, RAW$SEQN), ]
stopifnot(identical(RAW$SEQN, d$SEQN))

# ---------- 3. Per-variable cleaning ----------
LESS_RE <- "^(less|fewer) than (1|one)\\b|less than the lowest possible"
NA_RE <- "can'?t do|unable to do|too large for caliper|never found pulse|haven'?t started|no telephone|could not be assessed|can'?t grade"
CAT_NA_RE <- "can'?t grade|could not be assessed|never found pulse|too large for caliper|no telephone"
SPECIAL_RE <- "blank but applicable|don'?t know|^refused\\b|not ascertained|cannot be determined|could not be (determined|assessed)|varies or varied|^unknown$"
TIME_RE <- "how long|how old|age (when|at|started|1st|first|you)|years since|months since|weeks since|when did|ago"
eq <- function(x, v) !is.na(x) & abs(x - v) < 1e-9
clean_one <- function(x, cc, vlab) {
  single <- cc[!cc$is_range & !cc$is_blank & !is.na(cc$lo), ]
  rng <- cc[cc$is_range & !is.na(cc$lo), ]
  covered <- is.na(x) | Reduce(`|`, lapply(single$lo, function(v) eq(x, v)), FALSE) |
    Reduce(`|`, lapply(seq_len(nrow(rng)), function(r) !is.na(x) & x >= rng$lo[r] - 1e-9 & x <= rng$hi[r] + 1e-9), FALSE)
  n_unc <- sum(!covered)
  x[!covered] <- NA
  spc <- single[grepl(SPECIAL_RE, single$code_label, ignore.case = TRUE), ]
  n_spc <- sum(Reduce(`|`, lapply(spc$lo, function(v) eq(x, v)), FALSE))
  for (v in spc$lo) x[eq(x, v)] <- NA
  good_rng <- rng[!grepl(SPECIAL_RE, rng$code_label, ignore.case = TRUE) & !is.na(rng$count) & rng$count > 0, ]
  lab_single <- single[!grepl(SPECIAL_RE, single$code_label, ignore.case = TRUE), ]
  unlab <- sum(!nzchar(lab_single$code_label))
  # continuous: a range spanning >= 4 or with decimals (BMD, calcium, waist-hip ratio span < 4), or >= 5 unlabelled values
  cont <- any(good_rng$hi - good_rng$lo >= 4 | grepl("\\.", good_rng$code_raw)) || (unlab >= 5 && unlab >= 0.8 * nrow(lab_single))
  n_rec <- 0
  if (cont) {
    for (r in which(grepl(LESS_RE, lab_single$code_label, ignore.case = TRUE) |
                    grepl("less than one half hour", lab_single$code_label, ignore.case = TRUE))) {
      lb <- lab_single$code_label[r]
      v <- lab_single$lo[r]
      nv <- if (grepl("half hour", lb, ignore.case = TRUE)) 0.25
            else if (grepl("lowest possible", lb, ignore.case = TRUE)) min(good_rng$lo) / 2
            else 0.5
      n_rec <- n_rec + sum(eq(x, v))
      x[eq(x, v)] <- nv
    }
    for (r in which(grepl(NA_RE, lab_single$code_label, ignore.case = TRUE) |
                    (grepl("^never\\b|none/never", lab_single$code_label, ignore.case = TRUE) & grepl(TIME_RE, vlab, ignore.case = TRUE)))) {
      n_spc <- n_spc + sum(eq(x, lab_single$lo[r]))
      x[eq(x, lab_single$lo[r])] <- NA
    }
  } else {
    # categorical: "Unable to do" is a real outcome level and is kept; "Can't grade", "could not be assessed" etc. -> missing
    for (r in which(grepl(CAT_NA_RE, lab_single$code_label, ignore.case = TRUE))) {
      n_spc <- n_spc + sum(eq(x, lab_single$lo[r]))
      x[eq(x, lab_single$lo[r])] <- NA
    }
    lab_single <- lab_single[!grepl(CAT_NA_RE, lab_single$code_label, ignore.case = TRUE), ]
  }
  # drop skip-pattern hints such as "(HAS26)"
  labs <- setNames(trimws(sub("\\s*\\(.*$", "", lab_single$code_label)),
                   format(lab_single$lo, scientific = FALSE, trim = TRUE))
  list(x = x, cont = cont, n_uncovered = n_unc, n_special = n_spc, n_recoded = n_rec, labs = labs)
}

# ---------- 4. Per-variable counts ----------
cand <- V[V$excl == "", ]
XR <- list()
AUD <- list()
LABS <- list()
CONT <- list()
for (i in seq_len(nrow(cand))) {
  v <- cand$variable[i]
  f <- cand$file[i]
  cc <- C[C$file == f & C$variable == v, ]
  cl <- clean_one(RAW[[v]], cc, cand$label[i])
  x <- cl$x
  # counts only, no fitting or effect sizes; Blank = not applicable/not asked, Blank but applicable = applicable but missing
  raw <- RAW[[v]]
  single <- cc[!cc$is_range & !cc$is_blank & !is.na(cc$lo), ]
  cnt <- function(re) {
    s0 <- single[grepl(re, single$code_label, ignore.case = TRUE), ]
    if (!nrow(s0)) 0L else sum(Reduce(`|`, lapply(s0$lo, function(u) eq(raw, u)), FALSE))
  }
  ux <- unique(x[!is.na(x)])
  ctype <- if (cl$cont && length(ux) >= 5) "continuous" else "categorical"
  tb <- if (ctype == "categorical") table(x) else NULL
  lv <- if (ctype == "categorical") {
    paste(sprintf("%s=%d", ifelse(names(tb) %in% names(cl$labs), cl$labs[names(tb)], names(tb)), as.integer(tb)), collapse = " | ")
  } else ""
  AUD[[v]] <- data.frame(
    file = f, variable = v, label = cand$label[i], section = cand$section[i], subsection = cand$subsection[i],
    type = ctype, n_levels = if (ctype == "categorical") length(ux) else NA, n_distinct = length(ux),
    n_cohort = length(raw), n_blank = sum(is.na(raw)), n_blank_applicable = cnt("blank but applicable"),
    n_dk_refused = cnt("don'?t know|^refused|not ascertained|cannot be determined|could not be (determined|assessed)"),
    n_cant_do = cnt(NA_RE), n_uncovered = cl$n_uncovered, n_recoded = cl$n_recoded, n_valid = sum(!is.na(x)),
    n_valid_phase1 = sum(!is.na(x) & d$phase == 1), n_valid_phase2 = sum(!is.na(x) & d$phase == 2),
    min_level_n = if (ctype == "categorical") min(as.integer(tb)) else NA, levels = lv, stringsAsFactors = FALSE)
  XR[[v]] <- cl$x   # cleaned raw values, for 10_phewas.R
  LABS[[v]] <- cl$labs
  CONT[[v]] <- cl$cont
}
A <- do.call(rbind, AUD)
rownames(A) <- NULL
A$n_applicable <- A$n_cohort - A$n_blank
A$validity_applicable <- A$n_valid / A$n_applicable
A$validity_cohort <- A$n_valid / A$n_cohort
write.csv(A, file.path(ew, "phewas_audit_counts.csv"), row.names = FALSE)
saveRDS(list(SEQN = d$SEQN, x = XR, labs = LABS, cont = CONT), file.path(ew, "phewas_clean.rds"))
cat(sprintf("Counts for %d candidate variables -> phewas_audit_counts.csv\n", nrow(A)))
