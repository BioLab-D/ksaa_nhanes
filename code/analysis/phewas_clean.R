# EWAS cleaning and scan. EWAS_STAGE=audit (used by the release pipeline): only per-variable counts and cleaned phenotypes for phewas_audit.R / phewas.R; no models fitted.
# Default run stage: KSAA (gap, years) ~ variable + age + sex + race, unweighted OLS; continuous variables per SD (log-transformed first if skewed), with standardized beta.
#   Categorical variables (>=3 levels) use an overall F test; all variables also get partial R^2 (relative drop in residual SS when added), an effect size comparable across types.
# Missing codes: handled per variable using the codes read from the official codebook by phewas_dict.R:
#   (1) values not listed in the codebook -> missing ("Don't know" in several *S variables is 9998 in the data, not 9999 as in the codebook);
#   (2) labels such as Blank but applicable / Don't know / Refused / Could not be assessed -> missing;
#   (3) 7-fill "Less than 1 ..." in continuous variables -> 0.5 units ("less than one half hour" -> 0.25; "less than the lowest possible" -> half the lower bound);
#   (4) "Never" in continuous variables -> 0 (frequency items), but missing in "how long ago / at what age" items; "Says can't do", "Skinfold too large", etc. -> missing.
# Inclusion: numeric DuckDB columns; the demographic header shared by the three files taken once (adult); >= 1000 non-missing in the cohort; continuous >= 5 distinct values; categorical >= 30 per level (smaller levels set to missing, >= 2 levels left).
# Exclusions: sampling weights / design variables / administrative and process fields / tooth-surface and periodontal-site codes / SI-unit duplicates / the covariates themselves / codebook check failures, plus a manual list (EXCL_MANUAL, with reasons).
# Multiple testing: BH-FDR over all tests; Bonferroni and the Li-Ji effective number of tests also reported; FDR-significant variables greedily merged into independent signals at |r| > 0.7.
# Input results/analysis/analysis_data.rds, results/analysis/phewas/ewas_{vars,codes,dict_check}.csv, nhanes_lab/data/nhanes.duckdb
# Output results/analysis/phewas/: audit stage ewas_audit_counts.csv, ewas_clean.rds; run stage ewas_candidates.csv, ewas_{results,domains,top}{_main,}.csv (_main excludes fractures, falls and physical function), ewas_cor.rds, ewas_clean_top.rds
# Run by run_all.sh with EWAS_STAGE=audit; candidate list only: EWAS_STAGE=list Rscript code/analysis/phewas_clean.R
suppressMessages({library(DBI); library(duckdb)})
res <- here::here("results"); out <- file.path(res, "analysis"); ew <- file.path(out, "phewas")
STAGE <- Sys.getenv("EWAS_STAGE", "run")
fp <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))
d <- readRDS(file.path(out, "analysis_data.rds")); SDG <- sd(d$gap)
V <- read.csv(file.path(ew, "ewas_vars.csv")); C <- read.csv(file.path(ew, "ewas_codes.csv")); K <- read.csv(file.path(ew, "ewas_dict_check.csv"))
V$section[is.na(V$section)] <- ""; V$subsection[is.na(V$subsection)] <- ""; C$code_label[is.na(C$code_label)] <- ""
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
ex(V$variable %in% c("SEQN", "DMPFSEQ", "DMPSTAT", "HSSEX", "DMARETHN", "DMARACER", "DMAETHNR", "HSAGEIR", "HSAGEU", "HSAITMOR", "DMPCNTYR", "DMPFIPSR", "DMPCREGN",
                     "HSDOIMO", "HYAITMO", "MXPAXTMR", "HXPAXTMR", "MXPSESSR", "HXPSESSR", "MXPTIMO", "HXPTIMO", "MXPTIDW", "HXPTIDW"),
   "identifier, covariate or exam timing")
ex(grepl("NAME/SSN", V$sec), "identifier")
ex(grepl("DENTAL CARIES|PERIODONTAL ASSESSMENTS|RESTORATIONS AND TOOTH CONDITIONS|TYMPANOGRAM", V$sec), "tooth-, site- or point-level item")
ex(grepl("CONDITIONS AFFECTING|INTRODUCTORY INFORMATION|SCREENING QUESTIONS|PHLEBOTOMY SCREENING|\\| GENERAL INFORMATION$|RESPONDENT", V$sec) &
   !(V$variable %in% c("DMPPIR", "HSFSIZER", "HSHSIZER", "DMPMETRO")), "administrative / pre-analytic section")
ex(grepl("SI$", V$variable) & sub("SI$", "", V$variable) %in% V$variable, "SI-unit duplicate")      # variables with only an SI version (e.g. ALP) are kept
ex(grepl("check item|interviewer|examiner|technician|recorder|comment|\\bflag\\b|reliab|language of|version|proxy|who (answered|completed)|reason|why (not|did not)|\\bstatus\\b|session|date of|time of day|day of (the )?week|today's date|\\bmonth of\\b|\\byear of (exam|interview)|serial|machine|calibration|\\bcode\\b",
         V$label, ignore.case = TRUE), "administrative field (label)")
bad <- K[with(K, !ok_all & (n_uncovered_values > 0.01 * n_rows | abs(sum_counts - n_rows) > 0.01 * n_rows)), ]
ex(paste(V$file, V$variable) %in% paste(bad$file, bad$variable), "codebook/data mismatch > 1%")
# manual list (added after reviewing each candidate label; reasons given)
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
  mk(c("GUPR1ID", "GUPR1CNF", "GUPR1QLT", "GUPR1DX1", "GUPTDX1R"), "reader ID / confidence / quality; preliminary reads (final adjudicated GUPFDX1R kept)"),
  mk(c("DEPSTDX1"), "lesion diagnosis of first lesion only (conditional descriptor)"),
  mk(grep("^HFH", V$variable, value = TRUE), "family reference person, not the participant (age, sex, ID or duplicate of own item)"))
ex(V$variable %in% names(EXCL_MANUAL), "manual")
V$excl_note <- ifelse(V$variable %in% names(EXCL_MANUAL), EXCL_MANUAL[V$variable], "")
# wrongly excluded by the rules, restored manually: urban/rural, orientation (today's date / day of week / address), interviewer-observed mobility, marital status, diet changed for illness, edentulism, physician-rated health, insulin use, reason for no food
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
RAW <- RAW[match(d$SEQN, RAW$SEQN), ]; stopifnot(identical(RAW$SEQN, d$SEQN))

# ---------- 3. Per-variable cleaning ----------
LESS_RE <- "^(less|fewer) than (1|one)\\b|less than the lowest possible"
NA_RE <- "can'?t do|unable to do|too large for caliper|never found pulse|haven'?t started|no telephone|could not be assessed|can'?t grade"
CAT_NA_RE <- "can'?t grade|could not be assessed|never found pulse|too large for caliper|no telephone"
SPECIAL_RE <- "blank but applicable|don'?t know|^refused\\b|not ascertained|cannot be determined|could not be (determined|assessed)|varies or varied|^unknown$"
TIME_RE <- "how long|how old|age (when|at|started|1st|first|you)|years since|months since|weeks since|when did|ago"
eq <- function(x, v) !is.na(x) & abs(x - v) < 1e-9
clean_one <- function(x, cc, vlab) {
  single <- cc[!cc$is_range & !cc$is_blank & !is.na(cc$lo), ]; rng <- cc[cc$is_range & !is.na(cc$lo), ]
  covered <- is.na(x) | Reduce(`|`, lapply(single$lo, function(v) eq(x, v)), FALSE) |
    Reduce(`|`, lapply(seq_len(nrow(rng)), function(r) !is.na(x) & x >= rng$lo[r] - 1e-9 & x <= rng$hi[r] + 1e-9), FALSE)
  n_unc <- sum(!covered); x[!covered] <- NA
  spc <- single[grepl(SPECIAL_RE, single$code_label, ignore.case = TRUE), ]
  n_spc <- sum(Reduce(`|`, lapply(spc$lo, function(v) eq(x, v)), FALSE)); for (v in spc$lo) x[eq(x, v)] <- NA
  good_rng <- rng[!grepl(SPECIAL_RE, rng$code_label, ignore.case = TRUE) & !is.na(rng$count) & rng$count > 0, ]
  lab_single <- single[!grepl(SPECIAL_RE, single$code_label, ignore.case = TRUE), ]
  unlab <- sum(!nzchar(lab_single$code_label))
  # continuous: a range with integer span >= 4, or a range with decimals (BMD, calcium mmol/L and waist-hip ratio span < 4); or >= 5 unlabelled enumerated values
  cont <- any(good_rng$hi - good_rng$lo >= 4 | grepl("\\.", good_rng$code_raw)) || (unlab >= 5 && unlab >= 0.8 * nrow(lab_single))
  n_rec <- 0
  if (cont) {
    for (r in which(grepl(LESS_RE, lab_single$code_label, ignore.case = TRUE) | grepl("less than one half hour", lab_single$code_label, ignore.case = TRUE))) {
      lb <- lab_single$code_label[r]; v <- lab_single$lo[r]
      nv <- if (grepl("half hour", lb, ignore.case = TRUE)) 0.25 else if (grepl("lowest possible", lb, ignore.case = TRUE)) min(good_rng$lo) / 2 else 0.5
      n_rec <- n_rec + sum(eq(x, v)); x[eq(x, v)] <- nv
    }
    for (r in which(grepl(NA_RE, lab_single$code_label, ignore.case = TRUE) |
                    (grepl("^never\\b|none/never", lab_single$code_label, ignore.case = TRUE) & grepl(TIME_RE, vlab, ignore.case = TRUE)))) {
      n_spc <- n_spc + sum(eq(x, lab_single$lo[r])); x[eq(x, lab_single$lo[r])] <- NA
    }
  } else {
    # categorical: "Unable to do / Can't do" is a real outcome level and is kept; "Can't grade / could not be assessed / never found pulse" means not measurable -> missing
    for (r in which(grepl(CAT_NA_RE, lab_single$code_label, ignore.case = TRUE))) {
      n_spc <- n_spc + sum(eq(x, lab_single$lo[r])); x[eq(x, lab_single$lo[r])] <- NA
    }
    lab_single <- lab_single[!grepl(CAT_NA_RE, lab_single$code_label, ignore.case = TRUE), ]
  }
  labs <- setNames(trimws(sub("\\s*\\(.*$", "", lab_single$code_label)), format(lab_single$lo, scientific = FALSE, trim = TRUE))   # drop skip-pattern hints such as "(HAS26)"
  list(x = x, cont = cont, n_uncovered = n_unc, n_special = n_spc, n_recoded = n_rec, labs = labs)
}
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3

# ---------- 4. Per-variable models ----------
cand <- V[V$excl == "", ]
base <- data.frame(gap = d$gap, age = d$age, sex = factor(d$sex), race = factor(d$race))
R <- list(); X <- list(); XR <- list(); FILT <- character(); AUD <- list(); LABS <- list(); CONT <- list()   # AUD: per-variable counts for the audit stage (EWAS_STAGE=audit)
for (i in seq_len(nrow(cand))) {
  v <- cand$variable[i]; f <- cand$file[i]; cc <- C[C$file == f & C$variable == v, ]
  cl <- clean_one(RAW[[v]], cc, cand$label[i]); x <- cl$x
  if (STAGE == "audit") {   # counts only, no fitting or effect sizes; Blank = not applicable/not asked, Blank but applicable = applicable but missing
    raw <- RAW[[v]]; single <- cc[!cc$is_range & !cc$is_blank & !is.na(cc$lo), ]
    cnt <- function(re) { s0 <- single[grepl(re, single$code_label, ignore.case = TRUE), ]; if (!nrow(s0)) 0L else sum(Reduce(`|`, lapply(s0$lo, function(u) eq(raw, u)), FALSE)) }
    ux <- unique(x[!is.na(x)]); ctype <- if (cl$cont && length(ux) >= 5) "continuous" else "categorical"
    tb <- if (ctype == "categorical") table(x) else NULL
    lv <- if (ctype == "categorical") paste(sprintf("%s=%d", ifelse(names(tb) %in% names(cl$labs), cl$labs[names(tb)], names(tb)), as.integer(tb)), collapse = " | ") else ""
    AUD[[v]] <- data.frame(file = f, variable = v, label = cand$label[i], section = cand$section[i], subsection = cand$subsection[i],
      type = ctype, n_levels = if (ctype == "categorical") length(ux) else NA, n_distinct = length(ux),
      n_cohort = length(raw), n_blank = sum(is.na(raw)), n_blank_applicable = cnt("blank but applicable"),
      n_dk_refused = cnt("don'?t know|^refused|not ascertained|cannot be determined|could not be (determined|assessed)"),
      n_cant_do = cnt(NA_RE), n_uncovered = cl$n_uncovered, n_recoded = cl$n_recoded, n_valid = sum(!is.na(x)),
      n_valid_phase1 = sum(!is.na(x) & d$phase == 1), n_valid_phase2 = sum(!is.na(x) & d$phase == 2),
      min_level_n = if (ctype == "categorical") min(as.integer(tb)) else NA, levels = lv, stringsAsFactors = FALSE)
    XR[[v]] <- cl$x; LABS[[v]] <- cl$labs; CONT[[v]] <- cl$cont      # cleaned raw values, for phewas.R
    next
  }
  ok <- !is.na(x); n <- sum(ok)
  if (n < 1000) { FILT[v] <- sprintf("n %d < 1000 after cleaning", n); next }
  s <- cbind(base, x = x)[ok, ]
  type <- if (cl$cont) "continuous" else "categorical"
  if (type == "continuous" && length(unique(s$x)) < 5) type <- "categorical"
  logged <- FALSE
  if (type == "continuous") {
    if (min(s$x) >= 0 && skew(s$x) > 1) { c0 <- if (min(s$x) > 0) 0 else min(s$x[s$x > 0]) / 2; s$x <- log(s$x + c0); logged <- TRUE }
    s$x <- (s$x - mean(s$x)) / sd(s$x)
    s$xi <- qnorm((rank(s$x, ties.method = "average") - 0.5) / nrow(s))                 # sensitivity: rank-based inverse normal
  } else {
    tb <- table(s$x); keep <- names(tb)[tb >= 30]
    if (length(keep) < 2) { FILT[v] <- "fewer than 2 levels with >= 30"; next }
    s <- s[as.character(s$x) %in% keep, ]; n <- nrow(s)
    if (n < 1000) { FILT[v] <- sprintf("n %d < 1000 after dropping sparse levels", n); next }
    # reference level: labels such as No / None / Never / No difficulty / Normal first, otherwise the smallest code (NHANES III ordinal scales mostly use 1 = none / mildest)
    lv <- as.character(sort(as.numeric(keep))); ll <- cl$labs[lv]
    neg <- which(grepl("^(no|none|never|normal|not)\\b|no difficulty|no findings|no pain", ll, ignore.case = TRUE))
    ref <- if (length(neg)) lv[neg[1]] else lv[1]
    s$x <- relevel(factor(as.character(s$x), levels = lv), ref = ref)
  }
  cov <- paste(c("age", if (length(unique(s$sex)) > 1) "sex", if (length(unique(s$race)) > 1) "race"), collapse = " + ")   # single-sex items (reproductive health) are not adjusted for sex
  m0 <- lm(as.formula(paste("gap ~", cov)), s); m1 <- lm(as.formula(paste("gap ~ x +", cov)), s)
  a <- anova(m0, m1); p <- a$`Pr(>F)`[2]; pr2 <- (deviance(m0) - deviance(m1)) / deviance(m0)
  cf <- summary(m1)$coefficients; xr <- grep("^x", rownames(cf))
  lab_of <- function(code) { l <- cl$labs[code]; ifelse(is.na(l) | !nzchar(l), code, paste0(code, " ", l)) }
  row <- data.frame(file = f, variable = v, label = cand$label[i], section = cand$section[i], subsection = cand$subsection[i], type = type,
    levels = if (type == "categorical") nlevels(s$x) else NA, logged = logged, n = n, n_uncovered = cl$n_uncovered, n_special = cl$n_special, n_recoded = cl$n_recoded,
    contrast = NA, ref_level = NA, top_level = NA, beta_y = NA, lo_y = NA, hi_y = NA, std_beta = NA, partial_r = NA, partial_r2 = pr2, p = p, p_int = NA, beta_int = NA)
  if (type == "continuous" || nlevels(s$x) == 2) {
    ci <- confint(m1)[xr, ]
    row$beta_y <- cf[xr, 1]; row$lo_y <- ci[1]; row$hi_y <- ci[2]; row$std_beta <- cf[xr, 1] / SDG; row$partial_r <- sign(cf[xr, 1]) * sqrt(pr2)
    row$contrast <- if (type == "continuous") (if (logged) "per SD (log)" else "per SD") else sprintf("%s vs %s", lab_of(levels(s$x)[2]), lab_of(levels(s$x)[1]))
    if (type != "continuous") { row$ref_level <- lab_of(levels(s$x)[1]); row$top_level <- lab_of(levels(s$x)[2]) }
  } else {
    # multi-level: overall F test is primary; the level most different from the reference and its 95% CI are also recorded (descriptive; beta_y etc. left blank to avoid confusion with per-SD estimates)
    k <- which.max(abs(cf[xr, 1])); ci <- confint(m1)[xr[k], ]
    row$ref_level <- lab_of(levels(s$x)[1]); row$top_level <- lab_of(sub("^x", "", rownames(cf)[xr][k]))
    row$top_beta <- cf[xr[k], 1]; row$top_lo <- ci[1]; row$top_hi <- ci[2]
    row$contrast <- sprintf("%d levels; ref %s; largest %s %+.2f y", nlevels(s$x), row$ref_level, row$top_level, cf[xr[k], 1])
  }
  if (type == "continuous") { mi <- lm(as.formula(paste("gap ~ xi +", cov)), s); row$beta_int <- coef(mi)[["xi"]]; row$p_int <- summary(mi)$coefficients["xi", 4] }
  for (nm in c("top_beta", "top_lo", "top_hi")) if (is.null(row[[nm]])) row[[nm]] <- NA
  R[[v]] <- row
  # numeric version for the correlation matrix (row = cohort row; multi-level dichotomized as "non-reference level", used only for merging and the effective number of tests)
  z <- rep(NA_real_, nrow(d)); z[as.integer(rownames(s))] <- if (type == "continuous") s$x else as.numeric(s$x != levels(s$x)[1]); X[[v]] <- z
  XR[[v]] <- cl$x                                                                            # cleaned, untransformed raw values
}
if (STAGE == "audit") {
  A <- do.call(rbind, AUD); rownames(A) <- NULL
  A$n_applicable <- A$n_cohort - A$n_blank; A$validity_applicable <- A$n_valid / A$n_applicable; A$validity_cohort <- A$n_valid / A$n_cohort
  write.csv(A, file.path(ew, "ewas_audit_counts.csv"), row.names = FALSE)
  saveRDS(list(SEQN = d$SEQN, x = XR, labs = LABS, cont = CONT), file.path(ew, "ewas_clean.rds"))
  cat(sprintf("audit: counts for %d candidates within the rules written to ewas_audit_counts.csv (no effect sizes or P values)\n", nrow(A))); quit(save = "no")
}
R <- do.call(rbind, R); rownames(R) <- NULL

cand_out <- V[, c("file", "variable", "label", "section", "subsection", "excl", "excl_note")]
cand_out$analysed <- cand_out$variable %in% R$variable & cand_out$excl == ""
cand_out$filter <- ifelse(cand_out$variable %in% names(FILT) & cand_out$excl == "", FILT[cand_out$variable], "")
write.csv(cand_out, file.path(ew, "ewas_candidates.csv"), row.names = FALSE)
cat(sprintf("Numeric columns %d; excluded by rules %d; cleaned %d; not meeting inclusion after cleaning %d; analyzed %d (continuous %d / binary %d / multi-level %d)\n", nrow(V), sum(V$excl != ""), nrow(cand),
            length(FILT), nrow(R), sum(R$type == "continuous"), sum(R$type == "categorical" & R$levels == 2), sum(R$type == "categorical" & R$levels > 2)))
print(table(V$excl[V$excl != ""]))
if (STAGE == "list") { write.csv(R[, c("file", "variable", "label", "section", "subsection", "type", "levels", "n")], file.path(ew, "ewas_analysable_preview.csv"), row.names = FALSE); quit(save = "no") }

# ---------- 5. Domains ----------
LABDOM <- c(
  setNames(rep("Cardiometabolic biomarkers", 17), c("AAP", "ABP", "CHP", "TCP", "HDP", "LCP", "LPP", "TGP", "TRP", "G1P", "G2P", "SGP", "GHP", "I1P", "I2P", "C1P", "C2P")),
  setNames(rep("Inflammation", 2), c("CRP", "FBP")),
  setNames(rep("Haematology & iron status", 6), c("EPP", "FEP", "SFP", "TIP", "PXP", "FRP")),
  setNames(rep("Nutritional biomarkers", 12), c("ACP", "BCP", "BXP", "LUP", "LYP", "REP", "VAP", "VBP", "VCP", "VEP", "FOP", "RBP")),
  setNames(rep("Nutritional biomarkers", 1), "SEP"),
  setNames(rep("Bone & mineral biochemistry", 5), c("APPSI", "CAPSI", "SCP", "PSP", "ICPSI")),
  setNames(rep("Environmental exposure", 2), c("PBP", "UDP")))
dom_of <- function(file, sec, v) {
  if (file == "lab") {
    if (v %in% names(LABDOM)) return(LABDOM[[v]])
    if (v == "PHPINSU") return("Diabetes history")                          # pre-venipuncture question "currently taking insulin"
    if (grepl("HEMATOLOGY", sec)) return("Haematology & iron status")
    if (grepl("ANTIBODY", sec)) return("Infection & immunity")
    if (grepl("URINE|BIOCHEMISTRY", sec)) return("Liver, kidney & electrolytes")
    return("Other laboratory")
  }
  rules <- list(
    "Physical function & mobility" = c("PHYSICAL FUNCTIONING|PHYSICAL FUNCTION EVALUATION|LOCOMOTION", "^(PEP13E|HAY11|PEP3|HAS24|HAS25|HAH)"),
    "Musculoskeletal (joints, bone, fracture)" = c("MUSCULOSKELETAL|JOINT EXAMINATION|BONE DENSITOMETRY", "^(HAC1A|HAC1M|HAC3AR|HAC6|HAC7|HAM15L|HAM15M)$"),
    "Anthropometry & body composition" = c("BODY MEASUREMENTS", "^(HAM5S|HAM6S|HAM8S|HAM9S|HAM10S|HAM11|HAM12|HAM13|MAPC1|BMP|PEP12)"),
    "Physical activity" = c("EXERCISE", "^HAT"),
    "Tobacco & alcohol" = c("TOBACCO|ALCOHOL", "^(HAR|MAPA|MAPE|HFF1$|HFF2R|HFF3A)"),
    "Diet" = c("SECTION M. DIET|FOOD FREQUENCY", "^(HAM|HAN|MAPC8)"),
    "Cardiovascular & blood pressure" = c("HIGH BLOOD PRESSURE|CARDIOVASCULAR DISEASE|BLOOD PRESSURE MEASUREMENT|\\| BLOOD PRESSURE|\\| HEART", "^(HAC1C|HAC1D|HAC4B|HAE|HAF|HAZ|PEP6|PEPMN|PEP8)"),
    "Diabetes history" = c("SECTION D. DIABETES", "^(HAD|HAC4A|HAC5A|PHPINSU)"),
    "Respiratory" = c("RESPIRATORY|SPIROMETRY|\\| CHEST", "^(HAC1E|HAC1F|HAC1G|HAC1H|HAL|SPP|PEP7)"),
    "Oral health" = c("DENTAL", "^(HAQ[1-4]|DEP)"),
    "Sensory, retina & cognition" = c("VISION AND HEARING|COGNITIVE|FUNDUS", "^(HAP|HAA[4-6]|FPP|MAPB|MAPD|HAC1I)"),
    "Reproductive history" = c("REPRODUCTIVE", "^MAPF"),
    "Health care & general health" = c("HEALTH SERVICES|PHYSICIAN'S IMPRESSION", "^(HAB|PEP13)"),
    "Socioeconomic, social & household" = c("FAMILY QUESTIONNAIRE|OCCUPATION|SOCIAL SUPPORT|^ \\| GENERAL INFORMATION", "^(HF|HAS|HAV|DMP|HS)"))
  for (nm in names(rules)) if (grepl(rules[[nm]][2], v)) return(nm)
  for (nm in names(rules)) if (grepl(rules[[nm]][1], sec)) return(nm)
  "Other medical history & examination"
}
R$domain <- mapply(dom_of, R$file, paste(R$section, "|", R$subsection), R$variable)
# groups: physical function is separate (analyzed as an outcome); bone and joint (including bone-metabolism biochemistry) is another group
GROUP <- c("Musculoskeletal (joints, bone, fracture)" = "Bone & joints", "Physical function & mobility" = "Physical function",
           "Bone & mineral biochemistry" = "Bone & joints", "Anthropometry & body composition" = "Body size & composition",
           "Physical activity" = "Behaviour", "Tobacco & alcohol" = "Behaviour", "Diet" = "Behaviour",
           "Cardiometabolic biomarkers" = "Cardiometabolic & inflammatory", "Inflammation" = "Cardiometabolic & inflammatory",
           "Cardiovascular & blood pressure" = "Cardiometabolic & inflammatory", "Diabetes history" = "Cardiometabolic & inflammatory")
R$group <- ifelse(R$domain %in% names(GROUP), GROUP[R$domain], "Other")

# ---------- 6. Two scan scopes ----------
# the main scan excludes fractures, falls and physical function (analyzed as outcomes); the all-variable scan is a sensitivity analysis
FXF <- c("HAG5A", "HAG5B", "HAG5C", grep("^HAG(6|10)", R$variable, value = TRUE), "HAG14", "HAG15")
R$outcome_item <- R$domain == "Physical function & mobility" | R$variable %in% FXF
M <- do.call(cbind, X[R$variable])
t0 <- Sys.time(); CM <- suppressWarnings(cor(M, use = "pairwise.complete.obs")); CM[is.na(CM)] <- 0; diag(CM) <- 1
cat(sprintf("\nCorrelation matrix %dx%d (%.0f s)\n", ncol(M), ncol(M), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("Cleaning: %d values set to missing as not listed in the codebook (%d variables); missing codes set to missing %d; 7-fill recoded %d\n",
            sum(R$n_uncovered), sum(R$n_uncovered > 0), sum(R$n_special), sum(R$n_recoded)))

show <- function(x, k) {
  x <- head(x, k)
  for (i in seq_len(nrow(x))) with(x[i, ], cat(sprintf("%3d %-9s %-44s %-34s n %4d  %s  pr %.4f  P %s  FDR %s%s\n", rank, variable, substr(label, 1, 44), substr(domain, 1, 34), n,
    if (is.na(beta_y)) sprintf("%-40s", substr(contrast, 1, 40)) else sprintf("%+.2f y (%+.2f, %+.2f) β %+.3f %-12s", beta_y, lo_y, hi_y, std_beta, substr(contrast, 1, 12)),
    partial_r2, fp(p), fp(fdr), if (!is.na(lead) && lead != variable) paste0("  ~", lead) else "")))
}
# multiple testing (BH, Bonferroni, Li-Ji effective number of tests), greedy merging of independent signals at |r| > 0.7, domain and group summaries; all recomputed within the scope
finish <- function(R, title, suffix) {
  cat(sprintf("\n################ %s ################\n", title))
  R$fdr <- p.adjust(R$p, "BH"); R$p_bonf <- p.adjust(R$p, "bonferroni")
  cm <- CM[R$variable, R$variable]
  ev <- abs(eigen(cm, symmetric = TRUE, only.values = TRUE)$values); meff <- sum((ev >= 1) + (ev - floor(ev)))
  R$meff <- meff; R$p_meff <- pmin(1, R$p * meff)
  R <- R[order(R$p), ]; R$rank <- seq_len(nrow(R))
  R$lead <- NA_character_; leads <- character()
  for (j in which(R$fdr < 0.05)) {
    v <- R$variable[j]; r <- if (length(leads)) abs(cm[v, leads]) else numeric()
    if (length(r) && max(r) > 0.7) R$lead[j] <- leads[which.max(r)] else { leads <- c(leads, v); R$lead[j] <- v }
  }
  R$is_lead <- !is.na(R$lead) & R$lead == R$variable
  cat(sprintf("Tests %d (continuous %d / binary %d / multi-level %d); Li-Ji effective tests %.0f; FDR threshold P = %.2e; Meff threshold P = %.1e\n", nrow(R), sum(R$type == "continuous"),
              sum(R$type == "categorical" & R$levels == 2), sum(R$type == "categorical" & R$levels > 2), meff, max(c(0, R$p[R$fdr < 0.05])), 0.05 / meff))
  cat(sprintf("FDR < 0.05: %d (independent signals %d); Meff-corrected < 0.05: %d; Bonferroni < 0.05: %d\n", sum(R$fdr < 0.05), sum(R$is_lead), sum(R$p_meff < 0.05), sum(R$p_bonf < 0.05)))
  cat("\n--- Top 30 (by P; ~ = |r| > 0.7 with an earlier independent signal) ---\n"); show(R, 30)
  cat("\n--- Top 30 independent signals ---\n"); show(R[R$is_lead, ], 30)
  D <- do.call(rbind, lapply(split(R, R$domain), function(x) data.frame(domain = x$domain[1], group = x$group[1], n_tested = nrow(x), n_fdr = sum(x$fdr < 0.05),
    pct_fdr = 100 * mean(x$fdr < 0.05), n_lead = sum(x$is_lead), best_rank = min(x$rank), best_var = x$variable[which.min(x$rank)],
    max_partial_r2 = max(x$partial_r2), median_partial_r2 = median(x$partial_r2))))
  tot_sig <- sum(R$fdr < 0.05); D$enrich_or <- NA; D$enrich_p <- NA
  for (i in seq_len(nrow(D))) {
    a <- D$n_fdr[i]; b <- D$n_tested[i] - a; c <- tot_sig - a; dd <- (nrow(R) - D$n_tested[i]) - c
    ft <- fisher.test(matrix(c(a, b, c, dd), 2)); D$enrich_or[i] <- unname(ft$estimate); D$enrich_p[i] <- ft$p.value
  }
  D <- D[order(D$best_rank), ]
  cat("\n--- Domains (by best rank; enrichment = Fisher test of the FDR-significant proportion) ---\n")
  for (i in seq_len(nrow(D))) with(D[i, ], cat(sprintf("  %-40s tests %3d  FDR sig %3d (%4.1f%%)  independent %2d  best #%-3d %-9s max partial R² %.4f  enrichment OR %.2f P %s\n",
    domain, n_tested, n_fdr, pct_fdr, n_lead, best_rank, best_var, max_partial_r2, enrich_or, fp(enrich_p))))
  LD <- R[R$is_lead, ]; G <- sort(unique(R$group))
  cat("\n--- Group composition (top k independent signals by P; last column = all tests) and rank test (group P ranks vs the rest, Wilcoxon) ---\n")
  for (g in G) { ing <- R$group == g; wp <- if (sum(ing) > 1) wilcox.test(R$rank[ing], R$rank[!ing])$p.value else NA
    cat(sprintf("  %-32s top10 %2d top20 %2d top50 %2d top100 %3d | all %3d (%2.0f%%) | median rank %4.0f vs %4.0f P %s\n", g, sum(head(LD$group, 10) == g),
                sum(head(LD$group, 20) == g), sum(head(LD$group, 50) == g), sum(head(LD$group, 100) == g), sum(ing), 100 * mean(ing),
                median(R$rank[ing]), median(R$rank[!ing]), if (is.na(wp)) "—" else fp(wp))) }
  write.csv(R, file.path(ew, paste0("ewas_results", suffix, ".csv")), row.names = FALSE)
  write.csv(D, file.path(ew, paste0("ewas_domains", suffix, ".csv")), row.names = FALSE)
  write.csv(R[R$is_lead, ], file.path(ew, paste0("ewas_top", suffix, ".csv")), row.names = FALSE)
  invisible(R)
}
RM <- finish(R[!R$outcome_item, ], sprintf("Main scan (Fig. 4b): %d variables excluding fractures, falls and physical function", sum(!R$outcome_item)), "_main")
RA <- finish(R, sprintf("Sensitivity: all %d variables", nrow(R)), "")
saveRDS(CM, file.path(ew, "ewas_cor.rds"))
# independent signals passing the effective-number-of-tests correction (union of lead signals of the full and main scans): cleaned raw values + type, for reverse models (variable ~ KSAA)
TL <- union(RA$variable[RA$is_lead & RA$p_meff < 0.05], RM$variable[RM$is_lead & RM$p_meff < 0.05])
saveRDS(list(SEQN = d$SEQN, X = XR[TL], meta = RA[match(TL, RA$variable), ]), file.path(ew, "ewas_clean_top.rds"))

# ---------- 7. Robustness (all variables) ----------
cat("\n=========== Robustness ===========\n")
cc <- RA[RA$type == "continuous", ]
cat(sprintf("  %d continuous variables: same direction after rank-based inverse normal transformation %d / %d; originally FDR-significant with P < 0.05 under INT: %d / %d\n", nrow(cc),
            sum(sign(cc$beta_y) == sign(cc$beta_int)), nrow(cc), sum(cc$fdr < 0.05 & cc$p_int < 0.05), sum(cc$fdr < 0.05)))
top <- rbind(head(RA, 20), head(RM[!RM$variable %in% head(RA$variable, 20), ], 10))
cat("  Value check (top 20 of the full scan + further top 10 of the main scan; continuous: cleaned quantiles; categorical: counts per level) - residual missing codes?\n")
for (i in seq_len(nrow(top))) {
  v <- top$variable[i]; x <- clean_one(RAW[[v]], C[C$file == top$file[i] & C$variable == v, ], top$label[i])$x
  if (top$type[i] == "continuous") cat(sprintf("    %-9s min %g | p1 %g | median %g | p99 %g | max %g | missing codes removed %d unlisted values %d\n", v, min(x, na.rm = TRUE),
    quantile(x, .01, na.rm = TRUE), median(x, na.rm = TRUE), quantile(x, .99, na.rm = TRUE), max(x, na.rm = TRUE), top$n_special[i], top$n_uncovered[i]))
  else { tb <- table(x); cat(sprintf("    %-9s %s | missing codes removed %d\n", v, paste(sprintf("%s=%d", names(tb), tb), collapse = " "), top$n_special[i])) }
}
cat("  Comparison with the same measures in S16 (this scan beta / S16 beta, years per SD):\n")
s16 <- read.csv(file.path(out, "correlates.csv"))
for (pr in list(c("BDPFNBMD", "Femoral neck BMD"), c("BDPTOBMD", "Total hip BMD"), c("BMPBMI", "BMI"), c("GHP", "HbA1c"), c("BMPWAIST", "Waist circumference"), c("WCP", "White blood cell count")))
  if (pr[1] %in% RA$variable) cat(sprintf("    %-9s %+.2f(n %d)/ %s %+.2f(n %d)\n", pr[1], RA$beta_y[RA$variable == pr[1]], RA$n[RA$variable == pr[1]], pr[2],
    s16$beta_y[s16$factor == pr[2]][1], s16$n[s16$factor == pr[2]][1]))
cat("\n-> ewas_candidates.csv, ewas_{results,domains,top}{_main,}.csv, ewas_cor.rds, ewas_clean_top.rds\n")
