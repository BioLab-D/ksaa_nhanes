# PheWAS: survey-weighted scan of the phenotypes in the audit table (09_phewas_audit.R), KSAA -> phenotype.
# One model per phenotype: KSAA per SD + age, sex, race/ethnicity, education and PIR. Continuous: svyglm Gaussian
#   (log-transformed if skewed, then standardized); binary: quasibinomial; ordinal: svyolr (levels < 30 merged into a
#   neighbour); nominal: reverse joint test, P only.
# Design: WTPFEX6, SDPPSU6, SDPSTRA6; fasting-subsample measures use WTPFSD6; reproductive items in women only.
#   Subsets always via subset(design). P and CI use t with design df (49; 23 for a single phase).
# Replication: each phase fitted with its own design; replicated = same direction and P < 0.05 in both phases.
# Multiple testing: BH-FDR, Bonferroni and Li-Ji, separately for main, supplementary and extended-module sets.
# Input: results/analysis/{analysis_data.rds, fxf_vars.csv}, results/analysis/phewas/{phewas_clean.rds, phewas_audit_table.csv,
#   phewas_multilevel_plan.csv}
# Output: results/analysis/phewas/phewas_results.csv, phewas_y.rds, phewas_summary.txt
suppressMessages({library(survey); library(MASS); library(DBI); library(duckdb)})
options(survey.lonely.psu = "adjust")
res <- here::here("results")
out <- file.path(res, "analysis")
ew <- file.path(out, "phewas")
d <- readRDS(file.path(out, "analysis_data.rds"))
CL <- readRDS(file.path(ew, "phewas_clean.rds"))
stopifnot(identical(CL$SEQN, d$SEQN))
A <- read.csv(file.path(ew, "phewas_audit_table.csv"), stringsAsFactors = FALSE)
A <- A[A$in_pool & A$status %in% c("main", "supplementary", "exposure_module"), ]
PL <- read.csv(file.path(ew, "phewas_multilevel_plan.csv"), stringsAsFactors = FALSE)
CLS <- setNames(PL$class, PL$variable)
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
dbWriteTable(con, "coh", data.frame(SEQN = d$SEQN), temporary = TRUE)
W <- dbGetQuery(con,
    "select SEQN, WTPFSD6, WTPFEX1, WTPFSD1, SDPPSU1, SDPSTRA1, WTPFEX2, WTPFSD2, SDPPSU2, SDPSTRA2 from nhanes3_adult join coh using (SEQN)")
dbDisconnect(con, shutdown = TRUE)
W <- W[match(d$SEQN, W$SEQN), ]
for (nm in setdiff(names(W), "SEQN")) { d[[nm]] <- W[[nm]]; d[[nm]][is.na(d[[nm]])] <- 0 }
stopifnot(all(d$WTPFEX1[d$phase == 1] > 0), all(d$WTPFEX2[d$phase == 2] > 0))
SDG <- sd(d$gap)
d$z <- d$gap / SDG
d$sex <- factor(d$sex)
d$race <- factor(d$race)
d$edu3 <- factor(d$edu3)
COV <- c("age", "sex", "race", "edu3", "pir")
mkdes <- function(w, psu, str, rows) svydesign(id = as.formula(paste0("~", psu)), strata = as.formula(paste0("~", str)),
    weights = as.formula(paste0("~", w)), nest = TRUE, data = d[rows, ])
# designs: 6-year, phase 1, phase 2, fasting subsample; rows aligned with y
DS <- list(all = list(des = mkdes("WTPFEX6", "SDPPSU6", "SDPSTRA6", TRUE), rows = rep(TRUE, nrow(d))),
           p1 = list(des = mkdes("WTPFEX1", "SDPPSU1", "SDPSTRA1", d$phase == 1), rows = d$phase == 1),
           p2 = list(des = mkdes("WTPFEX2", "SDPPSU2", "SDPSTRA2", d$phase == 2), rows = d$phase == 2),
           fall = list(des = mkdes("WTPFSD6", "SDPPSU6", "SDPSTRA6", d$WTPFSD6 > 0), rows = d$WTPFSD6 > 0),
           fp1 = list(des = mkdes("WTPFSD1", "SDPPSU1", "SDPSTRA1", d$WTPFSD1 > 0), rows = d$WTPFSD1 > 0),
           fp2 = list(des = mkdes("WTPFSD2", "SDPPSU2", "SDPSTRA2", d$WTPFSD2 > 0), rows = d$WTPFSD2 > 0))
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3

# ---------- Composite scores: story recall immediate / delayed, 6 segments each -> 0-6; address orientation, 4 items -> 0-4 ----------
lab_of <- function(v) CL$labs[[v]]
score <- function(vs, re) {
  M <- sapply(vs, function(v) {
    x <- CL$x[[v]]
    lb <- lab_of(v)
    yes <- as.numeric(names(lb)[grepl(re, lb, ignore.case = TRUE)])
    ifelse(is.na(x), NA, as.numeric(x %in% yes))
  })
  ifelse(rowSums(is.na(M)) > 0, NA, rowSums(M))
}
COMP <- list(
  MAPB_story_recall_immediate_total = list(items = paste0("MAPB2", LETTERS[1:6]), re = "^yes",
                                           label = "Story recall, immediate (details recalled, 0-6)"),
  MAPD_story_recall_delayed_total = list(items = paste0("MAPD2", LETTERS[1:6]), re = "^yes",
                                         label = "Story recall, delayed (details recalled, 0-6)"),
  HAA_orientation_score = list(items = paste0("HAA6", LETTERS[1:4]), re = "^correct",
                               label = "Orientation: address items correct (0-4)"))
for (nm in names(COMP)) {
  it <- COMP[[nm]]
  stopifnot(all(it$items %in% A$variable))
  CL$x[[nm]] <- score(it$items, it$re)
  CL$cont[[nm]] <- TRUE
  CL$labs[[nm]] <- character()
  r <- A[A$variable == it$items[1], ]
  r$variable <- nm
  r$label <- it$label
  r$type <- "continuous"
  r$summarise_to <- ""
  r$reason <- paste("composite of", length(it$items), "items")
  A <- rbind(A[!(A$variable %in% it$items), ], r)
}

# ---------- Summaries of other repeated measures / multi-item tests, defined as in Table S10 ----------
# derived phenotypes must pass: continuous n >= 500 and >= 5 values; binary n >= 500 and >= 100 per class; ordinal n >= 500
gate <- function(y, type) {
  y <- y[!is.na(y)]
  if (type == "continuous") length(y) >= 500 && length(unique(y)) >= 5
  else if (length(unique(y)) == 2) length(y) >= 500 && min(table(y)) >= 100 else length(y) >= 500
}
add_pheno <- function(nm, y, type, label, tmpl, cls, src, note, labs = character()) {
  stopifnot(length(y) == nrow(d), !(nm %in% A$variable))
  r <- A[A$variable == tmpl, ]
  if (nrow(r) == 0) {
    r <- read.csv(file.path(ew, "phewas_audit_table.csv"), stringsAsFactors = FALSE)
    r <- r[r$variable == tmpl, ]
  }
  stopifnot(nrow(r) == 1)
  r$variable <- nm
  r$label <- label
  r$type <- type
  r$health_class <- cls
  r$summarise_to <- ""
  r$s8_version <- ""
  r$reason <- note
  r$status <- if (gate(y, type)) "main" else "supplementary"
  CL$x[[nm]] <<- y
  CL$cont[[nm]] <<- type == "continuous"
  CL$labs[[nm]] <<- labs
  A <<- rbind(A[!(A$variable %in% src), ], r)
}
xv <- function(v) CL$x[[v]]
lev_is <- function(v, re) {
  ll <- CL$labs[[v]][as.character(CL$x[[v]])]
  ifelse(is.na(CL$x[[v]]), NA, as.numeric(grepl(re, ll, ignore.case = TRUE)))
}
# repeated measures: gait speed from the faster trial (asserted equal to the Table S10 variable), mean steps, mean pulse
sp <- 2.4384 / suppressWarnings(pmin(xv("PFPTWLKA"), xv("PFPTWLKB"), na.rm = TRUE))
sp[!is.finite(sp)] <- NA
stopifnot(identical(is.na(sp), is.na(d$gait)), max(abs(sp - d$gait), na.rm = TRUE) < 1e-9)
add_pheno("PFP_gait_speed", sp, "continuous", "Gait speed (faster of two 8-ft walks, m/s)", "PFPTWLKA", "function",
          c("PFPTWLKA", "PFPTWLKB"), "repeated trials summarised (faster of two), identical to the Table S8 gait speed")
mn2 <- function(a, b) {
  m <- rowMeans(cbind(a, b), na.rm = TRUE)
  m[is.nan(m)] <- NA
  m
}
add_pheno("PFP_steps_8ft_mean", mn2(xv("PFPNSTPA"), xv("PFPNSTPB")), "continuous", "Steps to walk 8 ft (mean of two trials)",
          "PFPNSTPA", "function", c("PFPNSTPA", "PFPNSTPB"), "repeated trials summarised (mean)")
add_pheno("pulse_mean", mn2(xv("HAZA5R"), xv("PEP6DR")), "continuous", "Pulse rate (mean of household and examination readings)",
          "HAZA5R", "objective", c("HAZA5R", "PEP6DR"),
          "repeated measurement summarised (mean of household and MEC readings, as for blood pressure)")
# serial subtraction, 5 items -> number correct (17/14/11/8/5): 55 "can't do" and 99 "don't know" count as wrong, 88 and
#   blank as missing; the score is missing if any item is. Raw codes are read from DuckDB because cleaning drops 55/99.
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
dbWriteTable(con, "coh", data.frame(SEQN = d$SEQN), temporary = TRUE)
H18 <- dbGetQuery(con, "select SEQN, HAP18A, HAP18B, HAP18C, HAP18D, HAP18E from nhanes3_adult join coh using (SEQN)")
dbDisconnect(con, shutdown = TRUE)
H18 <- H18[match(d$SEQN, H18$SEQN), ]
KEY <- c(HAP18A = 17, HAP18B = 14, HAP18C = 11, HAP18D = 8, HAP18E = 5)
M18 <- sapply(names(KEY), function(v) {
  r <- H18[[v]]
  ifelse(is.na(r) | r == 88, NA, as.numeric(r == KEY[[v]]))
})
# valid values in the cleaned table match the raw codes
for (v in names(KEY)) stopifnot(all(xv(v)[!is.na(xv(v))] == H18[[v]][!is.na(xv(v))]))
add_pheno("HAP18_serial_subtraction", ifelse(rowSums(is.na(M18)) > 0, NA, rowSums(M18)), "continuous",
          "Serial subtraction of $3, answers correct (0-5)", "HAP18A", "objective", names(KEY),
          "5 items scored correct / incorrect ('can't do' and 'don't know' = incorrect)")
sc3 <- function(vs) {
  M <- sapply(vs, lev_is, re = "^correct")
  ifelse(rowSums(is.na(M)) > 0, NA, rowSums(M))
}
L03 <- c(`0` = "0 correct", `1` = "1 correct", `2` = "2 correct", `3` = "3 correct")
add_pheno("HAP17_word_registration", sc3(c("HAP17A1", "HAP17A2", "HAP17A3")), "categorical",
          "Three-word registration, first trial (words correct, 0-3)", "HAP17A2", "objective",
          c("HAP17A1", "HAP17A2", "HAP17A3"), "3 items summed (as for story recall); ordinal (fewer than 5 values)", L03)
add_pheno("HAP19_word_recall", sc3(c("HAP19A", "HAP19B", "HAP19C")), "categorical",
          "Three-word recall (words correct, 0-3)", "HAP19A", "objective",
          c("HAP19A", "HAP19B", "HAP19C"), "3 items summed (as for story recall); ordinal (fewer than 5 values)", L03)
# binary recodes matching Table S10 (asserted against the Table S10 variables): difficulty items HAH1-12 and PEP13E1A-E5A
#   -> any difficulty vs none; HAG14 -> >= 1 fall in 12 months; PFPTTAND -> tandem stand < 10 s.
FXV <- read.csv(file.path(out, "fxf_vars.csv"))
FXV <- FXV[match(d$SEQN, FXV$SEQN), ]
cmp <- function(a, b) {
  k <- !is.na(a) & !is.na(b)
  stopifnot(all(a[k] == b[k]))
  invisible(TRUE)
}
recode_bin <- function(v, y, lab0, lab1, note) {
  stopifnot(v %in% A$variable, length(y) == nrow(d))
  CL$x[[v]] <<- as.numeric(y)
  CL$cont[[v]] <<- FALSE
  CL$labs[[v]] <<- setNames(c(lab0, lab1), c("0", "1"))
  i <- which(A$variable == v)
  A$type[i] <<- "categorical"
  A$reason[i] <<- note
  A$status[i] <<- if (gate(as.numeric(y), "categorical")) "main" else "supplementary"
}
DIFF <- c(paste0("HAH", 1:12), paste0("PEP13E", 1:5, "A"))
for (v in DIFF) {
  ll <- CL$labs[[v]][as.character(CL$x[[v]])]
  stopifnot(any(grepl("^no difficulty", ll, ignore.case = TRUE)))
  y <- 1 - lev_is(v, "^no difficulty")
  if (v %in% c("HAH1", "HAH2", "HAH8", "HAH9")) {
    cmp(FXV[[paste0("dif_", v)]], y)
    stopifnot(identical(is.na(FXV[[paste0("dif_", v)]]), is.na(y)))
  }
  recode_bin(v, y, "No difficulty", "Any difficulty", "any difficulty (including unable) vs none, as in Table S8")
}
cmp(d$fall, as.numeric(xv("HAG14") >= 1))
stopifnot(identical(is.na(d$fall), is.na(xv("HAG14"))))
recode_bin("HAG14", d$fall, "No", "Yes", "at least one fall in the past 12 months, as in Table S8")
cmp(FXV$tand_fail, as.numeric(xv("PFPTTAND") < 10))
recode_bin("PFPTTAND", FXV$tand_fail, "No", "Yes", "unable to hold the tandem stand for 10 s, as in Table S8")
cat(sprintf("Main-analysis phenotypes %d (supplementary %d; extended module %d)\n",
    sum(A$status == "main"), sum(A$status == "supplementary"), sum(A$status == "exposure_module")))

# ---------- Outcome type and recoding per variable ----------
# third level = not applicable -> missing; binary
NA3 <- c(HAF2 = "Never walks uphill or hurries", HAL14D = "Doesn't work", HAL17 = "Never received allergy shots")
ORDMAP <- list(HAK4 = c(`3` = 0, `1` = 1, `2` = 2),          # never < 1-2 < 3+ times
               PFPSCOOT = c(`2` = 0, `1` = 1, `3` = 2))      # no < yes(uses arms) < unable
# mislabelled NOM/LAT/YNA in the plan; actually ordinal
ORD_AS_IS <- c("PFPUSTND", "PFPLSER", "PFPRSER", "ANP", "HZP", "MIP", "DEPDQL1", "DEPDQU1")
BIN_RULE <- list(
  DEPLPTYP = list(one = "partial", zero = "full", drop = "unable", note = "partial vs full denture; 'unable to evaluate' -> missing"),
  DEPUPTYP = list(one = "partial", zero = "full", drop = "unable", note = "partial vs full denture; 'unable to evaluate' -> missing"),
  GUPFDX1R = list(one = "gallstone|cholecystectomy", zero = "^normal", drop = "no conclusion|abnormal",
                  note = "gallbladder disease (stones or cholecystectomy) vs normal; inconclusive/other -> missing"),
  HAP2 = list(one = "eyeglasses|contact|both", zero = "neither", drop = "^$", note = "any vision correction vs none"))
# ordinal: levels < 30 merged into a neighbour (end levels inward, middle levels toward the smaller neighbour)
merge_sparse <- function(y, min_n = 30) {
  notes <- character()
  repeat {
    tb <- table(y)
    if (length(tb) <= 2 || min(tb) >= min_n) break
    i <- which.min(tb)
    L <- names(tb)
    j <- if (i == 1) 2 else if (i == length(tb)) i - 1 else if (tb[i - 1] <= tb[i + 1]) i - 1 else i + 1
    a <- L[min(i, j)]
    b <- L[max(i, j)]
    new <- paste(a, b, sep = "+")
    notes <- c(notes, sprintf("merged [%s]+[%s]", a, b))
    y <- as.character(y)
    y[y %in% c(a, b)] <- new
    lv <- L
    lv[lv %in% c(a, b)] <- new
    y <- factor(y, unique(lv), ordered = TRUE)
  }
  list(y = y, notes = notes)
}

prep <- function(v, r) {
  x <- CL$x[[v]]
  lb <- CL$labs[[v]]
  notes <- character()
  typ <- r$type
  y <- x
  if (typ == "continuous") {
    if (min(y, na.rm = TRUE) >= 0 && skew(y[!is.na(y)]) > 1) {
      c0 <- if (min(y, na.rm = TRUE) > 0) 0 else min(y[y > 0], na.rm = TRUE) / 2
      y <- log(y + c0)
      notes <- c(notes, "log")
    }
    y <- (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
    return(list(y = y, model = "gaussian", notes = notes, levels = ""))
  }
  ll <- lb[as.character(x)]
  ll[is.na(x)] <- NA
  k <- length(unique(x[!is.na(x)]))
  cls <- if (v %in% names(CLS)) CLS[[v]] else "ORD"
  if (v %in% names(BIN_RULE)) {
    b <- BIN_RULE[[v]]
    y <- ifelse(is.na(ll), NA, ifelse(grepl(b$one, ll, ignore.case = TRUE), 1, ifelse(grepl(b$zero, ll, ignore.case = TRUE), 0, NA)))
    notes <- b$note
    k <- 2
  } else if (v %in% names(NA3)) {
    y <- ifelse(is.na(ll), NA, ifelse(grepl("^yes", ll, ignore.case = TRUE), 1, ifelse(grepl("^no$", ll, ignore.case = TRUE), 0, NA)))
    notes <- sprintf("'%s' -> missing", NA3[[v]])
    k <- 2
  } else if (v %in% names(ORDMAP)) {
    m <- ORDMAP[[v]]
    y <- unname(m[as.character(x)])
    notes <- "levels re-ordered to a monotone scale"
    k <- length(m)
  } else if (cls == "LAT" && any(grepl("right|left", ll, ignore.case = TRUE))) {
    y <- ifelse(is.na(ll), NA, ifelse(grepl("^no\\b", ll, ignore.case = TRUE), 0, ifelse(grepl("both", ll, ignore.case = TRUE), 2, 1)))
    notes <- "none < unilateral < bilateral"
    k <- 3
  } else if (cls == "NOM" && !(v %in% ORD_AS_IS) && k > 2) {
    y <- factor(ll)
    return(list(y = y, model = "nominal", notes = "nominal: reverse joint test only",
                levels = paste(sprintf("%s=%d", names(table(y)), as.integer(table(y))), collapse = " | ")))
  }
  # reference level: No/None/Never/Normal first, otherwise the smallest code
  if (k == 2) {
    u <- sort(unique(y[!is.na(y)]))
    if (!identical(as.numeric(u), c(0, 1))) {
      l2 <- lb[as.character(u)]
      neg <- which(grepl("^(no|none|never|normal|not)\\b|no difficulty|no findings|no pain|^correct|^yes, recalled", l2, ignore.case = TRUE))
      ref <- if (length(neg)) u[neg[1]] else u[1]
      y <- as.numeric(y != ref)
      notes <- c(notes, sprintf("1 = %s", lb[as.character(setdiff(u, ref))]))
    }
    tb <- table(y)
    return(list(y = y, model = "binomial", notes = notes, levels = sprintf("0=%d | 1=%d", tb[["0"]], tb[["1"]])))
  }
  yo <- factor(y, sort(unique(y[!is.na(y)])), ordered = TRUE)
  levels(yo) <- ifelse(is.na(lb[levels(yo)]) | !nzchar(lb[levels(yo)]), levels(yo), lb[levels(yo)])
  ms <- merge_sparse(yo)
  yo <- ms$y
  notes <- c(notes, ms$notes)
  if (nlevels(yo) == 2) {
    y <- as.numeric(yo) - 1
    tb <- table(y)
    return(list(y = y, model = "binomial", notes = c(notes, "2 levels after merge"),
                levels = sprintf("0=%d | 1=%d", tb[["0"]], tb[["1"]])))
  }
  list(y = yo, model = "ordinal", notes = notes, levels = paste(sprintf("%s=%d", levels(yo), as.integer(table(yo))), collapse = " | "))
}

fit <- function(model, des, cov) {
  v <- des$variables
  # drop covariates with one level in the subset (sex for single-sex items)
  cov <- cov[sapply(cov, function(cv) length(unique(v[[cv]])) > 1)]
  df <- degf(des)
  tp <- function(b, se) 2 * pt(-abs(b / se), df)
  f <- as.formula(paste("y ~ z +", paste(cov, collapse = " + ")))
  err <- ""
  note <- ""
  used <- model   # used = model actually fitted (a failed ordinal fit falls back to binary)
  o <- tryCatch({
    if (model == "gaussian") {
      m <- svyglm(f, des)
      s <- summary(m)$coefficients
      c(s["z", 1], s["z", 2], tp(s["z", 1], s["z", 2]))
    } else if (model == "binomial") {
      m <- svyglm(f, des, family = quasibinomial())
      s <- summary(m)$coefficients
      c(s["z", 1], s["z", 2], tp(s["z", 1], s["z", 2]))
    } else if (model == "ordinal") {
      r <- tryCatch({
        m <- svyolr(f, des)
        b <- coef(m)[["z"]]
        se <- sqrt(vcov(m)["z", "z"])
        c(b, se, tp(b, se))
      }, error = function(e) NULL)
      if (is.null(r)) {
        des$variables$y <- as.numeric(as.integer(des$variables$y) > 1)
        note <- "svyolr singular -> refitted as binomial: any level above the lowest vs the lowest"
        used <- "binomial (ordinal fallback)"
        m <- svyglm(f, des, family = quasibinomial())
        s <- summary(m)$coefficients
        r <- c(s["z", 1], s["z", 2], tp(s["z", 1], s["z", 2]))
      }
      r
    } else {
      m <- svyglm(as.formula(paste("z ~ y +", paste(cov, collapse = " + "))), des)
      c(NA, NA, regTermTest(m, ~y)$p[1])
    }
  }, error = function(e) {
    err <<- conditionMessage(e)
    c(NA, NA, NA)
  })
  list(b = o[1], se = o[2], p = o[3], err = err, note = note, df = df, used = used)
}

R <- list()
YN <- list()
run <- function(key, y, model, keep, cov) {
  D <- DS[[key]]
  n <- sum(keep & D$rows)
  # phase-specific items
  if (n < 50) {
    return(list(fit = list(b = NA, se = NA, p = NA, err = if (n == 0) "not asked in this phase" else sprintf("n = %d < 50", n),
                           note = "", df = NA, used = model), n = n))
  }
  des <- D$des
  des$variables$y <- y[D$rows]
  ds <- subset(des, keep[D$rows])
  list(fit = fit(model, ds, cov), n = n)
}
for (i in seq_len(nrow(A))) {
  v <- A$variable[i]
  r <- A[i, ]
  pr <- prep(v, r)
  women <- r$health_class == "reproductive_exposure"
  keep <- !is.na(pr$y) & complete.cases(d[, COV]) & (!women | d$sex == "2")
  # fasting-subsample measures (WTPFSD6): glucose, insulin and C-peptide (fasting only in the morning session; not HbA1c)
  #   and variables measured only in the morning subsample (LDL)
  fast <- (grepl("DIABETES TESTING", r$section) && v != "GHP") || (sum(keep) >= 50 && mean(d$WTPFSD6[keep] > 0) > 0.95)
  cov <- if (women) setdiff(COV, "sex") else COV
  a1 <- run(if (fast) "fp1" else "p1", pr$y, pr$model, keep, cov)
  a2 <- run(if (fast) "fp2" else "p2", pr$y, pr$model, keep, cov)
  # items collected in one phase only (judged from valid n per phase) use that phase's design
  key <- if (a1$n >= 50 && a2$n < 50) (if (fast) "fp1" else "p1")
         else if (a2$n >= 50 && a1$n < 50) (if (fast) "fp2" else "p2")
         else (if (fast) "fall" else "all")
  single <- key != (if (fast) "fall" else "all")
  a <- run(key, pr$y, pr$model, keep, cov)
  ex <- pr$model %in% c("binomial", "ordinal")
  tr <- function(b) if (ex) exp(b) else b
  f <- a$fit
  f1 <- a1$fit
  f2 <- a2$fit
  R[[v]] <- data.frame(
    variable = v, phenotype = r$label, set = r$status, health_class = r$health_class, model = pr$model,
    n = a$n, n_phase1 = a1$n, n_phase2 = a2$n, levels = pr$levels,
    weight = if (single) sprintf("%s (single-phase item)", c(p1 = "WTPFEX1", p2 = "WTPFEX2", fp1 = "WTPFSD1", fp2 = "WTPFSD2")[[key]])
             else if (fast) "WTPFSD6 (phase: WTPFSD1/2)" else "WTPFEX6 (phase: WTPFEX1/2)",
    population = if (women) "women" else "all", df = f$df,
    measure = c(gaussian = "beta (SD of phenotype per SD KSAA)", binomial = "OR per SD KSAA", ordinal = "OR per SD KSAA (proportional odds)",
                `binomial (ordinal fallback)` = "OR per SD KSAA (any level above the lowest vs the lowest)", nominal = "joint test only")[[f$used]],
    model_used = f$used,
    # CI and P from the t distribution with the same df
    est = tr(f$b), lo = tr(f$b - qt(0.975, f$df) * f$se), hi = tr(f$b + qt(0.975, f$df) * f$se), p = f$p,
    est_phase1 = tr(f1$b), p_phase1 = f1$p, est_phase2 = tr(f2$b), p_phase2 = f2$p,
    # single-phase items cannot be replicated across phases
    replicated = if (single) NA else (!is.na(f1$b) && !is.na(f2$b) && sign(f1$b) == sign(f2$b) && f1$p < 0.05 && f2$p < 0.05),
    notes = paste(c(pr$notes, f$note,
        if (single) "single-phase item: main analysis uses that phase's design; no cross-phase replication" else NULL),
        collapse = "; "), error = f$err, section = r$section, subsection = r$subsection, stringsAsFactors = FALSE)
  YN[[v]] <- if (pr$model == "nominal") NULL else {
    z <- rep(NA_real_, nrow(d))
    z[keep] <- as.numeric(pr$y)[keep]
    z
  }
  if (i %% 50 == 0) cat(sprintf("  %d / %d\n", i, nrow(A)))
}
R <- do.call(rbind, R)
rownames(R) <- NULL
FASTV <- sort(R$variable[R$set == "main" & grepl("WTPFSD", R$weight)])
cat("Fasting-subsample variables (main set):", paste(FASTV, collapse = ", "), "\n")
# six glucose-tolerance/insulin measures + LDL; HbA1c (GHP) is not among them
stopifnot(setequal(FASTV, c("C1P", "C2P", "G1P", "G2P", "I1P", "I2P", "LCP")))

# ---------- Multiple testing: within each set ----------
liji <- function(M) {
  M <- M[, colSums(!is.na(M)) > 0, drop = FALSE]
  C <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  C[is.na(C)] <- 0
  diag(C) <- 1
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  ev[ev < 0] <- 0
  sum(as.numeric(ev >= 1) + (ev - floor(ev)))
}
for (s in unique(R$set)) {
  k <- R$set == s & !is.na(R$p)
  R$fdr[k] <- p.adjust(R$p[k], "BH")
  R$bonf[k] <- p.adjust(R$p[k], "bonferroni")
  vs <- R$variable[k & R$variable %in% names(YN)]
  M <- do.call(cbind, YN[vs])
  me <- liji(M)
  R$meff[k] <- me
  R$p_liji[k] <- pmin(1, R$p[k] * me)
}
# independent signals: main-analysis FDR < 5% hits clumped in order of P at |r| > 0.7 (pairwise-complete r)
k <- which(R$set == "main" & !is.na(R$fdr) & R$fdr < 0.05 & R$variable %in% names(YN))
k <- k[order(R$p[k])]
R$clump_lead <- NA_character_
leads <- character()
for (i in k) {
  v <- R$variable[i]
  hit <- NA
  for (l in leads) {
    r <- suppressWarnings(cor(YN[[v]], YN[[l]], use = "pairwise.complete.obs"))
    if (!is.na(r) && abs(r) > 0.7) { hit <- l; break }
  }
  if (is.na(hit)) { leads <- c(leads, v); R$clump_lead[i] <- v } else R$clump_lead[i] <- hit
}
saveRDS(YN, file.path(ew, "phewas_y.rds"))
R <- R[order(R$set != "main", R$p), ]
write.csv(R, file.path(ew, "phewas_results.csv"), row.names = FALSE)

sink(file.path(ew, "phewas_summary.txt"))
cat("PheWAS (audit table v1.2; KSAA -> phenotype; weighted; covariates age/sex/race/education/PIR)\n\n")
for (s in c("main", "supplementary", "exposure_module")) {
  k <- R$set == s
  cat(sprintf("[%s] phenotypes %d (failed fits %d); models: %s; FDR<5%% %d; Bonferroni %d; Li-Ji meff %.0f -> %d; replicated in both phases %d\n",
              s, sum(k), sum(k & is.na(R$p)), paste(names(table(R$model[k])), table(R$model[k]), collapse = " "),
              sum(R$fdr[k] < 0.05, na.rm = TRUE), sum(R$bonf[k] < 0.05, na.rm = TRUE), R$meff[k][1],
              sum(R$p_liji[k] < 0.05, na.rm = TRUE), sum(R$replicated[k], na.rm = TRUE)))
}
cat(sprintf("\nMain-analysis independent signals (FDR<5%% merged at |r|>0.7): %d; leads: %s\n",
    length(unique(na.omit(R$clump_lead))), paste(unique(na.omit(R$clump_lead)), collapse = " ")))
cat("\nMain analysis FDR<5% by health class:\n")
print(table(R$health_class[R$set == "main" & R$fdr < 0.05]))
cat("\nMain analysis top 30:\n")
k <- which(R$set == "main")[1:30]
for (i in k) cat(sprintf("  %-9s %-46s %-9s n=%-5d %s %s  P=%.1e FDR=%.1e  ph1 %.2f/%.1e ph2 %.2f/%.1e %s\n", R$variable[i],
    substr(R$phenotype[i], 1, 46), R$model[i], R$n[i],
  ifelse(is.na(R$est[i]), "   -  ", sprintf("%6.2f", R$est[i])), ifelse(is.na(R$lo[i]), "", sprintf("(%.2f-%.2f)", R$lo[i],
      R$hi[i])), R$p[i], R$fdr[i], R$est_phase1[i], R$p_phase1[i], R$est_phase2[i], R$p_phase2[i], ifelse(R$replicated[i],
      "REP", "")))
cat("\nFailed fits:\n")
print(R[is.na(R$p), c("variable", "phenotype", "model", "n", "error")], row.names = FALSE)
cat("\nFallback/recoding notes:\n")
print(table(R$notes[grepl("singular|merged", R$notes)]))
sink()
cat(readLines(file.path(ew, "phewas_summary.txt")), sep = "\n")
