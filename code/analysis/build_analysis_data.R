# Build the analysis data set (analysis_data.rds) and update work/oof_pheno.rds with Phenotypic Age recomputed in this cohort.
# Input: results/ksaa_oof.csv, results/predictions/m4_*.csv, results/work/*.csv and the NHANES III DuckDB
# Output: results/analysis/analysis_data.rds, results/work/oof_pheno.rds
.fresh_rng <- function() { RNGkind("default", "default", "default")   # each part starts from a fresh-session RNG state
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }

# ======== Part 1: analysis_data.rds (4,875 participants: outcomes, PhenoAge, covariates, K–L, BMD, function, per-knee readouts) ========
.fresh_rng(); local({
suppressMessages(library(duckdb))
res <- here::here("results"); sp <- file.path(res, "work"); out <- file.path(res, "analysis")
dir.create(out, showWarnings = FALSE)
q <- function(x){ x <- suppressWarnings(as.numeric(x)); for(k in 1:6){ s8 <- as.numeric(strrep("8", k)); lo <- x[!is.na(x) & x < s8]
  if(any(x == s8, na.rm = TRUE) && length(lo) && s8 >= max(lo) * 1.5) x[x == s8] <- NA }; x }   # all-8 code is treated as missing only when far above the other values
na_at <- function(x, codes) { x <- as.numeric(x); x[x %in% codes] <- NA; x }                      # missing codes given explicitly
yn <- function(x) ifelse(x == 1, 1L, ifelse(x == 2, 0L, NA_integer_))                               # 1 = yes, 2 = no, 8/9 = missing
prev <- function(v) sprintf("%-10s n %4d  missing %4d  positive %5.1f%%", v, sum(!is.na(d[[v]])), sum(is.na(d[[v]])), 100 * mean(d[[v]], na.rm = TRUE))

# ---------------- Cohort and outcomes ----------------
d <- read.csv(file.path(res, "ksaa_oof.csv"))
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN); stopifnot(!anyNA(i), all(d$age == sl$age[i]))
for(v in c("sex", "race", "phase", "pir", "bmi", "fnbmd", "WTPFEX6", "WTPFEX2", "SDPPSU6", "SDPSTRA6", "SDPPSU2", "SDPSTRA2")) d[[v]] <- sl[[v]][i]
d$t <- sl$PERMTH_EXM[i]; d$e <- sl$MORTSTAT[i]; d$ucod <- sl$UCOD_LEADING[i]; d$gap <- d$gap_corr   # time zero = MEC examination (PERMTH_EXM)
d$e_cvd <- as.integer(d$e == 1 & d$ucod %in% c(1, 5)); d$e_cancer <- as.integer(d$e == 1 & d$ucod %in% 2)   # UCOD_LEADING: 1 heart disease, 5 cerebrovascular, 2 cancer; the 5 deaths without a cause count as neither
cat(sprintf("Cohort n %d; deaths %d (no cause %d); CVD deaths %d; cancer deaths %d\n", nrow(d), sum(d$e), sum(d$e == 1 & is.na(d$ucod)), sum(d$e_cvd), sum(d$e_cancer)))
cat(sprintf("surv_labels bmi range %s, fnbmd range %s (confirm cleaned)\n", paste(range(d$bmi, na.rm = TRUE), collapse = "–"), paste(round(range(d$fnbmd, na.rm = TRUE), 3), collapse = "–")))

# ---------------- NHANES III source variables ----------------
con <- dbConnect(duckdb::duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
dbWriteTable(con, "ids", data.frame(SEQN = d$SEQN), temporary = TRUE)
get <- function(tb, v) { x <- dbGetQuery(con, sprintf("SELECT i.SEQN, %s FROM ids i LEFT JOIN %s t USING(SEQN)", paste0("t.", v, collapse = ","), tb)); x[match(d$SEQN, x$SEQN), ] }
PA <- list(c("HAT1S", "HAT1MET", NA), c("HAT3S", "HAT2MET", "HAT2"), c("HAT5S", "HAT4MET", "HAT4"), c("HAT7S", "HAT6MET", "HAT6"), c("HAT9S", "HAT8MET", "HAT8"),
           c("HAT11S", "HAT10MET", "HAT10"), c("HAT13S", "HAT12MET", "HAT12"), c("HAT15S", "HAT14MET", "HAT14"), c("HAT17S", "HAT16MET", "HAT16"),
           c("HAT20S", "HAT19MET", "HAT18"), c("HAT22S", "HAT21MET", "HAT18"), c("HAT24S", "HAT23MET", "HAT18"), c("HAT26S", "HAT25MET", "HAT18"))
AD <- get("nhanes3_adult", unique(c("HFA8R", "HAR1", "HAR3", "HAD1", "HAD6", "HAE2", "HAE5A", "HAC1A", "HAC1C", "HAC1D", "HAC1O", "HAF10", na.omit(unlist(PA)))))
EX <- get("nhanes3_exam", c("MAPE1", "MAPE2", "PEPMNK1R", "PEPMNK5R", "BMPWAIST", "PFPTTAND"))
LB <- get("nhanes3_lab", c("AMP", "CEP", "SGP", "CRP", "LMPPCNT", "MVPSI", "RWP", "APPSI", "WCP", "TCP", "HDP", "TGP", "GHP", "FBP", "PHPFAST", "G1P", "I1P", "PHPINSU"))
dbDisconnect(con, shutdown = TRUE)

# ---------------- PhenoAge ----------------
lim <- c(AMP = 8, CEP = 25, SGP = 800, CRP = 50, LMPPCNT = 100, MVPSI = 150, RWP = 50, APPSI = 2000, WCP = 300)
for(v in names(lim)) { x <- LB[[v]]; x[x >= lim[v]] <- NA; d[[v]] <- x }
pheno <- function(x, use_crp = TRUE){
  xb <- -19.9067 - 0.0336 * (x$AMP * 10) + 0.0095 * (x$CEP * 88.4) + 0.1953 * (x$SGP * 0.0555) - 0.0120 * x$LMPPCNT + 0.0268 * x$MVPSI +
    0.3306 * x$RWP + 0.00188 * x$APPSI + 0.0554 * x$WCP + 0.0804 * x$age
  if(use_crp) xb <- xb + 0.0954 * log(x$CRP)
  g <- 0.0076927; k <- (exp(120 * g) - 1) / g
  141.50225 + (log(0.00553 * k) + xb) / 0.09165   # log(1 - M) = -exp(xb) * k, written this way to avoid overflow
}
d$phenoage <- pheno(d); d$pa_accel <- resid(lm(phenoage ~ age, d, na.action = na.exclude))
# PhenoAge acceleration is winsorized at the 1st/99th percentiles (long right tail); pa_accel_raw is kept for descriptive tables
d$pa_accel_raw <- d$pa_accel
PA_W <- quantile(d$pa_accel_raw, c(.01, .99), na.rm = TRUE)
d$pa_accel <- pmin(pmax(d$pa_accel_raw, PA_W[1]), PA_W[2])
cat(sprintf("\nPhenoAge acceleration winsorized: 1%% / 99%% percentiles %.2f / %.2f years; %d set to the lower and %d to the upper limit; SD %.2f -> %.2f\n", PA_W[1], PA_W[2],
    sum(d$pa_accel_raw < PA_W[1], na.rm = TRUE), sum(d$pa_accel_raw > PA_W[2], na.rm = TRUE), sd(d$pa_accel_raw, na.rm = TRUE), sd(d$pa_accel, na.rm = TRUE)))
old <- readRDS(file.path(sp, "oof_pheno.rds")); j <- match(d$SEQN, old$SEQN)
cat(sprintf("\nPhenoAge computable %d / %d; max difference vs old oof_pheno.rds in shared participants %.2e; newly computable %d (SEQN %s)\n", sum(!is.na(d$phenoage)), nrow(d),
    max(abs(d$phenoage - old$phenoage[j]), na.rm = TRUE), sum(!is.na(d$phenoage) & is.na(j)), paste(d$SEQN[!is.na(d$phenoage) & is.na(j)], collapse = ", ")))

# ---------------- Sociodemographics and lifestyle ----------------
ed <- na_at(AD$HFA8R, c(88, 99)); d$edu3 <- cut(ed, c(-Inf, 11, 12, Inf), labels = c("<12", "12", ">12"))
d$smk3 <- factor(ifelse(AD$HAR1 == 2, "never", ifelse(AD$HAR1 == 1 & AD$HAR3 == 2, "former", ifelse(AD$HAR1 == 1 & AD$HAR3 == 1, "current", NA))), levels = c("never", "former", "current"))
# Alcohol: <12 drinks in lifetime = never; >=12 in lifetime but <12 in the past year = former/light (includes current light drinkers); >=12 in the past year = current
d$alc3 <- factor(ifelse(EX$MAPE1 == 2, "never", ifelse(EX$MAPE1 == 1 & EX$MAPE2 == 2, "former_light", ifelse(EX$MAPE1 == 1 & EX$MAPE2 == 1, "current", NA))), levels = c("never", "former_light", "current"))
# Leisure-time physical activity: sessions in the past month x MET; moderate 3-5.9, vigorous >= 6.
# Active = moderate >= 20/month (~5/week) or vigorous >= 12/month (~3/week); inactive = none
bad <- rep(FALSE, nrow(d)); mod <- vig <- tot <- mets <- rep(0, nrow(d)); nomet <- 0
for(p in PA){ f <- as.numeric(AD[[p[1]]]); m <- as.numeric(AD[[p[2]]]); gate <- if(is.na(p[3])) rep(1, nrow(d)) else as.numeric(AD[[p[3]]])
  bad <- bad | (!is.na(f) & f >= 7777) | (!is.na(gate) & gate %in% c(8, 9)); f[is.na(f) | f >= 7777] <- 0
  bad <- bad | (f > 0 & !is.na(m) & m >= 88); m[!is.na(m) & m >= 88] <- NA   # HAT19MET = 8888 (3 people): intensity unknown, set to missing
  nomet <- nomet + sum(f > 0 & is.na(m)); m[is.na(m)] <- 3
  tot <- tot + f; mets <- mets + f * m; mod <- mod + f * (m >= 3 & m < 6); vig <- vig + f * (m >= 6) }
d$pa_sessions <- ifelse(bad, NA, tot); d$pa_metsess <- ifelse(bad, NA, mets)
d$pa3 <- factor(ifelse(bad, NA, ifelse(tot == 0, "inactive", ifelse(mod >= 20 | vig >= 12, "active", "insufficient"))), levels = c("active", "insufficient", "inactive"))
cat(sprintf("Physical activity: sentinel frequency or gate %d; frequency without MET (set to 3) %d activity-person; monthly sessions P50/P90/max %s\n", sum(bad), nomet, paste(quantile(d$pa_sessions, c(.5, .9, 1), na.rm = TRUE), collapse = "/")))

# ---------------- Cardiometabolic factors and medical history ----------------
d$sbp <- na_at(EX$PEPMNK1R, 888); d$dbp <- na_at(EX$PEPMNK5R, 888); d$waist <- na_at(EX$BMPWAIST, 88888); d$fast_h <- na_at(LB$PHPFAST, 88888)
d$tc <- na_at(LB$TCP, 888); d$hdl <- na_at(LB$HDP, 888); d$tg <- na_at(LB$TGP, 8888); d$hba1c <- na_at(LB$GHP, 8888); d$fib <- na_at(LB$FBP, 8888)
or1 <- function(...) { m <- cbind(...); ifelse(rowSums(m == 1, na.rm = TRUE) > 0, 1L, ifelse(rowSums(is.na(m)) == 0, 0L, NA_integer_)) }
d$htn <- or1(yn(AD$HAE2), as.integer(d$sbp >= 140), as.integer(d$dbp >= 90))                       # self-reported physician diagnosis or measured >=140/90
d$dm <- or1(yn(AD$HAD1), as.integer(d$hba1c >= 6.5))                                               # self-report or HbA1c >=6.5%
d$dm[is.na(d$dm) & AD$HAD1 == 2 & is.na(d$hba1c)] <- 0L                                            # self-reported no with HbA1c missing -> 0
bpmed <- ifelse(AD$HAE2 == 2 & is.na(AD$HAE5A), 0L, yn(AD$HAE5A))                                   # HAE5A is asked only of those reporting hypertension
d$dyslip <- or1(as.integer(d$tc >= 240), as.integer(d$hdl < 40))
d$cvd_hist <- or1(yn(AD$HAF10), yn(AD$HAC1C), yn(AD$HAC1D))                                        # myocardial infarction, heart failure, stroke
d$cancer_hist <- yn(AD$HAC1O)                                                                       # non-skin cancer ("other cancer")
d$arth <- yn(AD$HAC1A); d$obese <- as.integer(d$bmi >= 30)
# eGFR: NHANES III creatinine calibrated to the standardized method (0.960 x Scr - 0.184; Selvin 2007), then CKD-EPI 2021 (no race term)
scr <- 0.960 * d$CEP - 0.184; fem <- d$sex == 2; kap <- ifelse(fem, 0.7, 0.9); alp <- ifelse(fem, -0.241, -0.302)
d$egfr <- 142 * pmin(scr / kap, 1)^alp * pmax(scr / kap, 1)^-1.200 * 0.9938^d$age * ifelse(fem, 1.012, 1)
# Metabolic syndrome (>= 3 of 5 ATP III criteria). Fasted >= 8 h: TG >= 150 and glucose >= 100; otherwise TG >= 175 and HbA1c >= 5.7.
# Diabetes and BP medication count toward their criteria
fasting <- !is.na(d$fast_h) & d$fast_h >= 8
cr <- cbind(waist = as.integer(ifelse(fem, d$waist > 88, d$waist > 102)),
            tg    = as.integer(ifelse(fasting, d$tg >= 150, d$tg >= 175)),
            hdl   = as.integer(ifelse(fem, d$hdl < 50, d$hdl < 40)),
            bp    = or1(as.integer(d$sbp >= 130), as.integer(d$dbp >= 85), bpmed),
            glu   = or1(d$dm, as.integer(ifelse(fasting, d$SGP >= 100, d$hba1c >= 5.7))))
nmet <- rowSums(cr == 1, na.rm = TRUE); nmis <- rowSums(is.na(cr)); d$mets <- ifelse(nmet >= 3, 1L, ifelse(nmet + nmis < 3, 0L, NA_integer_))
# HOMA-IR = G1P (mg/dL) x I1P (uU/mL) / 405, in those fasted >= 8 h and not using insulin; missing codes G1P 88888, I1P 888888
d$fpg <- na_at(LB$G1P, 88888); d$insulin <- na_at(LB$I1P, 888888); ins <- (yn(AD$HAD6) %in% 1L) | (LB$PHPINSU %in% 1)
d$lhoma <- ifelse(fasting & !ins, log(d$fpg * d$insulin / 405), NA); d$linsulin <- ifelse(fasting & !ins, log(d$insulin), NA)
cat(sprintf("HOMA-IR: fasting >=8 h, no insulin and both values valid %d (insulin users %d); PHPINSU values %s; median HOMA-IR %.2f\n", sum(!is.na(d$lhoma)), sum(ins),
    paste(names(table(LB$PHPINSU, useNA = "ifany")), table(LB$PHPINSU, useNA = "ifany"), sep = ":", collapse = " "), exp(median(d$lhoma, na.rm = TRUE))))
cat(sprintf("Fasting >=8 h: %d (%.0f%%); metabolic syndrome criteria met %s\n", sum(fasting), 100 * mean(fasting), paste(sprintf("%s %.0f%%", colnames(cr), 100 * colMeans(cr, na.rm = TRUE)), collapse = "  ")))

# ---------------- Musculoskeletal and inflammation ----------------
kl <- read.csv(file.path(res, "work", "kl_per_knee.csv")); cat(sprintf("K–L values: %s\n", paste(names(table(kl$kl)), table(kl$kl), sep = ":", collapse = "  ")))
klp <- aggregate(kl ~ seqn, kl[kl$kl %in% 0:4, ], max); d$klmax <- klp$kl[match(d$SEQN, klp$seqn)]
d$oa <- as.integer(d$klmax >= 2); d$severe_oa <- as.integer(d$klmax >= 3)
KP <- read.csv(file.path(sp, "knee_pain.csv")); d$pain <- yn(q(KP$HAG22[match(d$SEQN, KP$SEQN)]))
DX <- read.csv(file.path(sp, "dxa_all.csv")); j <- match(d$SEQN, DX$SEQN); d$fn_bmd <- q(DX$BDPFNBMD)[j]; d$tot_bmd <- (q(DX$BDPTOBMC) / q(DX$BDPTOARE))[j]
P <- read.csv(file.path(sp, "physfun.csv")); j <- match(d$SEQN, P$SEQN)
d$gait <- 2.4384 / pmin(q(P$PFPTWLKA), q(P$PFPTWLKB), na.rm = TRUE)[j]; d$gait[!is.finite(d$gait)] <- NA
d$tandem <- na_at(EX$PFPTTAND, 8888)
d$lcrp <- log(d$CRP); d$wbc <- d$WCP; d$glucose <- d$SGP
# Fracture, falls, functional limitation, poverty (for the baseline table by quartile)
qn <- function(x){ x <- suppressWarnings(as.numeric(x)); for(k in 1:6) for(dg in c("8", "9")){ s <- as.numeric(strrep(dg, k)); lo <- x[!is.na(x) & x < s]
  if(any(x == s, na.rm = TRUE) && length(lo) && s >= max(lo) * 1.5) x[x == s] <- NA }; x }
FR <- read.csv(file.path(sp, "fracture.csv")); jf <- match(d$SEQN, FR$SEQN); jp <- match(d$SEQN, P$SEQN)
d$fx_hip <- yn(FR$HAG5A[jf]); fxm <- cbind(d$fx_hip, yn(FR$HAG5B[jf]), yn(FR$HAG5C[jf]))
d$fx_any <- as.integer(rowSums(fxm, na.rm = TRUE) > 0); d$fx_any[rowSums(!is.na(fxm)) == 0] <- NA
d$fall <- as.integer(qn(FR$HAG14)[jf] > 0)
hh <- sapply(c("HAH1", "HAH2", "HAH8", "HAH9"), function(v) { x <- P[[v]][jp]; as.integer(ifelse(x %in% c(8, 9), NA, x) >= 2) })
d$anylimit <- as.integer(rowSums(hh, na.rm = TRUE) > 0); d$anylimit[rowSums(!is.na(hh)) == 0] <- NA
d$poor <- as.integer(d$pir < 1)

# ---------------- Left/right knee readouts (local heads of the bilateral model; same fold weights and bias correction) ----------------
# The local heads also see the other knee through cross-attention, so gap_L / gap_R are not independent single-knee estimates
A5 <- c(resnet18 = "resnet18", resnet50 = "resnet50", efficientnetb0 = "efficientnetb0", convnexttiny = "convnexttiny", densenet121 = "densenet121")
FM <- sapply(A5, function(a) read.csv(file.path(res, "predictions", sprintf("m4_%s_folds.csv", a)))$inner_val_mae); W <- apply(FM, 1, function(m) (1/m) / sum(1/m))
for(sd_ in c("L", "R")){
  M <- sapply(A5, function(a){ x <- read.csv(file.path(res, "predictions", sprintf("m4_%s_oof.csv", a))); x[[paste0("pred_", sd_)]][match(d$SEQN, x$SEQN)] }); stopifnot(!anyNA(M))
  raw <- rep(NA_real_, nrow(d)); for(k in 0:4){ i <- d$fold == k; raw[i] <- as.numeric(M[i, , drop = FALSE] %*% W[, k + 1]) }
  g <- rep(NA_real_, nrow(d)); for(k in 0:4){ te <- d$fold == k; f <- lm(I(raw - age) ~ age, data.frame(raw = raw[!te], age = d$age[!te])); g[te] <- raw[te] - d$age[te] - predict(f, data.frame(age = d$age[te])) }
  d[[paste0("gap_", sd_)]] <- g }
# check: the mean of the two knee readouts should be close to the primary readout
M <- sapply(A5, function(a){ x <- read.csv(file.path(res, "predictions", sprintf("m4_%s_oof.csv", a))); x$pred_age[match(d$SEQN, x$SEQN)] })
raw <- rep(NA_real_, nrow(d)); for(k in 0:4){ i <- d$fold == k; raw[i] <- as.numeric(M[i, , drop = FALSE] %*% W[, k + 1]) }
cat(sprintf("\nRebuilt primary raw prediction vs ksaa_oof.csv max difference %.2e; correlation of left/right knee gap with primary gap %.3f / %.3f; left-right correlation %.3f\n",
    max(abs(raw - d$pred_raw)), cor(d$gap_L, d$gap), cor(d$gap_R, d$gap), cor(d$gap_L, d$gap_R)))

# ---------------- Checks: missingness, prevalence, direction with age ----------------
cat("\nBinary variables:\n"); for(v in c("htn", "dm", "dyslip", "cvd_hist", "cancer_hist", "arth", "obese", "mets", "oa", "severe_oa", "pain")) cat("  ", prev(v), "\n")
cat("\nCategorical variables:\n"); for(v in c("edu3", "smk3", "alc3", "pa3")) cat(sprintf("  %-5s %s  missing %d\n", v, paste(names(table(d[[v]])), table(d[[v]]), sep = ":", collapse = "  "), sum(is.na(d[[v]]))))
cat("\nContinuous variables (n, mean, SD, min, max):\n")
for(v in c("sbp", "dbp", "tc", "hdl", "tg", "hba1c", "egfr", "waist", "fib", "lcrp", "wbc", "glucose", "fn_bmd", "tot_bmd", "gait", "tandem", "pa_metsess", "phenoage", "pa_accel_raw", "pa_accel"))
  cat(sprintf("  %-10s n %4d  %8.2f  %7.2f  %8.2f  %8.2f\n", v, sum(!is.na(d[[v]])), mean(d[[v]], na.rm = TRUE), sd(d[[v]], na.rm = TRUE), min(d[[v]], na.rm = TRUE), max(d[[v]], na.rm = TRUE)))
cat("\nDirection with age (OR per 10 years; should be plausible):\n")
for(v in c("htn", "dm", "cvd_hist", "cancer_hist", "oa", "mets")) cat(sprintf("  %-12s %.2f\n", v, exp(10 * coef(glm(d[[v]] ~ d$age, family = binomial))[2])))
cat(sprintf("  inactive vs other %.2f; current drinking vs other %.2f; current smoking vs other %.2f\n", exp(10 * coef(glm(I(pa3 == "inactive") ~ age, binomial, d))[2]),
    exp(10 * coef(glm(I(alc3 == "current") ~ age, binomial, d))[2]), exp(10 * coef(glm(I(smk3 == "current") ~ age, binomial, d))[2])))

saveRDS(d, file.path(out, "analysis_data.rds")); cat(sprintf("\n-> %s (%d × %d)\n", file.path(out, "analysis_data.rds"), nrow(d), ncol(d)))
})

# ======== Part 2: replace work/oof_pheno.rds with the PhenoAge values from part 1 (original kept as oof_pheno_v1cohort4872.rds) ========
.fresh_rng(); local({
# Downstream scripts use SEQN, phenoage and pa_accel (the first column containing "accel"); other columns are carried over by SEQN.
res <- here::here("results"); sp <- file.path(res, "work")
bk <- file.path(sp, "oof_pheno_v1cohort4872.rds"); old <- readRDS(file.path(sp, "oof_pheno.rds"))
if (!file.exists(bk)) { stopifnot(nrow(old) == 4872); saveRDS(old, bk); cat("backup ->", bk, "\n") } else old <- readRDS(bk)
v <- readRDS(file.path(res, "analysis", "analysis_data.rds"))
# pa_accel: winsorized (analyses); pa_accel_raw: raw (descriptive tables)
new <- data.frame(SEQN = v$SEQN, phenoage = v$phenoage, pa_accel = v$pa_accel, zpa = as.numeric(scale(v$pa_accel)), pa_accel_raw = v$pa_accel_raw)
keep <- setdiff(names(old), c("phenoage", "pa_accel", "zpa", "phenoage_nc", "znc", "pa_accel_raw"))
new <- merge(new, old[, keep, drop = FALSE], by = "SEQN", all.x = TRUE)
stopifnot(nrow(new) == 4875, grep("accel", names(new), value = TRUE)[1] == "pa_accel")
i <- match(old$SEQN, new$SEQN)
cat(sprintf("New oof_pheno.rds: %d participants, PhenoAge non-missing %d (old %d); in shared participants max difference phenoage %.1e, pa_accel %.3f (residuals re-estimated in the new sample)\n",
    nrow(new), sum(!is.na(new$phenoage)), sum(!is.na(old$phenoage)), max(abs(new$phenoage[i] - old$phenoage), na.rm = TRUE), max(abs(new$pa_accel[i] - old$pa_accel), na.rm = TRUE)))
saveRDS(new, file.path(sp, "oof_pheno.rds")); cat("->", file.path(sp, "oof_pheno.rds"), "\n")
})
