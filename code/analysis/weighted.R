# Survey-weighted primary analyses (regression models); standardized risks and discrimination are in weighted_boot.R.
# Output (results/analysis/weighted/*.csv): table2, tableS4, tableS5, tableS7, tableS10 (+ _gait_extra), figure3a, figure3_cells,
# figure3_within, figureS3, model1_full_cohort
source(here::here("code", "analysis", "weighted_common.R"))
suppressMessages({library(mice); library(mitools)})

# ---------------- Table 2: three outcomes × Models 1–3, per SD and by quartile, P for trend ----------------
cat("=========== Table 2 (weighted) ===========\n")
T2 <- NULL
for (o in names(EV)) for (mn in names(MODS)) { des <- dom6(F$cc); ev <- EV[[o]]
  a <- wcox(des, paste("z +", MODS[[mn]]), ev = ev)
  mq <- svycoxph(as.formula(sprintf("Surv(t, %s) ~ q + %s", ev, MODS[[mn]])), design = des)
  qq <- t(sapply(paste0("q", QL[-1]), function(k) ci3(coef(mq)[[k]], sqrt(vcov(mq)[k, k]))))
  tr <- wcox(des, paste("qmed +", MODS[[mn]]), term = "qmed", ev = ev)
  T2 <- rbind(T2, data.frame(outcome = o, model = mn, n = a[["n"]], events = a[["events"]], hr_sd = a[["hr"]], lo_sd = a[["lo"]], hi_sd = a[["hi"]], p_sd = a[["p"]],
    hr_q2 = qq[1, 1], lo_q2 = qq[1, 2], hi_q2 = qq[1, 3], hr_q3 = qq[2, 1], lo_q3 = qq[2, 2], hi_q3 = qq[2, 3], hr_q4 = qq[3, 1], lo_q4 = qq[3, 2], hi_q4 = qq[3, 3], p_trend = tr[["p"]])) }
write.csv(T2, file.path(W, "table2.csv"), row.names = FALSE)
for (i in seq_len(nrow(T2))) with(T2[i, ], cat(sprintf("  %-8s %s  SD %s  Q2 %s  Q3 %s  Q4 %s  P trend %.3g\n", outcome, model, fci(hr_sd, lo_sd, hi_sd), fci(hr_q2, lo_q2, hi_q2), fci(hr_q3, lo_q3, hi_q3), fci(hr_q4, lo_q4, hi_q4), p_trend)))

# ---------------- Table S4: Fine–Gray (weight = survey weight x fgwt after finegray expansion) ----------------
cat("\n=========== Table S4 Fine–Gray (weighted) ===========\n")
C <- F[F$cc, ]; C$id <- seq_len(nrow(C))
FG <- NULL
for (et in c("cvd", "cancer")) {
  fg <- finegray(Surv(t, cause) ~ ., data = C[, c("id", "t", "cause", "z", "q", "qmed", "W6", "P6", "S6", M3V)], etype = et)
  fg$wfg <- fg$W6 * fg$fgwt
  des <- svydesign(id = ~P6, strata = ~S6, weights = ~wfg, nest = TRUE, data = fg)
  for (mn in names(MODS)) {
    fit <- function(rhs, term) { m <- svycoxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", rhs)), design = des); b <- coef(m)[[term]]; se <- sqrt(vcov(m)[term, term]); c(ci3(b, se), p = 2 * pnorm(-abs(b / se))) }
    a <- fit(paste("z +", MODS[[mn]]), "z"); mq <- svycoxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~ q +", MODS[[mn]])), design = des)
    qq <- t(sapply(paste0("q", QL[-1]), function(k) ci3(coef(mq)[[k]], sqrt(vcov(mq)[k, k])))); tr <- fit(paste("qmed +", MODS[[mn]]), "qmed")
    FG <- rbind(FG, data.frame(outcome = et, model = mn, hr_sd = a[["hr"]], lo_sd = a[["lo"]], hi_sd = a[["hi"]], hr_q2 = qq[1, 1], lo_q2 = qq[1, 2], hi_q2 = qq[1, 3],
      hr_q3 = qq[2, 1], lo_q3 = qq[2, 2], hi_q3 = qq[2, 3], hr_q4 = qq[3, 1], lo_q4 = qq[3, 2], hi_q4 = qq[3, 3], p_trend = tr[["p"]])) } }
write.csv(FG, file.path(W, "tableS4.csv"), row.names = FALSE)
for (i in seq_len(nrow(FG))) with(FG[i, ], cat(sprintf("  %-7s %s  SD %s  Q4 %s  P trend %.3g\n", outcome, model, fci(hr_sd, lo_sd, hi_sd), fci(hr_q4, lo_q4, hi_q4), p_trend)))

# ---------------- Table S5: sensitivity analyses (weighted Model 3; phase 2 rows use the phase 2 design) ----------------
cat("\n=========== Table S5 (weighted) ===========\n")
spw <- SPW
mt <- read.csv(file.path(spw, "metal_scan.csv")); mt$SEQN <- as.integer(sub("K.*", "", mt$roi)); art <- unique(mt$SEQN[mt$frac_top > 1e-3])
resc <- unique(as.integer(sub("K.*", "", c(read.csv(file.path(spw, "boxes_rot.csv"))$stem, read.csv(file.path(spw, "boxes_lowconf.csv"))$stem))))
X <- read.csv(file.path(spw, "xray_components.csv")); q8 <- function(x) ifelse(x %in% c(8, 9, 88, 99, 888, 999), NA, x)
jr <- suppressWarnings(pmax(q8(X$XRPJRR), q8(X$XRPJRL), na.rm = TRUE)); jr[!is.finite(jr)] <- NA; F$jrep <- as.integer(jr == 1)[match(F$SEQN, X$SEQN)]
F$ab5 <- cut(F$age, c(59, 64, 69, 74, 79, 84, 120))
dd <- d; dd$pooled <- resid(lm(pred_raw ~ age, dd)); fmn <- tapply(dd$gap, dd$fold, mean); dd$fcent <- dd$gap - fmn[as.character(dd$fold)]
SDX <- sapply(c("gap_raw", "pooled", "fcent"), function(v) sd(dd[[v]]))
for (v in c("pooled", "fcent")) F[[v]] <- dd[[v]][match(F$SEQN, dd$SEQN)]
D6 <- svydesign(id = ~P6, strata = ~S6, weights = ~W6, nest = TRUE, data = F)        # rebuilt with the new columns
F2 <- F[F$phase_all == 2 & !is.na(F$W2) & F$W2 > 0, ]; D2 <- svydesign(id = ~P2, strata = ~S2, weights = ~W2, nest = TRUE, data = F2)
rhs3 <- function(sub, age_term = "age") { keep <- COV[!sapply(M3V, function(v) length(unique(F[sub, v])) < 2)]; keep[keep == "age"] <- age_term; paste(keep, collapse = " + ") }
wrow <- function(lab, sub, des_fun = dom6, x = "z", k = 1, age_term = "age") { a <- wcox(des_fun(sub), paste(x, "+", rhs3(sub, age_term)), term = x, k = k)
  data.frame(label = lab, n = a[["n"]], deaths = a[["events"]], hr = a[["hr"]], lo = a[["lo"]], hi = a[["hi"]]) }
cc <- F$cc; p2 <- cc & F$phase_all == 2 & !is.na(F$klmax)
S5 <- rbind(
  wrow("Complete-case sample", cc),
  wrow("Excluding deaths within 2 years", cc & !(F$e == 1 & F$t <= 24)),
  wrow("Excluding baseline cardiovascular disease", cc & F$cvd_hist == 0),
  wrow("Excluding baseline cancer", cc & F$cancer_hist == 0),
  wrow("Excluding baseline cardiovascular disease or cancer", cc & F$cvd_hist == 0 & F$cancer_hist == 0),
  wrow("Excluding films with bright artefacts", cc & !(F$SEQN %in% art)),
  wrow("Excluding rescued detections", cc & !(F$SEQN %in% resc)),
  wrow("Phase 2 with radiographic reading", p2, dom2),
  wrow("Excluding severe OA (K–L 3–4)", p2 & F$klmax <= 2, dom2),
  wrow("Excluding radiographic OA (K–L ≥2)", p2 & F$klmax <= 1, dom2),
  wrow("Excluding joint replacement", p2 & !is.na(F$jrep) & F$jrep == 0, dom2),
  wrow("Age modelled as a natural spline (5 df)", cc, age_term = "ns(age, df = 5)"),
  wrow("Stratified by 5-year age band", cc, age_term = "age + strata(ab5)"),
  wrow("Uncorrected age gap (predicted minus chronological age)", cc, x = "gap_raw", k = SDX[["gap_raw"]]),
  wrow("KSAA with the bias correction fitted in the whole cohort", cc, x = "pooled", k = SDX[["pooled"]]),
  wrow("KSAA centred within each outer fold", cc, x = "fcent", k = SDX[["fcent"]]))
u <- coxph(as.formula(paste("Surv(t, e) ~ z +", MODS[["m3"]])), F[cc, ]); ub <- coef(u)[["z"]]; use <- sqrt(vcov(u)["z", "z"])
S5 <- rbind(S5, data.frame(label = "Unweighted", n = u$n, deaths = u$nevent, t(ci3(ub, use))))
# Multiple imputation (mice, 20 imputations, Nelson–Aalen cumulative hazard as a predictor; design built on the cohort in each data set)
IV <- c("edu3", "smk3", "alc3", "pa3", "bmi", "htn", "dm", "tc", "hdl", "egfr", "cvd_hist", "cancer_hist")
di <- d[, c("t", "e", "e_cvd", "e_cancer", "gap", "z", "q", "age", "sex", "race", IV, "WTPFEX6", "SDPPSU6", "SDPSTRA6")]
for (v in c("sex", "race", "htn", "dm", "cvd_hist", "cancer_hist")) di[[v]] <- factor(di[[v]])
di$H <- nelsonaalen(di, t, e)
meth <- make.method(di); pred <- make.predictorMatrix(di); pred[, c("t", "z", "q", "WTPFEX6", "SDPPSU6", "SDPSTRA6")] <- 0
imp <- mice(di, m = 20, maxit = 10, method = meth, predictorMatrix = pred, printFlag = FALSE, seed = 20260914)
M3f <- gsub("factor\\((sex|race)\\)", "\\1", MODS[["m3"]])
fits <- lapply(1:20, function(i) { x <- complete(imp, i); svycoxph(as.formula(paste("Surv(t, e) ~ z +", M3f)), design = svydesign(id = ~SDPPSU6, strata = ~SDPSTRA6, weights = ~WTPFEX6, nest = TRUE, data = x)) })
mi <- summary(MIcombine(fits), alpha = 0.05); mz <- mi["z", ]
S5 <- rbind(S5, data.frame(label = "Multiple imputation, Model 3 covariates", n = nrow(d), deaths = sum(d$e), hr = exp(mz$results), lo = exp(mz$`(lower`), hi = exp(mz$`upper)`)))
hr <- S5$hr[1]; lo <- S5$lo[1]
ev_rr <- function(h) { rr <- (1 - 0.5^sqrt(h)) / (1 - 0.5^sqrt(1 / h)); rr + sqrt(rr * (rr - 1)) }
S5$evalue <- NA; S5$evalue[1] <- ev_rr(hr); S5$evalue_ci <- NA; S5$evalue_ci[1] <- ev_rr(lo)
write.csv(S5, file.path(W, "tableS5.csv"), row.names = FALSE)
for (i in seq_len(nrow(S5))) with(S5[i, ], cat(sprintf("  %-58s %5d / %5d  %s\n", label, deaths, n, fci(hr, lo, hi))))
cat(sprintf("  E-value %.2f (confidence limit %.2f)\n", S5$evalue[1], S5$evalue_ci[1]))

# ---------------- Table S7: KSAA with K–L / BMD in Model 3 (K–L uses the phase 2 design) ----------------
cat("\n=========== Table S7 (weighted) ===========\n")
J <- NULL
for (v in c("klmax", "fn_bmd", "tot_bmd")) { sub <- cc & !is.na(F[[v]]); df <- if (v == "klmax") dom2 else dom6; if (v == "klmax") sub <- sub & F$phase_all == 2
  F$m <- (F[[v]] - mean(F[[v]][sub])) / sd(F[[v]][sub]); D6 <- svydesign(id = ~P6, strata = ~S6, weights = ~W6, nest = TRUE, data = F)
  F2 <- F[F$phase_all == 2 & !is.na(F$W2) & F$W2 > 0, ]; D2 <- svydesign(id = ~P2, strata = ~S2, weights = ~W2, nest = TRUE, data = F2)
  for (o in c("allcause", "cvd")) { a <- wcox(df(sub), paste("z +", MODS[["m3"]]), ev = EV[[o]]); j <- svycoxph(as.formula(sprintf("Surv(t, %s) ~ z + m + %s", EV[[o]], MODS[["m3"]])), design = df(sub))
    kz <- ci3(coef(j)[["z"]], sqrt(vcov(j)["z", "z"])); km <- ci3(coef(j)[["m"]], sqrt(vcov(j)["m", "m"]))
    J <- rbind(J, data.frame(marker = v, outcome = o, n = a[["n"]], events = a[["events"]], hr_alone = a[["hr"]], lo_alone = a[["lo"]], hi_alone = a[["hi"]],
      hr_ksaa = kz[["hr"]], lo_ksaa = kz[["lo"]], hi_ksaa = kz[["hi"]], hr_marker = km[["hr"]], lo_marker = km[["lo"]], hi_marker = km[["hi"]])) } }
write.csv(J, file.path(W, "tableS7.csv"), row.names = FALSE)
for (i in seq_len(nrow(J))) with(J[i, ], cat(sprintf("  %-8s %-8s %5d/%5d  KSAA alone %s | joint KSAA %s, marker %s\n", marker, outcome, events, n, fci(hr_alone, lo_alone, hi_alone), fci(hr_ksaa, lo_ksaa, hi_ksaa), fci(hr_marker, lo_marker, hi_marker))))

# ---------------- Figure 3A: KSAA and PhenoAge acceleration in one model (age + sex), Wald test of the difference between the two log HRs ----------------
cat("\n=========== Figure 3A (weighted) ===========\n")
sp <- F$inc & !is.na(F$pa_accel)
F$zg <- (F$gap - mean(F$gap[sp])) / sd(F$gap[sp]); F$zm <- (F$pa_accel - mean(F$pa_accel[sp])) / sd(F$pa_accel[sp])
D6 <- svydesign(id = ~P6, strata = ~S6, weights = ~W6, nest = TRUE, data = F)
WD <- NULL
for (o in names(EV)) { m <- svycoxph(as.formula(sprintf("Surv(t, %s) ~ zg + zm + age + factor(sex)", EV[[o]])), design = dom6(sp)); b <- coef(m); V <- vcov(m)
  dif <- b[["zg"]] - b[["zm"]]; sed <- sqrt(V["zg", "zg"] + V["zm", "zm"] - 2 * V["zg", "zm"])
  WD <- rbind(WD, data.frame(outcome = o, n = m$n, events = m$nevent, t(setNames(ci3(b[["zg"]], sqrt(V["zg", "zg"])), c("hr_ksaa", "lo_ksaa", "hi_ksaa"))),
    t(setNames(ci3(b[["zm"]], sqrt(V["zm", "zm"])), c("hr_pa", "lo_pa", "hi_pa"))), ratio = exp(dif), ratio_lo = exp(dif - 1.96 * sed), ratio_hi = exp(dif + 1.96 * sed), p_diff = 2 * pnorm(-abs(dif / sed)))) }
write.csv(WD, file.path(W, "figure3a.csv"), row.names = FALSE)
for (i in seq_len(nrow(WD))) with(WD[i, ], cat(sprintf("  %-8s KSAA %s | PhenoAge %s | ratio %s  P %.2g\n", outcome, fci(hr_ksaa, lo_ksaa, hi_ksaa), fci(hr_pa, lo_pa, hi_pa), fci(ratio, ratio_lo, ratio_hi), p_diff)))

# ---------------- Figure 3C–E: cross-classified quartiles; KSAA per SD within each PhenoAge quartile (unweighted quartile cut points) ----------------
cat("\n=========== Figure 3C–E (weighted) ===========\n")
qg <- quantile(F$gap[sp], 0:4 / 4); qm <- quantile(F$pa_accel[sp], 0:4 / 4)
F$qk <- cut(F$gap, qg, include.lowest = TRUE, labels = QL); F$qp <- cut(F$pa_accel, qm, include.lowest = TRUE, labels = QL)
F$cell <- relevel(factor(ifelse(sp, paste0(F$qk, "_", F$qp), NA)), ref = "Q1_Q1")
D6 <- svydesign(id = ~P6, strata = ~S6, weights = ~W6, nest = TRUE, data = F)
JC <- NULL; JW <- NULL
for (o in names(EV)) { m <- svycoxph(as.formula(sprintf("Surv(t, %s) ~ cell + age + factor(sex)", EV[[o]])), design = dom6(sp))
  for (lv in levels(F$cell)) { k <- paste0("cell", lv); r <- if (lv == "Q1_Q1") c(hr = 1, lo = NA, hi = NA) else ci3(coef(m)[[k]], sqrt(vcov(m)[k, k]))
    JC <- rbind(JC, data.frame(outcome = o, y = sub("_.*", "", lv), x = sub(".*_", "", lv), n = sum(F$cell == lv, na.rm = TRUE), deaths = sum(F[[EV[[o]]]][which(F$cell == lv)]), hr = r[["hr"]], hr_lo = r[["lo"]], hr_hi = r[["hi"]])) }
  for (pq in QL) { a <- wcox(dom6(sp & F$qp %in% pq), "zg + age + factor(sex)", term = "zg", ev = EV[[o]], k = sd(F$gap[sp]) / SDG * 0 + 1)
    JW <- rbind(JW, data.frame(outcome = o, stratum = pq, hr = a[["hr"]], lo = a[["lo"]], hi = a[["hi"]])) } }
qhet <- function(w) { b <- log(w$hr); se <- (log(w$hi) - log(w$lo)) / (2 * qnorm(.975)); wt <- 1 / se^2; pchisq(sum(wt * (b - sum(wt * b) / sum(wt))^2), length(b) - 1, lower.tail = FALSE) }
write.csv(JC, file.path(W, "figure3_cells.csv"), row.names = FALSE); write.csv(JW, file.path(W, "figure3_within.csv"), row.names = FALSE)
for (o in c("allcause", "cvd")) { w <- JW[JW$outcome == o, ]; c44 <- JC[JC$outcome == o & JC$y == "Q4" & JC$x == "Q4", ]
  cat(sprintf("  %-8s KSAA per SD within strata %s–%s, P heterogeneity %.2f; Q4×Q4 %s\n", o, f2(min(w$hr)), f2(max(w$hr)), qhet(w), fci(c44$hr, c44$hr_lo, c44$hr_hi))) }

# ---------------- Figure S3: subgroups (weighted Model 3; radiographic OA uses the phase 2 design) ----------------
cat("\n=========== Figure S3 (weighted) ===========\n")
DX0 <- read.csv(file.path(SPW, "dxa_all.csv"))[, c("SEQN", "BDPFNBMD")]
q8b <- function(x) { x <- suppressWarnings(as.numeric(x)); for (k in 1:6) { s8 <- as.numeric(strrep("8", k)); lo <- x[!is.na(x) & x < s8]; if (any(x == s8, na.rm = TRUE) && length(lo) && s8 >= max(lo) * 1.5) x[x == s8] <- NA }; x }
DX0$BDPFNBMD <- q8b(DX0$BDPFNBMD)
REFD <- subset(merge(read.csv(file.path(SPW, "demo_all.csv")), DX0, by = "SEQN"), age >= 20 & age <= 29 & sex == 2 & race == 1 & !is.na(BDPFNBMD))
F$who <- cut((F$fn_bmd - mean(REFD$BDPFNBMD)) / sd(REFD$BDPFNBMD), c(-Inf, -2.5, -1, Inf), labels = c("osteo", "low", "normal"))
D6 <- svydesign(id = ~P6, strata = ~S6, weights = ~W6, nest = TRUE, data = F)
F2 <- F[F$phase_all == 2 & !is.na(F$W2) & F$W2 > 0, ]; D2 <- svydesign(id = ~P2, strata = ~S2, weights = ~W2, nest = TRUE, data = F2)
GR <- list(list("Sex", "factor(sex)", list(Men = cc & F$sex == 1, Women = cc & F$sex == 2), dom6),
  list("Age", NULL, list("60–74 years" = cc & F$age < 75, "≥75 years" = cc & F$age >= 75), dom6),
  list("Race and ethnicity", "factor(race)", list("Non-Hispanic White" = cc & F$race == 1, "Non-Hispanic Black" = cc & F$race == 2, "Mexican American" = cc & F$race == 3), dom6),
  list("BMI", NULL, list("<25" = cc & F$bmi < 25, "25–<30" = cc & F$bmi >= 25 & F$bmi < 30, "≥30" = cc & F$bmi >= 30), dom6),
  list("Radiographic OA (phase 2)", NULL, list("K–L <2" = p2 & F$klmax <= 1, "K–L ≥2" = p2 & F$klmax >= 2), dom2),
  list("Femoral neck BMD (WHO)", NULL, list(Normal = cc & F$who %in% "normal", "Low bone mass" = cc & F$who %in% "low", Osteoporosis = cc & F$who %in% "osteo"), dom6))
SG <- NULL
for (o in c("allcause", "cvd")) for (g in GR) { r <- do.call(rbind, lapply(names(g[[3]]), function(lv) { a <- wcox(g[[4]](g[[3]][[lv]]), paste("z +", paste(setdiff(COV, g[[2]]), collapse = " + ")), ev = EV[[o]])
    data.frame(outcome = o, fam = g[[1]], lev = lv, n = a[["n"]], hr = a[["hr"]], lo = a[["lo"]], hi = a[["hi"]]) }))
  r$p_het <- qhet(r); SG <- rbind(SG, r) }
ref <- do.call(rbind, lapply(c("allcause", "cvd"), function(o) { a <- wcox(dom6(cc), paste("z +", MODS[["m3"]]), ev = EV[[o]]); data.frame(outcome = o, fam = "All participants", lev = "", n = a[["n"]], hr = a[["hr"]], lo = a[["lo"]], hi = a[["hi"]], p_het = NA) }))
SG <- rbind(ref, SG); write.csv(SG, file.path(W, "figureS3.csv"), row.names = FALSE)
for (i in seq_len(nrow(SG))) with(SG[i, ], cat(sprintf("  %-8s %-26s %-20s %5d  %s  %s\n", outcome, fam, lev, n, fci(hr, lo, hi), if (is.na(p_het)) "" else sprintf("P %.3f", p_het))))

# ---------------- Table S10: fracture, falls and physical function (weighted logistic / linear) ----------------
cat("\n=========== Table S10 (weighted) ===========\n")
OUT <- c(fx_hip = "logistic", fx_frag = "logistic", fx_spi = "logistic", fx_wri = "logistic", fx_any = "logistic", fall = "logistic", gait = "linear",
         tand_fail = "logistic", dif_HAH1 = "logistic", dif_HAH2 = "logistic", dif_HAH8 = "logistic", dif_HAH9 = "logistic")
wglm <- function(des, y, rhs, type) { m <- svyglm(as.formula(paste(y, "~ z +", rhs)), design = des, family = if (type == "logistic") quasibinomial() else gaussian())
  b <- coef(m)[["z"]]; se <- sqrt(vcov(m)["z", "z"]); q <- qt(0.975, degf(des))
  if (type == "logistic") exp(c(b, b - q * se, b + q * se)) else 100 * c(b, b - q * se, b + q * se) }
S11 <- NULL
for (y in names(OUT)) { s <- cc & !is.na(F[[y]]); sb <- s & !is.na(F$tot_bmd)
  r <- sapply(names(MODS), function(mn) wglm(dom6(s), y, MODS[[mn]], OUT[[y]]))
  a <- wglm(dom6(sb), y, MODS[["m3"]], OUT[[y]]); b <- wglm(dom6(sb), y, paste(MODS[["m3"]], "+ tot_bmd"), OUT[[y]])
  ret <- if (OUT[[y]] == "linear") b[1] / a[1] else if (a[2] <= 1 && a[3] >= 1) NA else log(b[1]) / log(a[1])
  S11 <- rbind(S11, data.frame(y = y, type = OUT[[y]], n = sum(s), events = if (OUT[[y]] == "logistic") sum(F[[y]][s]) else NA,
    m1 = r[1, "m1"], m1_lo = r[2, "m1"], m1_hi = r[3, "m1"], m2 = r[1, "m2"], m2_lo = r[2, "m2"], m2_hi = r[3, "m2"], m3 = r[1, "m3"], m3_lo = r[2, "m3"], m3_hi = r[3, "m3"],
    n_b = sum(sb), events_b = if (OUT[[y]] == "logistic") sum(F[[y]][sb]) else NA, m3b0 = a[1], m3b0_lo = a[2], m3b0_hi = a[3], m3b = b[1], m3b_lo = b[2], m3b_hi = b[3], retained = ret)) }
write.csv(S11, file.path(W, "tableS10.csv"), row.names = FALSE)
for (i in seq_len(nrow(S11))) with(S11[i, ], cat(sprintf("  %-9s %4s/%4d  M1 %s  M3 %s | DXA subset %s -> +BMD %s  retained %s\n", y, ifelse(is.na(events), "—", events), n,
  fci(m1, m1_lo, m1_hi), fci(m3, m3_lo, m3_hi), fci(m3b0, m3b0_lo, m3b0_hi), fci(m3b, m3b_lo, m3b_hi), ifelse(is.na(retained), "—", sprintf("%.0f%%", 100 * retained)))))

cat("\n-> results/analysis/weighted/\n")

# ---------------- Figure 1D: weighted HR in the full cohort, age + sex (as for the Figure 2 curves) ----------------
M1F <- do.call(rbind, lapply(names(EV), function(o) { a <- wcox(dom6(F$inc), "z + age + factor(sex)", ev = EV[[o]]); data.frame(outcome = o, n = a[["n"]], events = a[["events"]], hr_sd = a[["hr"]], lo_sd = a[["lo"]], hi_sd = a[["hi"]]) }))
write.csv(M1F, file.path(W, "model1_full_cohort.csv"), row.names = FALSE)
for (i in seq_len(nrow(M1F))) with(M1F[i, ], cat(sprintf("  Figure 1D full cohort, age + sex %-8s %s\n", outcome, fci(hr_sd, lo_sd, hi_sd))))

# ---------------- Table S10 note: additional adjustment of gait speed (weighted; knee pain, or phase 2 K–L grade) ----------------
wgait <- function(des, rhs) { m <- svyglm(as.formula(paste("gait ~ z +", rhs)), design = des); b <- coef(m)[["z"]]; se <- sqrt(vcov(m)["z", "z"]); q <- qt(0.975, degf(des)); 100 * c(b, b - q * se, b + q * se) }
sp <- F$cc & !is.na(F$gait) & !is.na(F$pain); sk <- F$cc & !is.na(F$gait) & !is.na(F$klmax)
gp <- wgait(dom6(sp), paste(MODS[["m3"]], "+ pain")); gk <- wgait(dom2(sk), paste(MODS[["m3"]], "+ klmax"))
GX <- data.frame(item = c("gait_m3_pain", "gait_m3_kl"), n = c(sum(sp), sum(sk)), est = c(gp[1], gk[1]), lo = c(gp[2], gk[2]), hi = c(gp[3], gk[3]))
write.csv(GX, file.path(W, "tableS10_gait_extra.csv"), row.names = FALSE)
for (i in 1:2) with(GX[i, ], cat(sprintf("  gait speed %-13s n %d  %.1f (%.1f to %.1f) cm/s\n", item, n, est, lo, hi)))
