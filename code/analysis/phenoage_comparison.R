# KSAA and PhenoAge acceleration: head-to-head comparison (part 1) and combined PhenoAge + KSAA score (part 2).
.fresh_rng <- function() { RNGkind("default", "default", "default")   # reset the RNG to the state of a fresh R session at the start of each part
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }

# ======== R6: head-to-head comparison of KSAA and PhenoAge acceleration (joint model + Wald test of the coefficient difference) ========
.fresh_rng(); local({
# R6: head-to-head comparison of KSAA and PhenoAge acceleration (joint model + Wald test of the coefficient difference)
# Design follows Chen et al. Cell Rep Methods 2025;5:101115, Fig. 5D.
# Sample = participants with pa_accel (n = 4,522); each measure standardized within this sample; adjusted for age + sex.
#   Panel D: both in one model, Wald test of the difference between the two coefficients (joint covariance, accounts for their correlation)
# Input results/analysis/analysis_data.rds; output results/analysis/phenoage_wald.csv
suppressMessages(library(survival))
res <- here::here("results"); out <- file.path(res, "analysis")
d <- readRDS(file.path(out, "analysis_data.rds"))
fp <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))

s <- d[!is.na(d$pa_accel), ]
sdz <- function(x) (x - mean(x)) / sd(x)
s$zg <- sdz(s$gap); s$zm <- sdz(s$pa_accel)
EV <- c("All-cause mortality" = "e", "CVD mortality" = "e_cvd", "Cancer mortality" = "e_cancer")
cat(sprintf("Sample n = %d; KSAA SD = %.2f y, PhenoAge acceleration SD = %.2f y; r = %.3f\n\n",
            nrow(s), sd(s$gap), sd(s$pa_accel), cor(s$gap, s$pa_accel)))

# ---------- Panel D: joint model + Wald test of the difference between the two coefficients ----------
cat("\n=========== D. Both in one model (per SD), testing the difference between the two log HRs ===========\n")
WD <- NULL
for (o in names(EV)) {
  j <- coxph(as.formula(sprintf("Surv(t, %s) ~ zg + zm + age + factor(sex)", EV[o])), s)
  b <- coef(j); V <- vcov(j)
  dif <- b[["zg"]] - b[["zm"]]; sed <- sqrt(V["zg", "zg"] + V["zm", "zm"] - 2 * V["zg", "zm"])
  z <- dif / sed; p <- 2 * pnorm(-abs(z))
  WD <- rbind(WD, data.frame(outcome = o, n = nrow(s), events = sum(s[[EV[o]]]),
    hr_ksaa = exp(b[["zg"]]), lo_ksaa = exp(b[["zg"]] - 1.96 * sqrt(V["zg", "zg"])), hi_ksaa = exp(b[["zg"]] + 1.96 * sqrt(V["zg", "zg"])),
    hr_pa = exp(b[["zm"]]), lo_pa = exp(b[["zm"]] - 1.96 * sqrt(V["zm", "zm"])), hi_pa = exp(b[["zm"]] + 1.96 * sqrt(V["zm", "zm"])),
    ratio = exp(dif), ratio_lo = exp(dif - 1.96 * sed), ratio_hi = exp(dif + 1.96 * sed), z = z, p_diff = p))
}
write.csv(WD, file.path(out, "phenoage_wald.csv"), row.names = FALSE)
for (i in seq_len(nrow(WD))) with(WD[i, ], cat(sprintf("  %-20s n %4d events %4d | KSAA %.2f (%.2f–%.2f) | PhenoAge %.2f (%.2f–%.2f) | ratio %.2f (%.2f–%.2f) Wald P %s\n",
  outcome, n, events, hr_ksaa, lo_ksaa, hi_ksaa, hr_pa, lo_pa, hi_pa, ratio, ratio_lo, ratio_hi, fp(p_diff))))
cat(sprintf("\n  Note: ratio = KSAA HR / PhenoAge HR (same model); the Wald test uses the joint covariance and accounts for their correlation (r = %.2f)\n", cor(s$gap, s$pa_accel)))
cat("\n-> phenoage_wald.csv\n")
})

# ======== R7: PhenoAge alone vs combined PhenoAge + KSAA score (HR per SD / by category, C and R²_D, bootstrap CI of the improvement) ========
.fresh_rng(); local({
# R7: PhenoAge alone vs combined PhenoAge + KSAA score (HR per SD / by category, C and R²_D, bootstrap CI of the improvement)
# Structure follows Chen et al. Cell Rep Methods 2025;5:101115, Fig. 5A–C. A KSAA-only score is also computed, for reference only.
# Score = linear combination of the marker terms in a Surv ~ markers + age + sex model, standardized to SD = 1 and entered as a single variable.
# Combination weights are estimated in the same participants (apparent values, no optimism correction), so improvements need CIs: each bootstrap resample re-estimates the weights and the quantile cut points.
# Sample = participants with pa_accel, n = 4,522. Input results/analysis/analysis_data.rds; outputs phenoage_combined_scores.csv, phenoage_combined_cindex.csv, phenoage_combined_improvement.csv
suppressMessages({library(survival); library(parallel)})
RNGkind("L'Ecuyer-CMRG"); set.seed(20260915); NC <- max(1, detectCores() - 1); B <- 200
res <- here::here("results"); out <- file.path(res, "analysis")
d <- readRDS(file.path(out, "analysis_data.rds"))

s <- d[!is.na(d$pa_accel), ]
EV <- c("All-cause mortality" = "e", "CVD mortality" = "e_cvd", "Cancer mortality" = "e_cancer")
SCORE <- list("Phenotypic Age" = "pa_accel", "KSAA" = "gap", "Phenotypic Age + KSAA" = c("pa_accel", "gap"))
ADJ <- "age + factor(sex)"
cat(sprintf("Sample n = %d; KSAA SD %.2f y, PhenoAge acceleration SD %.2f y, r = %.3f\n", nrow(s), sd(s$gap), sd(s$pa_accel), cor(s$gap, s$pa_accel)))

mkscore <- function(dat, ev, vars) {   # coefficients estimated in dat; only the marker part of the linear predictor
  m <- coxph(as.formula(sprintf("Surv(t, %s) ~ %s + %s", ev, paste(vars, collapse = " + "), ADJ)), dat)
  z <- as.numeric(as.matrix(dat[, vars, drop = FALSE]) %*% coef(m)[vars])
  (z - mean(z)) / sd(z)
}
hr_sd_of <- function(dat, ev) {        # dat must already have an sc column
  m <- summary(coxph(as.formula(sprintf("Surv(t, %s) ~ sc + %s", ev, ADJ)), dat))$coefficients["sc", ]
  b <- as.numeric(m["coef"]); se <- as.numeric(m["se(coef)"]); c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
}
catvar <- function(x, scheme) {
  if (scheme == "q5") { q <- cut(x, quantile(x, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
    factor(ifelse(q == 5, "top", ifelse(q %in% 2:4, "mid", NA)), c("mid", "top"))
  } else { p <- quantile(x, c(.25, .75, .95))
    factor(ifelse(x >= p[3], "top", ifelse(x >= p[1] & x <= p[2], "mid", NA)), c("mid", "top")) }
}
hr_cat_of <- function(dat, ev, scheme) {
  dat$cg <- catvar(dat$sc, scheme); ss <- dat[!is.na(dat$cg), ]
  m <- summary(coxph(as.formula(sprintf("Surv(t, %s) ~ cg + %s", ev, ADJ)), ss))$coefficients["cgtop", ]
  b <- as.numeric(m["coef"]); se <- as.numeric(m["se(coef)"])
  c(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), ntop = sum(ss$cg == "top"), nref = sum(ss$cg == "mid"))
}

# ---------- HR per SD and by category (all three scores) ----------
S <- NULL
for (o in names(EV)) { ev <- EV[[o]]
  for (sn in names(SCORE)) {
    s$sc <- mkscore(s, ev, SCORE[[sn]]); h <- hr_sd_of(s, ev)
    row <- data.frame(outcome = o, score = sn, n = nrow(s), events = sum(s[[ev]]),
                      hr_sd = h[["hr"]], lo_sd = h[["lo"]], hi_sd = h[["hi"]])
    for (scm in c("q5", "p95")) { k <- hr_cat_of(s, ev, scm)
      row[[paste0("hr_", scm)]] <- k[["hr"]]; row[[paste0("lo_", scm)]] <- k[["lo"]]; row[[paste0("hi_", scm)]] <- k[["hi"]]
      row[[paste0("ntop_", scm)]] <- k[["ntop"]]; row[[paste0("nref_", scm)]] <- k[["nref"]] }
    S <- rbind(S, row)
  }
}
write.csv(S, file.path(out, "phenoage_combined_scores.csv"), row.names = FALSE)

# ---------- Discrimination + explained variation ----------
# Royston–Sauerbrei R²_D (R.D of survival::royston), like C a whole-model quantity, computed in the same resample; relative gain = combined / PhenoAge − 1 (%)
cfit <- function(dat, ev) {
  f <- lapply(SCORE, function(v) coxph(as.formula(sprintf("Surv(t, %s) ~ %s + %s", ev, paste(v, collapse = " + "), ADJ)), dat))
  c(setNames(sapply(f, function(m) concordance(m)$concordance), c("C_pa", "C_ks", "C_cb")),
    setNames(sapply(f, function(m) royston(m)[["R.D"]]), c("R2_pa", "R2_ks", "R2_cb")))
}
C <- NULL
for (o in names(EV)) { ev <- EV[[o]]
  c0 <- cfit(s, ev)
  bt <- simplify2array(mclapply(1:B, function(i) cfit(s[sample.int(nrow(s), replace = TRUE), ], ev), mc.cores = NC))
  q <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
  rel <- function(a, b) 100 * (a / b - 1)
  C <- rbind(C, data.frame(outcome = o, n = nrow(s), events = sum(s[[ev]]),
    C_pa = c0[["C_pa"]], C_pa_lo = q(bt["C_pa", ])[1], C_pa_hi = q(bt["C_pa", ])[2],
    C_ks = c0[["C_ks"]], C_ks_lo = q(bt["C_ks", ])[1], C_ks_hi = q(bt["C_ks", ])[2],
    C_cb = c0[["C_cb"]], C_cb_lo = q(bt["C_cb", ])[1], C_cb_hi = q(bt["C_cb", ])[2],
    dC_vs_pa = c0[["C_cb"]] - c0[["C_pa"]], dC_vs_pa_lo = q(bt["C_cb", ] - bt["C_pa", ])[1], dC_vs_pa_hi = q(bt["C_cb", ] - bt["C_pa", ])[2],
    relC = rel(c0[["C_cb"]], c0[["C_pa"]]), relC_lo = q(rel(bt["C_cb", ], bt["C_pa", ]))[1], relC_hi = q(rel(bt["C_cb", ], bt["C_pa", ]))[2],
    R2_pa = c0[["R2_pa"]], R2_pa_lo = q(bt["R2_pa", ])[1], R2_pa_hi = q(bt["R2_pa", ])[2],
    R2_ks = c0[["R2_ks"]], R2_ks_lo = q(bt["R2_ks", ])[1], R2_ks_hi = q(bt["R2_ks", ])[2],
    R2_cb = c0[["R2_cb"]], R2_cb_lo = q(bt["R2_cb", ])[1], R2_cb_hi = q(bt["R2_cb", ])[2],
    dR2 = c0[["R2_cb"]] - c0[["R2_pa"]], dR2_lo = q(bt["R2_cb", ] - bt["R2_pa", ])[1], dR2_hi = q(bt["R2_cb", ] - bt["R2_pa", ])[2],
    relR2 = rel(c0[["R2_cb"]], c0[["R2_pa"]]), relR2_lo = q(rel(bt["R2_cb", ], bt["R2_pa", ]))[1], relR2_hi = q(rel(bt["R2_cb", ], bt["R2_pa", ]))[2]))
}
rownames(C) <- NULL
write.csv(C, file.path(out, "phenoage_combined_cindex.csv"), row.names = FALSE)

# ---------- Bootstrap CI of the improvement (combined / PhenoAge): each resample re-estimates the weights and cut points ----------
one_boot <- function(x, ev) {
  x$sc <- mkscore(x, ev, SCORE[[1]])
  a_sd <- hr_sd_of(x, ev)[["hr"]]; a_95 <- hr_cat_of(x, ev, "p95")[["hr"]]; a_q5 <- hr_cat_of(x, ev, "q5")[["hr"]]
  x$sc <- mkscore(x, ev, SCORE[[3]])
  b_sd <- hr_sd_of(x, ev)[["hr"]]; b_95 <- hr_cat_of(x, ev, "p95")[["hr"]]; b_q5 <- hr_cat_of(x, ev, "q5")[["hr"]]
  c(per_sd = b_sd / a_sd, p95 = b_95 / a_95, q5 = b_q5 / a_q5)
}
IMP <- NULL
for (o in names(EV)) { ev <- EV[[o]]
  obs <- one_boot(s, ev)
  bt <- simplify2array(mclapply(1:B, function(i) tryCatch(one_boot(s[sample.int(nrow(s), replace = TRUE), ], ev),
                                                          error = function(e) c(per_sd = NA, p95 = NA, q5 = NA)), mc.cores = NC))
  for (k in c("per_sd", "p95", "q5")) { v <- bt[k, ]; qq <- quantile(v, c(.025, .975), na.rm = TRUE)
    IMP <- rbind(IMP, data.frame(outcome = o, metric = k, ratio = obs[[k]], ratio_lo = qq[1], ratio_hi = qq[2],
                                 pct = 100 * (obs[[k]] - 1), pct_lo = 100 * (qq[1] - 1), pct_hi = 100 * (qq[2] - 1),
                                 n_ok = sum(!is.na(v)))) }
}
write.csv(IMP, file.path(out, "phenoage_combined_improvement.csv"), row.names = FALSE)

cat("\n=========== A. Discrimination (Harrell C, age + sex + score) ===========\n")
for (i in seq_len(nrow(C))) with(C[i, ], cat(sprintf("  %-20s PhenoAge %.4f | combined %.4f | ΔC %+.4f (%+.4f, %+.4f)   [KSAA alone %.4f, for reference]\n",
  outcome, C_pa, C_cb, dC_vs_pa, dC_vs_pa_lo, dC_vs_pa_hi, C_ks)))
for (i in seq_len(nrow(C))) with(C[i, ], cat(sprintf("  %-20s relative C %+.2f%% (%+.2f%%, %+.2f%%) | R²_D PhenoAge %.4f -> combined %.4f  relative %+.1f%% (%+.1f%%, %+.1f%%)   [KSAA alone %.4f]\n",
  outcome, relC, relC_lo, relC_hi, R2_pa, R2_cb, relR2, relR2_lo, relR2_hi, R2_ks)))
cat("\n=========== B/C. HR per SD and HR for the top 5% ===========\n")
for (i in seq_len(nrow(S))) with(S[i, ], cat(sprintf("  %-20s %-22s per SD %.2f (%.2f–%.2f) | >=95th percentile %.2f (%.2f–%.2f)\n",
  outcome, score, hr_sd, lo_sd, hi_sd, hr_p95, lo_p95, hi_p95)))
cat("\n=========== Improvement (combined vs PhenoAge), bootstrap CI ===========\n")
for (i in seq_len(nrow(IMP))) with(IMP[i, ], cat(sprintf("  %-20s %-7s %+.1f%% (%+.1f%%, %+.1f%%)%s\n",
  outcome, metric, pct, pct_lo, pct_hi, if (pct_lo <= 0 && pct_hi >= 0) "   <- CI includes 0, improvement not significant" else "")))
cat(sprintf("\n  Combination weights and quantile cut points re-estimated within each resample (B = %d); the combined score is apparent, without optimism correction\n", B))
cat("\n-> phenoage_combined_scores.csv, phenoage_combined_cindex.csv, phenoage_combined_improvement.csv\n")
})
