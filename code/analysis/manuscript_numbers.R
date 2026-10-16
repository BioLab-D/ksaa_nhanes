# Numbers needed by the manuscript that downstream scripts do not print directly (ksaa_oof.csv, 4,875 participants):
#  1 follow-up summary; 2 continuous and median-split KSAA x Models 1–3 (n/deaths/ΔC, quartiles and log-rank); 3 C increment over the PhenoAge model;
#  4 effect modification by race and poverty (log HRs compared across stratified fits, no product terms); 5 trivial-baseline MAE and relative reduction by the raw ensemble
suppressMessages(library(survival))
res <- here::here("results"); sp <- file.path(res,"work")
d <- read.csv(file.path(res,"ksaa_oof.csv")); sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN)
for(v in c("sex","race","bmi","fnbmd","phase","pir","WTPFEX6")) d[[v]] <- sl[[v]][i]
d$t <- sl$PERMTH_EXM[i]; d$e <- sl$MORTSTAT[i]; d$gap <- d$gap_corr; d$z <- d$gap/sd(d$gap)
V <- read.csv(file.path(sp,"validate_vars.csv")); j <- match(d$SEQN, V$SEQN); q8 <- function(x) ifelse(x %in% c(8,9), NA, x)
d$dm <- as.integer(q8(V$HAD1[j])==1); d$arth <- as.integer(q8(V$HAC1A[j])==1); d$stroke <- as.integer(q8(V$HAC1D[j])==1)
d$mi <- as.integer(q8(V$HAF10[j])==1); d$smk <- as.integer(q8(V$HAR1[j])==1)
f3 <- function(s) sprintf("%.3f (%.3f–%.3f)", s[1], s[3], s[4])

cat("=========== 1. Follow-up ===========\n")
rk <- survfit(Surv(t, 1-e) ~ 1, d)
cat(sprintf("n %d; deaths %d (%.1f%%); person-years %.0f; median observed follow-up %.1f years (max %.1f); reverse-KM median potential follow-up %.1f years\n",
    nrow(d), sum(d$e), 100*mean(d$e), sum(d$t)/12, median(d$t)/12, max(d$t)/12, summary(rk)$table["median"]/12))

cat("\n=========== 2. Table 3: knee age gap and all-cause mortality ===========\n")
M <- list(m1="age+factor(sex)+factor(race)", m2="age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi", m3="age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi+fnbmd")
med <- median(d$gap); d$hi <- as.integer(d$gap > med)
cat(sprintf("Median cut point %.2f years\n", med))
for(g in 0:1){ x <- d[d$hi==g,]; cat(sprintf("  %s  n %d  deaths %d (%.1f%%)  gap %.2f ± %.2f\n", c("Low (<= median)","High (> median)")[g+1], nrow(x), sum(x$e), 100*mean(x$e), mean(x$gap), sd(x$gap))) }
for(nm in names(M)){ vars <- all.vars(as.formula(paste("~",M[[nm]]))); s <- d[complete.cases(d[,c("t","e","gap",vars)]),]
  m0 <- coxph(as.formula(paste("Surv(t,e)~",M[[nm]])), s); my <- coxph(as.formula(paste("Surv(t,e)~gap+",M[[nm]])), s)
  ms <- coxph(as.formula(paste("Surv(t,e)~z+",M[[nm]])), s); mh <- coxph(as.formula(paste("Surv(t,e)~hi+",M[[nm]])), s)
  cat(sprintf("  %s n %d deaths %d | per year %s | per SD %s P %.1e | high vs low %s P %.1e | ΔC (continuous) %+.4f\n", nm, nrow(s), sum(s$e),
      f3(summary(my)$conf.int["gap",]), f3(summary(ms)$conf.int["z",]), summary(ms)$coef["z",5], f3(summary(mh)$conf.int["hi",]), summary(mh)$coef["hi",5],
      concordance(ms)$concordance - concordance(m0)$concordance)) }
d$q <- cut(d$gap, quantile(d$gap, 0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))
sq <- summary(coxph(Surv(t,e)~q+age+factor(sex)+factor(race), d))$conf.int
cat(sprintf("  Quartiles (Model 1): Q2 %s  Q3 %s  Q4 %s  log-rank P %.1e\n", f3(sq[1,]), f3(sq[2,]), f3(sq[3,]), pchisq(survdiff(Surv(t,e)~q, d)$chisq, 3, lower.tail=FALSE)))
cat(sprintf("  Gap range by quartile: %s\n", paste(levels(d$q), collapse="  ")))

cat("\n=========== 3. C increment over the PhenoAge model (same n) ===========\n")
ph <- readRDS(file.path(sp,"oof_pheno.rds")); d$pz <- as.numeric(scale(ph$pa_accel[match(d$SEQN, ph$SEQN)]))
s <- d[complete.cases(d[,c("t","e","age","sex","race","pz","z")]),]
c0 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race), s))$concordance
c1 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race)+pz, s))$concordance
c2 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race)+pz+z, s))$concordance
c3 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race)+z, s))$concordance
cat(sprintf("  n %d: age+sex+race C %.4f; + knee age gap %.4f; + PhenoAge %.4f; + both %.4f (over PhenoAge %+.4f)\n", nrow(s), c0, c3, c1, c2, c2-c1))

cat("\n=========== 4. Effect modification: log HRs from stratified fits ===========\n")
strat <- function(sub, lab){ x <- d[sub,]; m <- coxph(Surv(t,e)~z+age+factor(sex), x); b <- coef(m)["z"]; se <- sqrt(vcov(m)["z","z"])
  cat(sprintf("  %-22s n %4d  HR %s\n", lab, nrow(x), f3(summary(m)$conf.int["z",]))); c(b=b, se=se) }
R3 <- rbind(W=strat(d$race==1,"Non-Hispanic White"), B=strat(d$race==2,"Non-Hispanic Black"), M=strat(d$race==3,"Mexican American"))
w <- 1/R3[,2]^2; bb <- sum(w*R3[,1])/sum(w); Q <- sum(w*(R3[,1]-bb)^2)
zt <- function(a,b) { z <- (R3[a,1]-R3[b,1])/sqrt(R3[a,2]^2+R3[b,2]^2); sprintf("z %.2f P %.3f", z, 2*pnorm(-abs(z))) }
cat(sprintf("  Heterogeneity across three groups Cochran Q %.2f, df 2, P %.4f; White vs Black %s; White vs Mexican American %s\n", Q, pchisq(Q,2,lower.tail=FALSE), zt("W","B"), zt("W","M")))
P2 <- rbind(lo=strat(!is.na(d$pir) & d$pir<1,"PIR < 1"), hi=strat(!is.na(d$pir) & d$pir>4,"PIR > 4"))
z <- (P2[1,1]-P2[2,1])/sqrt(sum(P2[,2]^2)); cat(sprintf("  PIR<1 vs >4: z %.2f P %.3f\n", z, 2*pnorm(-abs(z))))

cat("\n=========== 5. Trivial baseline ===========\n")
bm <- sapply(0:4, function(k){ tr <- d$fold!=k; te <- d$fold==k; mean(abs(median(d$age[tr]) - d$age[te])) })
rm <- sapply(0:4, function(k){ te <- d$fold==k; mean(abs(d$pred_raw[te]-d$age[te])) })
cat(sprintf("  Training-fold median prediction MAE %.3f ± %.3f; raw ensemble %.3f; relative reduction %.1f%%\n", mean(bm), sd(bm), mean(rm), 100*(1-mean(rm)/mean(bm))))
cat(sprintf("  Pooled OOF: raw r %.3f; corrected r %.3f\n", cor(d$pred_raw, d$age), cor(d$pred_corr, d$age)))
