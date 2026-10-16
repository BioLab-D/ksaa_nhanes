# 稿件需要、而下游脚本没有直接打印的数字(ksaa_oof.csv,4,875 人):
#  1 随访汇总;2 KSAA 连续与中位数分组 × 模型 1–3(人数/死亡/ΔC,四分位与 log-rank);3 PhenoAge 模型上的 C 增量;
#  4 种族与贫困的效应修饰(分层拟合后比较 log HR,不用乘积项);5 平凡基线 MAE 与原始集成的相对降幅
suppressMessages(library(survival))
res <- here::here("results"); sp <- file.path(res,"work")
d <- read.csv(file.path(res,"ksaa_oof.csv")); sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN)
for(v in c("sex","race","bmi","fnbmd","phase","pir","WTPFEX6")) d[[v]] <- sl[[v]][i]
d$t <- sl$PERMTH_EXM[i]; d$e <- sl$MORTSTAT[i]; d$gap <- d$gap_corr; d$z <- d$gap/sd(d$gap)
V <- read.csv(file.path(sp,"validate_vars.csv")); j <- match(d$SEQN, V$SEQN); q8 <- function(x) ifelse(x %in% c(8,9), NA, x)
d$dm <- as.integer(q8(V$HAD1[j])==1); d$arth <- as.integer(q8(V$HAC1A[j])==1); d$stroke <- as.integer(q8(V$HAC1D[j])==1)
d$mi <- as.integer(q8(V$HAF10[j])==1); d$smk <- as.integer(q8(V$HAR1[j])==1)
f3 <- function(s) sprintf("%.3f (%.3f–%.3f)", s[1], s[3], s[4])

cat("=========== 1. 随访 ===========\n")
rk <- survfit(Surv(t, 1-e) ~ 1, d)
cat(sprintf("n %d;死亡 %d (%.1f%%);人年 %.0f;观察随访中位 %.1f 年(最大 %.1f);反向 KM 潜在随访中位 %.1f 年\n",
    nrow(d), sum(d$e), 100*mean(d$e), sum(d$t)/12, median(d$t)/12, max(d$t)/12, summary(rk)$table["median"]/12))

cat("\n=========== 2. 表 3:膝龄差与全因死亡 ===========\n")
M <- list(m1="age+factor(sex)+factor(race)", m2="age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi", m3="age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi+fnbmd")
med <- median(d$gap); d$hi <- as.integer(d$gap > med)
cat(sprintf("中位数切点 %.2f 岁\n", med))
for(g in 0:1){ x <- d[d$hi==g,]; cat(sprintf("  %s  n %d  死亡 %d (%.1f%%)  gap %.2f ± %.2f\n", c("低风险(≤中位数)","高风险(>中位数)")[g+1], nrow(x), sum(x$e), 100*mean(x$e), mean(x$gap), sd(x$gap))) }
for(nm in names(M)){ vars <- all.vars(as.formula(paste("~",M[[nm]]))); s <- d[complete.cases(d[,c("t","e","gap",vars)]),]
  m0 <- coxph(as.formula(paste("Surv(t,e)~",M[[nm]])), s); my <- coxph(as.formula(paste("Surv(t,e)~gap+",M[[nm]])), s)
  ms <- coxph(as.formula(paste("Surv(t,e)~z+",M[[nm]])), s); mh <- coxph(as.formula(paste("Surv(t,e)~hi+",M[[nm]])), s)
  cat(sprintf("  %s n %d 死亡 %d | 每年 %s | 每 SD %s P %.1e | 高 vs 低 %s P %.1e | ΔC(连续)%+.4f\n", nm, nrow(s), sum(s$e),
      f3(summary(my)$conf.int["gap",]), f3(summary(ms)$conf.int["z",]), summary(ms)$coef["z",5], f3(summary(mh)$conf.int["hi",]), summary(mh)$coef["hi",5],
      concordance(ms)$concordance - concordance(m0)$concordance)) }
d$q <- cut(d$gap, quantile(d$gap, 0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))
sq <- summary(coxph(Surv(t,e)~q+age+factor(sex)+factor(race), d))$conf.int
cat(sprintf("  四分位(模型 1):Q2 %s  Q3 %s  Q4 %s  log-rank P %.1e\n", f3(sq[1,]), f3(sq[2,]), f3(sq[3,]), pchisq(survdiff(Surv(t,e)~q, d)$chisq, 3, lower.tail=FALSE)))
cat(sprintf("  各四分位 gap 范围:%s\n", paste(levels(d$q), collapse="  ")))

cat("\n=========== 3. PhenoAge 模型上的 C 增量(同一 n)===========\n")
ph <- readRDS(file.path(sp,"oof_pheno.rds")); d$pz <- as.numeric(scale(ph$pa_accel[match(d$SEQN, ph$SEQN)]))
s <- d[complete.cases(d[,c("t","e","age","sex","race","pz","z")]),]
c0 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race), s))$concordance
c1 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race)+pz, s))$concordance
c2 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race)+pz+z, s))$concordance
c3 <- concordance(coxph(Surv(t,e)~age+factor(sex)+factor(race)+z, s))$concordance
cat(sprintf("  n %d:年龄+性别+种族 C %.4f;+ 膝龄差 %.4f;+ PhenoAge %.4f;+ 两者 %.4f(在 PhenoAge 之上 %+.4f)\n", nrow(s), c0, c3, c1, c2, c2-c1))

cat("\n=========== 4. 效应修饰:分层拟合比 log HR ===========\n")
strat <- function(sub, lab){ x <- d[sub,]; m <- coxph(Surv(t,e)~z+age+factor(sex), x); b <- coef(m)["z"]; se <- sqrt(vcov(m)["z","z"])
  cat(sprintf("  %-22s n %4d  HR %s\n", lab, nrow(x), f3(summary(m)$conf.int["z",]))); c(b=b, se=se) }
R3 <- rbind(W=strat(d$race==1,"非西班牙裔白人"), B=strat(d$race==2,"非西班牙裔黑人"), M=strat(d$race==3,"墨西哥裔美国人"))
w <- 1/R3[,2]^2; bb <- sum(w*R3[,1])/sum(w); Q <- sum(w*(R3[,1]-bb)^2)
zt <- function(a,b) { z <- (R3[a,1]-R3[b,1])/sqrt(R3[a,2]^2+R3[b,2]^2); sprintf("z %.2f P %.3f", z, 2*pnorm(-abs(z))) }
cat(sprintf("  三组异质性 Cochran Q %.2f,df 2,P %.4f;白人 vs 黑人 %s;白人 vs 墨西哥裔 %s\n", Q, pchisq(Q,2,lower.tail=FALSE), zt("W","B"), zt("W","M")))
P2 <- rbind(lo=strat(!is.na(d$pir) & d$pir<1,"PIR < 1"), hi=strat(!is.na(d$pir) & d$pir>4,"PIR > 4"))
z <- (P2[1,1]-P2[2,1])/sqrt(sum(P2[,2]^2)); cat(sprintf("  PIR<1 vs >4:z %.2f P %.3f\n", z, 2*pnorm(-abs(z))))

cat("\n=========== 5. 平凡基线 ===========\n")
bm <- sapply(0:4, function(k){ tr <- d$fold!=k; te <- d$fold==k; mean(abs(median(d$age[tr]) - d$age[te])) })
rm <- sapply(0:4, function(k){ te <- d$fold==k; mean(abs(d$pred_raw[te]-d$age[te])) })
cat(sprintf("  训练折中位数预测 MAE %.3f ± %.3f;原始集成 %.3f;相对降幅 %.1f%%\n", mean(bm), sd(bm), mean(rm), 100*(1-mean(rm)/mean(bm))))
cat(sprintf("  合并 OOF:原始 r %.3f;校正后 r %.3f\n", cor(d$pred_raw, d$age), cor(d$pred_corr, d$age)))
