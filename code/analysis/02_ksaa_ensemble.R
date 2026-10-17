suppressMessages(library(survival))
res <- here::here("results")
sp  <- here::here("results", "work")
o <- read.csv(file.path(res,"predictions/cv_resnet18_oof.csv"))
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(o$SEQN, sl$SEQN)
stopifnot(!anyNA(i), all(o$age==sl$age[i]))
d <- data.frame(SEQN=o$SEQN, age=o$age, t=sl$PERMTH_EXM[i], e=sl$MORTSTAT[i], sex=sl$sex[i], fold=o$fold)
d$race <- sl$race[i]; d$bmi <- sl$bmi[i]; d$fnbmd <- sl$fnbmd[i]
V <- read.csv(file.path(sp,"validate_vars.csv")); j <- match(d$SEQN, V$SEQN)
q8 <- function(x) ifelse(x %in% c(8,9), NA, x)
d$dm <- as.integer(q8(V$HAD1[j])==1); d$arth <- as.integer(q8(V$HAC1A[j])==1); d$stroke <- as.integer(q8(V$HAC1D[j])==1)
d$mi <- as.integer(q8(V$HAF10[j])==1); d$smk <- as.integer(q8(V$HAR1[j])==1)
al <- function(f){ x <- read.csv(f); x$pred_age[match(d$SEQN, x$SEQN)] }
# each pre-specified backbone was trained once (split seed 20260831); none was selected on out-of-fold performance
pf <- function(a, s="oof") file.path(res, "predictions", sprintf("cv_%s_%s.csv", a, s))
A5 <- c("ResNet-18"="resnet18", "ResNet-50"="resnet50", "EfficientNet-B0"="efficientnetb0", "ConvNeXt-Tiny"="convnexttiny",
        "DenseNet-121"="densenet121")
cfg <- lapply(A5, function(a) al(pf(a)))
FM <- sapply(A5, function(a){   # internal-validation MAE by fold
  f <- read.csv(pf(a, "folds")); stopifnot(identical(f$fold, 0:4)); f$inner_val_mae
})
M <- do.call(cbind, cfg); stopifnot(!anyNA(M)); n <- nrow(d); K <- 0:4; age <- d$age

cat("=========== 1. Within-fold weights (prop. to 1/internal-validation MAE; fold k's held-out participants use only fold k's weights) ===========\n")
# Fold k uses only W[, k]: the five-fold mean weight w would include information from the held-out participants and is only reported.
# Equal weights are a sensitivity analysis.
W <- apply(FM, 1, function(m) (1/m)/sum(1/m)); rownames(W) <- names(cfg)
colnames(W) <- paste0("fold",K); w <- rowMeans(W)
print(round(cbind(W, mean=w),4))
d$pred_raw <- NA_real_
for(k in K){ i <- d$fold==k; d$pred_raw[i] <- as.numeric(M[i,,drop=FALSE] %*% W[,k+1]) }
d$pred_eq <- rowMeans(M)

cat("\n=========== 2. Bias correction (reference set = the four folds other than the held-out fold) ===========\n")
d$pred_corr <- d$pred_corr_s <- NA
AB <- matrix(NA,5,2,dimnames=list(paste0("fold",K),c("b(intercept)","a(slope)"))); ABs <- AB
for(k in K){
  te <- d$fold==k; tr <- !te; r <- d[tr,]; r$gap <- r$pred_raw - r$age
  f <- lm(gap ~ age, r); AB[k+1,] <- coef(f)
  d$pred_corr[te] <- d$pred_raw[te] - (coef(f)[2]*age[te] + coef(f)[1])
  rs <- r[r$t >= 120,]; f2 <- lm(gap ~ age, rs); ABs[k+1,] <- coef(f2)
  d$pred_corr_s[te] <- d$pred_raw[te] - (coef(f2)[2]*age[te] + coef(f2)[1])
}
cat("dbias = a*x + b, fitted in each fold's reference set:\n"); print(round(AB,3))
cat(sprintf("Sensitivity (reference: survivors >=10 y, n~%d/fold): a %.3f +/- %.3f, b %.2f +/- %.2f; r between the two corrections %.5f\n",
    round(sum(d$t>=120)*0.8), mean(ABs[,2]), sd(ABs[,2]), mean(ABs[,1]), sd(ABs[,1]), cor(d$pred_corr, d$pred_corr_s)))
d$gap_raw <- d$pred_raw - age; d$gap <- d$pred_corr - age

cat("\n=========== 3. Prediction performance (mean +/- SD over the five folds) ===========\n")
met <- function(p, i){
  e <- p[i]-age[i]
  c(MAE=mean(abs(e)), CS5=100*mean(abs(e)<=5), r=cor(p[i],age[i]), rs=cor(p[i],age[i],method="spearman"), gap_vs_age=cor(e,age[i]))
}
tab <- function(p){
  m <- sapply(K, function(k) met(p, d$fold==k))
  sprintf("MAE %.3f +/- %.3f  CS5 %.1f +/- %.1f%%  r %.3f +/- %.3f  rs %.3f +/- %.3f  gap~age r %+.3f",
          m[1,] |> mean(), sd(m[1,]), mean(m[2,]), sd(m[2,]), mean(m[3,]), sd(m[3,]), mean(m[4,]), sd(m[4,]), mean(m[5,]))
}
for(nm in names(cfg)) cat(sprintf("  %-36s %s\n", nm, tab(cfg[[nm]])))
cat(sprintf("  %-36s %s\n", "Ensemble (raw)", tab(d$pred_raw)))
cat(sprintf("  %-36s %s\n", "  sensitivity: equal weights (raw)", tab(d$pred_eq)))
cat(sprintf("  %-36s %s\n", "Ensemble (bias-corrected)", tab(d$pred_corr)))
cat(sprintf("  %-36s %s\n", "  sensitivity: healthy-reference corr.", tab(d$pred_corr_s)))
cat("  MAE by age group (before / after correction):\n")
d$ab <- cut(age, c(59,64,69,74,79,84,120), labels=c("60-64","65-69","70-74","75-79","80-84","85+"))
for(a in levels(d$ab)){
  i <- d$ab==a
  cat(sprintf("    %-6s n=%4d  %.2f / %.2f   mean gap %+.2f / %+.2f\n", a, sum(i), mean(abs(d$gap_raw[i])), mean(abs(d$gap[i])),
              mean(d$gap_raw[i]), mean(d$gap[i])))
}

cat("\n=========== 4. Corrected gap and all-cause mortality ===========\n")
d$gap_sd <- d$gap/sd(d$gap); cat(sprintf("  gap SD = %.2f years; n=%d, deaths %d\n", sd(d$gap), n, sum(d$e)))
run <- function(fm_add, lab, dat=d){
  vars <- all.vars(as.formula(paste("~",fm_add))); s <- dat[complete.cases(dat[,c("t","e","gap",vars)]),]
  m1 <- coxph(as.formula(paste("Surv(t,e)~gap+",fm_add)), s); m0 <- coxph(as.formula(paste("Surv(t,e)~",fm_add)), s)
  m2 <- coxph(as.formula(paste("Surv(t,e)~gap_sd+",fm_add)), s)
  c1 <- summary(m1)$conf.int["gap",]; c2 <- summary(m2)$conf.int["gap_sd",]
  cat(sprintf("  %-52s n=%4d  HR per year %.3f (%.3f-%.3f)  HR per SD %.3f (%.3f-%.3f)  P=%.1e  dC %+.4f\n", lab, nrow(s),
              c1[1],c1[3],c1[4], c2[1],c2[3],c2[4], summary(m1)$coef["gap",5], concordance(m1)$concordance-concordance(m0)$concordance))
  invisible(m1)
}
m1 <- run("age+factor(sex)+factor(race)", "Model 1: age + sex + race")
run("age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi", "Model 2: + BMI, smoking, diabetes, arthritis, stroke, MI")
run("age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi+fnbmd", "Model 3: + femoral neck BMD")
cat(sprintf("  Schoenfeld PH test (model 1, gap term) P = %.3f\n", cox.zph(m1)$table["gap","p"]))
cat("  Comparison: HR per year of uncorrected gap (model 1)")
s <- summary(coxph(Surv(t,e)~gap_raw+age+factor(sex)+factor(race), d))$conf.int["gap_raw",]
cat(sprintf(" %.3f (%.3f-%.3f)\n", s[1],s[3],s[4]))
d$q <- cut(d$gap, quantile(d$gap, 0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))
mq <- coxph(Surv(t,e)~q+age+factor(sex)+factor(race), d); sq <- summary(mq)$conf.int
cat("  Quartiles (Q1 reference, model 1):")
for(k in 1:3) cat(sprintf(" Q%d %.3f (%.3f-%.3f)", k+1, sq[k,1],sq[k,3],sq[k,4]))
cat(sprintf("   log-rank P = %.1e\n", pchisq(survdiff(Surv(t,e)~q, d)$chisq, 3, lower.tail=FALSE)))
cat("  Mean corrected gap by age group:")
for(a in levels(d$ab)) cat(sprintf(" %s %+.2f", a, mean(d$gap[d$ab==a])))
cat("\n")

out <- d[,c("SEQN","fold","age","pred_raw","pred_corr","gap_raw","gap","q")]; names(out)[7] <- "gap_corr"
write.csv(out, file.path(res,"ksaa_oof.csv"), row.names=FALSE)
cat("-> ksaa_oof.csv\n")
