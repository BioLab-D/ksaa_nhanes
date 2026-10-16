# KSAA:五个骨干的折外预测按折内权重集成 + 训练折年龄偏倚校正;输出 results/ksaa_oof.csv 与性能 / 下游 Cox 日志
# (日志中 "Ensemble (raw)" / "Ensemble (bias-corrected)" 两行由 tableS1_performance.R 解析)
.fresh_rng <- function() { RNGkind("default", "default", "default")   # 每段开始时恢复新 R 进程的随机数状态
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }

# ======== 膝龄最终模型,按 Mammo-AGE 流程(Yan et al., Nat Commun 2025):五个事先指定骨干的五折 OOF 折内加权集成 + 偏倚校正 ========
.fresh_rng(); local({
# 膝龄最终模型,按 Mammo-AGE 流程(Yan et al., Nat Commun 2025):五个事先指定骨干的五折 OOF 折内加权集成 + 偏倚校正
#  权重 ∝ 1/训练折内部验证 MAE(折内归一),只用于该折留出者
#  偏倚校正(Smith 2019):Δbias = ŷ − x 对 x 回归,ŷ_corr = ŷ − (a·x + b),gap = ŷ_corr − x;系数取自留出折以外四折
#    (本队列无健康子集;敏感性分析用随访 ≥10 年存活者作参考)
#  输出 OOF CSV、性能与下游 Cox(模型 3 加股骨颈 BMD,对应 Mammo-AGE 的乳腺密度)日志
suppressMessages(library(survival))
res <- here::here("results"); S <- file.path(res,"server_knee_yolo_20260906")
sp  <- here::here("results", "work")
V2 <- Sys.getenv("V2","0")=="1"   # V2=1 读膝适配版 results/predictions/(4,875 人)
d <- if(V2){ o <- read.csv(file.path(res,"predictions/m4_resnet18_oof.csv")); sl0 <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(o$SEQN, sl0$SEQN)
  stopifnot(!anyNA(i), all(o$age==sl0$age[i]))
  data.frame(SEQN=o$SEQN, age=o$age, t=sl0$PERMTH_EXM[i], e=sl0$MORTSTAT[i], sex=sl0$sex[i], fold=o$fold)
} else read.csv(file.path(res,"mammo/mammo_ensemble_oof.csv"))[,c("SEQN","age","t","e","sex","fold")]
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN); d$race <- sl$race[i]; d$bmi <- sl$bmi[i]; d$fnbmd <- sl$fnbmd[i]
V <- read.csv(file.path(sp,"validate_vars.csv")); j <- match(d$SEQN, V$SEQN); q8 <- function(x) ifelse(x %in% c(8,9), NA, x)
d$dm <- as.integer(q8(V$HAD1[j])==1); d$arth <- as.integer(q8(V$HAC1A[j])==1); d$stroke <- as.integer(q8(V$HAC1D[j])==1)
d$mi <- as.integer(q8(V$HAF10[j])==1); d$smk <- as.integer(q8(V$HAR1[j])==1)
al <- function(f){ x <- read.csv(f); x$pred_age[match(d$SEQN, x$SEQN)] }
# 五个事先指定骨干各训练一次(种子 20260831),不按 OOF 表现事后挑选
MV <- Sys.getenv("MV","0")=="1"   # MV=1 读官方多视图 448 版
if(V2){ v2f <- function(a, s="oof") file.path(res, "predictions", sprintf("m4_%s_%s.csv", a, s))
  A5 <- c("ResNet-18"="resnet18", "ResNet-50"="resnet50", "EfficientNet-B0"="efficientnetb0", "ConvNeXt-Tiny"="convnexttiny", "DenseNet-121"="densenet121")
  cfg <- lapply(A5, function(a) al(v2f(a)))   # 五骨干同一配方
  FM <- sapply(A5, function(a){ f <- read.csv(v2f(a, "folds")); stopifnot(identical(f$fold, 0:4)); f$inner_val_mae })
  cat("V2=1:膝适配版 v2(m4 队列,N = 4);权重用训练折内部验证 MAE(v2/m4_*_folds.csv)\n")
} else if(MV){ MVd <- file.path(res,"mammo_mv"); mvf <- function(a, lr=FALSE) file.path(MVd, sprintf("off_448%s_%s_oof.csv", if(lr) "lr" else "", a))
  cfg <- list("ResNet-18"=al(mvf("resnet18")), "ResNet-50"=al(mvf("resnet50")), "EfficientNet-B0"=al(mvf("efficientnet_b0")),
              "ConvNeXt-Tiny"=al(mvf("convnext_tiny", lr=TRUE)), "DenseNet-121"=al(mvf("densenet121")))   # ConvNeXt 用按内部验证调低 lr 的重跑版
  FM <- sapply(sub("_oof","_foldmae", c(mvf("resnet18"),mvf("resnet50"),mvf("efficientnet_b0"),mvf("convnext_tiny",TRUE),mvf("densenet121"))), function(f) read.csv(f)$mae)
  cat("MV=1:官方多视图 448 版;权重用训练折内 10% 内部验证 MAE(Mammo-AGE 原义)\n")
} else cfg <- list(
  "ResNet-18"       = al(file.path(res,"mammo/mammo_resnet18_oof.csv")),
  "ResNet-50"       = al(file.path(res,"mammo/mammo_resnet50_oof.csv")),
  "EfficientNet-B0" = al(file.path(res,"mammo/mammo_efficientnet_b0_oof.csv")),
  "ConvNeXt-Tiny"   = al(file.path(res,"mammo/mammo_convnext_tiny_oof.csv")),
  "DenseNet-121"    = al(file.path(res,"mammo/mammo_densenet121_oof.csv")))
M <- do.call(cbind, cfg); stopifnot(!anyNA(M)); n <- nrow(d); K <- 0:4; age <- d$age

cat("=========== 1. 折内权重(∝ 1/内部验证 MAE,折内归一;第 k 折留出者只用第 k 折的权重)===========\n")
# 防泄漏:第 k 折只用 W[,k](不用五折平均,其内部验证集含留出者);五折平均 w 只作报告;等权集成作敏感性。
# MV=0 路径用其余四折的 OOF MAE 定权。
W <- if(MV || V2) apply(FM, 1, function(m) (1/m)/sum(1/m)) else sapply(K, function(k){ i <- d$fold!=k; m <- apply(M[i,],2,function(p) mean(abs(p-age[i]))); (1/m)/sum(1/m) })
if(MV || V2) rownames(W) <- names(cfg)
colnames(W) <- paste0("fold",K); w <- rowMeans(W)
print(round(cbind(W, mean=w),4))
d$pred_raw <- NA_real_; for(k in K){ i <- d$fold==k; d$pred_raw[i] <- as.numeric(M[i,,drop=FALSE] %*% W[,k+1]) }
d$pred_eq <- rowMeans(M)

cat("\n=========== 2. 偏倚校正(参考集 = 留出折以外四折)===========\n")
d$pred_corr <- d$pred_corr_s <- NA; AB <- matrix(NA,5,2,dimnames=list(paste0("fold",K),c("b(intercept)","a(slope)"))); ABs <- AB
for(k in K){ te <- d$fold==k; tr <- !te; r <- d[tr,]; r$gap <- r$pred_raw - r$age
  f <- lm(gap ~ age, r); AB[k+1,] <- coef(f); d$pred_corr[te] <- d$pred_raw[te] - (coef(f)[2]*age[te] + coef(f)[1])
  rs <- r[r$t >= 120,]; f2 <- lm(gap ~ age, rs); ABs[k+1,] <- coef(f2); d$pred_corr_s[te] <- d$pred_raw[te] - (coef(f2)[2]*age[te] + coef(f2)[1]) }
cat("Δbias = a·x + b,各折参考集拟合:\n"); print(round(AB,3))
cat(sprintf("敏感性(参考集限随访≥10 年存活者,n≈%d/折):a %.3f ± %.3f,b %.2f ± %.2f;两种校正后年龄相关 %.5f\n",
    round(sum(d$t>=120)*0.8), mean(ABs[,2]), sd(ABs[,2]), mean(ABs[,1]), sd(ABs[,1]), cor(d$pred_corr, d$pred_corr_s)))
d$gap_raw <- d$pred_raw - age; d$gap <- d$pred_corr - age

cat("\n=========== 3. 预测性能(五折均值 ± SD;Mammo-AGE Table 2 口径)===========\n")
met <- function(p, i){ e <- p[i]-age[i]; c(MAE=mean(abs(e)), CS5=100*mean(abs(e)<=5), r=cor(p[i],age[i]), rs=cor(p[i],age[i],method="spearman"), gap_vs_age=cor(e,age[i])) }
tab <- function(p){ m <- sapply(K, function(k) met(p, d$fold==k)); sprintf("MAE %.3f ± %.3f  CS5 %.1f ± %.1f%%  r %.3f ± %.3f  rs %.3f ± %.3f  gap~age r %+.3f", m[1,] |> mean(), sd(m[1,]), mean(m[2,]), sd(m[2,]), mean(m[3,]), sd(m[3,]), mean(m[4,]), sd(m[4,]), mean(m[5,])) }
for(nm in names(cfg)) cat(sprintf("  %-36s %s\n", nm, tab(cfg[[nm]])))
cat(sprintf("  %-36s %s\n", "Ensemble (raw)", tab(d$pred_raw)))
cat(sprintf("  %-36s %s\n", "  sensitivity: equal weights (raw)", tab(d$pred_eq)))
cat(sprintf("  %-36s %s\n", "Ensemble (bias-corrected)", tab(d$pred_corr)))
cat(sprintf("  %-36s %s\n", "  sensitivity: healthy-reference corr.", tab(d$pred_corr_s)))
cat("  按年龄段 MAE(校正前 / 后):\n"); d$ab <- cut(age, c(59,64,69,74,79,84,120), labels=c("60-64","65-69","70-74","75-79","80-84","85+"))
for(a in levels(d$ab)){ i <- d$ab==a; cat(sprintf("    %-6s n=%4d  %.2f / %.2f   gap 均值 %+.2f / %+.2f\n", a, sum(i), mean(abs(d$gap_raw[i])), mean(abs(d$gap[i])), mean(d$gap_raw[i]), mean(d$gap[i]))) }

cat("\n=========== 4. 下游:校正后 gap → 全因死亡 ===========\n")
d$gap_sd <- d$gap/sd(d$gap); cat(sprintf("  gap SD = %.2f 岁;n=%d,死亡 %d\n", sd(d$gap), n, sum(d$e)))
run <- function(fm_add, lab, dat=d){ vars <- all.vars(as.formula(paste("~",fm_add))); s <- dat[complete.cases(dat[,c("t","e","gap",vars)]),]
  m1 <- coxph(as.formula(paste("Surv(t,e)~gap+",fm_add)), s); m0 <- coxph(as.formula(paste("Surv(t,e)~",fm_add)), s)
  m2 <- coxph(as.formula(paste("Surv(t,e)~gap_sd+",fm_add)), s); c1 <- summary(m1)$conf.int["gap",]; c2 <- summary(m2)$conf.int["gap_sd",]
  cat(sprintf("  %-52s n=%4d  每年 HR %.3f (%.3f-%.3f)  每 SD HR %.3f (%.3f-%.3f)  P=%.1e  ΔC %+.4f\n", lab, nrow(s), c1[1],c1[3],c1[4], c2[1],c2[3],c2[4], summary(m1)$coef["gap",5], concordance(m1)$concordance-concordance(m0)$concordance)); invisible(m1) }
m1 <- run("age+factor(sex)+factor(race)", "模型1 实龄+性别+种族")
run("age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi", "模型2 +BMI+吸烟+糖尿病+关节炎+卒中+心梗")
run("age+factor(sex)+factor(race)+bmi+smk+dm+arth+stroke+mi+fnbmd", "模型3 +股骨颈 BMD(模态混杂)")
cat(sprintf("  Schoenfeld PH 检验(模型1,gap 项)P = %.3f\n", cox.zph(m1)$table["gap","p"]))
cat("  对照:未校正 gap 每年 HR(同模型1)"); s <- summary(coxph(Surv(t,e)~gap_raw+age+factor(sex)+factor(race), d))$conf.int["gap_raw",]; cat(sprintf(" %.3f (%.3f-%.3f)\n", s[1],s[3],s[4]))
d$q <- cut(d$gap, quantile(d$gap, 0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))
mq <- coxph(Surv(t,e)~q+age+factor(sex)+factor(race), d); sq <- summary(mq)$conf.int
cat("  四分位(Q1 参照,模型1):"); for(k in 1:3) cat(sprintf(" Q%d %.3f (%.3f-%.3f)", k+1, sq[k,1],sq[k,3],sq[k,4])); cat(sprintf("   log-rank P = %.1e\n", pchisq(survdiff(Surv(t,e)~q, d)$chisq, 3, lower.tail=FALSE)))
cat("  校正后 gap 按年龄段的均值(独立性检查,ANCOVA 思路):"); for(a in levels(d$ab)) cat(sprintf(" %s %+.2f", a, mean(d$gap[d$ab==a]))); cat("\n")
if (file.exists(file.path(res, "age_20260831_oof.csv"))) {   # [release] legacy v1-model comparison (diagnostic only; inputs not shipped)
main <- rowMeans(sapply(c(20260831,11,22), function(s) al(file.path(res,sprintf("age_%d_oof.csv",s))))); acc_main <- resid(lm(main~age, na.action=na.exclude)); cat(sprintf("  与旧主模型(r18 三种子)残差法加速的相关 %.3f;其每 SD HR:", cor(d$gap, acc_main, use="complete.obs")))   # 新增 3 人无旧模型预测
s <- summary(coxph(Surv(t,e)~z+age+factor(sex)+factor(race), transform(d, z=as.numeric(scale(acc_main)))))$conf.int[1,]; cat(sprintf(" %.3f (%.3f-%.3f)\n", s[1],s[3],s[4]))
}

out <- d[,c("SEQN","fold","age","pred_raw","pred_corr","gap_raw","gap","q")]; names(out)[7] <- "gap_corr"
write.csv(out, file.path(res,"ksaa_oof.csv"), row.names=FALSE)
cat("-> ksaa_oof.csv\n")
})
