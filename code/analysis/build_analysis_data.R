# 构建统一分析数据集 analysis_data.rds(第 1 段),并把本队列重算的 PhenoAge 写回 work/oof_pheno.rds(第 2 段)。
.fresh_rng <- function() { RNGkind("default", "default", "default")   # 每段开始时恢复新 R 进程的随机数状态
  if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv()) }

# ======== 构建统一分析数据集(analysis_data.rds):主表型队列 4,875 人 + 结局、PhenoAge、协变量、K–L、BMD、功能与左右膝读出 ========
.fresh_rng(); local({
# 构建统一分析数据集(analysis_data.rds):主表型队列 4,875 人 + 结局、PhenoAge、协变量、K–L、BMD、功能与左右膝读出
# PhenoAge 在本队列上重算。
# 由 run_all.sh 运行
suppressMessages(library(duckdb))
res <- here::here("results"); sp <- file.path(res, "work"); out <- file.path(res, "analysis")
dir.create(out, showWarnings = FALSE)
q <- function(x){ x <- suppressWarnings(as.numeric(x)); for(k in 1:6){ s8 <- as.numeric(strrep("8", k)); lo <- x[!is.na(x) & x < s8]
  if(any(x == s8, na.rm = TRUE) && length(lo) && s8 >= max(lo) * 1.5) x[x == s8] <- NA }; x }   # 仅当全 8 远高于其余观测时视为哨兵
na_at <- function(x, codes) { x <- as.numeric(x); x[x %in% codes] <- NA; x }                      # 按变量字段宽度显式指定哨兵
yn <- function(x) ifelse(x == 1, 1L, ifelse(x == 2, 0L, NA_integer_))                               # 1 = 是,2 = 否,8/9 = 缺失
prev <- function(v) sprintf("%-10s n %4d  缺失 %4d  阳性 %5.1f%%", v, sum(!is.na(d[[v]])), sum(is.na(d[[v]])), 100 * mean(d[[v]], na.rm = TRUE))

# ---------------- 队列与结局 ----------------
d <- read.csv(file.path(res, "ksaa_oof.csv"))
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN); stopifnot(!anyNA(i), all(d$age == sl$age[i]))
for(v in c("sex", "race", "phase", "pir", "bmi", "fnbmd", "WTPFEX6", "WTPFEX2", "SDPPSU6", "SDPSTRA6", "SDPPSU2", "SDPSTRA2")) d[[v]] <- sl[[v]][i]
d$t <- sl$PERMTH_EXM[i]; d$e <- sl$MORTSTAT[i]; d$ucod <- sl$UCOD_LEADING[i]; d$gap <- d$gap_corr   # 时间零点 = MEC 体检(PERMTH_EXM)
d$e_cvd <- as.integer(d$e == 1 & d$ucod %in% c(1, 5)); d$e_cancer <- as.integer(d$e == 1 & d$ucod %in% 2)   # NCHS 分组:1 心脏病、5 脑血管、2 恶性肿瘤;无死因的 5 名死者两者均记 0
cat(sprintf("队列 n %d;死亡 %d(无死因 %d);CVD 死亡 %d;癌症死亡 %d\n", nrow(d), sum(d$e), sum(d$e == 1 & is.na(d$ucod)), sum(d$e_cvd), sum(d$e_cancer)))
cat(sprintf("surv_labels 的 bmi 范围 %s、fnbmd 范围 %s(确认已清洗)\n", paste(range(d$bmi, na.rm = TRUE), collapse = "–"), paste(round(range(d$fnbmd, na.rm = TRUE), 3), collapse = "–")))

# ---------------- NHANES III 原始变量 ----------------
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
  141.50225 + (log(0.00553 * k) + xb) / 0.09165   # log(1 − M) 解析化简为 −exp(xb)·k,避免浮点溢出
}
d$phenoage <- pheno(d); d$pa_accel <- resid(lm(phenoage ~ age, d, na.action = na.exclude))
# PhenoAge 加速按 1% / 99% 分位缩尾后用于一切分析(右尾极长);原始值 pa_accel_raw 仅供描述性表格。
d$pa_accel_raw <- d$pa_accel
PA_W <- quantile(d$pa_accel_raw, c(.01, .99), na.rm = TRUE)
d$pa_accel <- pmin(pmax(d$pa_accel_raw, PA_W[1]), PA_W[2])
cat(sprintf("\nPhenoAge 加速缩尾:1%% / 99%% 分位 %.2f / %.2f 年;截到下限 %d 人、上限 %d 人;SD %.2f -> %.2f\n", PA_W[1], PA_W[2],
    sum(d$pa_accel_raw < PA_W[1], na.rm = TRUE), sum(d$pa_accel_raw > PA_W[2], na.rm = TRUE), sd(d$pa_accel_raw, na.rm = TRUE), sd(d$pa_accel, na.rm = TRUE)))
old <- readRDS(file.path(sp, "oof_pheno.rds")); j <- match(d$SEQN, old$SEQN)
cat(sprintf("\nPhenoAge 可算 %d / %d;与旧 oof_pheno.rds 共有者最大差 %.2e;新增可算者 %d(SEQN %s)\n", sum(!is.na(d$phenoage)), nrow(d),
    max(abs(d$phenoage - old$phenoage[j]), na.rm = TRUE), sum(!is.na(d$phenoage) & is.na(j)), paste(d$SEQN[!is.na(d$phenoage) & is.na(j)], collapse = ", ")))

# ---------------- 社会人口与生活方式 ----------------
ed <- na_at(AD$HFA8R, c(88, 99)); d$edu3 <- cut(ed, c(-Inf, 11, 12, Inf), labels = c("<12", "12", ">12"))
d$smk3 <- factor(ifelse(AD$HAR1 == 2, "never", ifelse(AD$HAR1 == 1 & AD$HAR3 == 2, "former", ifelse(AD$HAR1 == 1 & AD$HAR3 == 1, "current", NA))), levels = c("never", "former", "current"))
# 饮酒:一生 <12 杯 = never;一生 ≥12 杯但过去一年 <12 杯 = former/light(含当前少量饮酒者);过去一年 ≥12 杯 = current
d$alc3 <- factor(ifelse(EX$MAPE1 == 2, "never", ifelse(EX$MAPE1 == 1 & EX$MAPE2 == 2, "former_light", ifelse(EX$MAPE1 == 1 & EX$MAPE2 == 1, "current", NA))), levels = c("never", "former_light", "current"))
# 休闲体力活动:过去一个月各项次数 × MET;中等 3–5.9、高强度 ≥6。活跃 = 中等 ≥20 次/月(≈5 次/周)或高强度 ≥12 次/月(≈3 次/周);不活跃 = 零次
bad <- rep(FALSE, nrow(d)); mod <- vig <- tot <- mets <- rep(0, nrow(d)); nomet <- 0
for(p in PA){ f <- as.numeric(AD[[p[1]]]); m <- as.numeric(AD[[p[2]]]); gate <- if(is.na(p[3])) rep(1, nrow(d)) else as.numeric(AD[[p[3]]])
  bad <- bad | (!is.na(f) & f >= 7777) | (!is.na(gate) & gate %in% c(8, 9)); f[is.na(f) | f >= 7777] <- 0
  bad <- bad | (f > 0 & !is.na(m) & m >= 88); m[!is.na(m) & m >= 88] <- NA   # HAT19MET 有 8888 哨兵(3 人):强度未知,体力活动记缺失
  nomet <- nomet + sum(f > 0 & is.na(m)); m[is.na(m)] <- 3
  tot <- tot + f; mets <- mets + f * m; mod <- mod + f * (m >= 3 & m < 6); vig <- vig + f * (m >= 6) }
d$pa_sessions <- ifelse(bad, NA, tot); d$pa_metsess <- ifelse(bad, NA, mets)
d$pa3 <- factor(ifelse(bad, NA, ifelse(tot == 0, "inactive", ifelse(mod >= 20 | vig >= 12, "active", "insufficient"))), levels = c("active", "insufficient", "inactive"))
cat(sprintf("体力活动:频次或门控为哨兵者 %d;有频次但无 MET(按 3 计)%d 项次;月次数 P50/P90/max %s\n", sum(bad), nomet, paste(quantile(d$pa_sessions, c(.5, .9, 1), na.rm = TRUE), collapse = "/")))

# ---------------- 心血管代谢与病史 ----------------
d$sbp <- na_at(EX$PEPMNK1R, 888); d$dbp <- na_at(EX$PEPMNK5R, 888); d$waist <- na_at(EX$BMPWAIST, 88888); d$fast_h <- na_at(LB$PHPFAST, 88888)
d$tc <- na_at(LB$TCP, 888); d$hdl <- na_at(LB$HDP, 888); d$tg <- na_at(LB$TGP, 8888); d$hba1c <- na_at(LB$GHP, 8888); d$fib <- na_at(LB$FBP, 8888)
or1 <- function(...) { m <- cbind(...); ifelse(rowSums(m == 1, na.rm = TRUE) > 0, 1L, ifelse(rowSums(is.na(m)) == 0, 0L, NA_integer_)) }
d$htn <- or1(yn(AD$HAE2), as.integer(d$sbp >= 140), as.integer(d$dbp >= 90))                       # 自报医生诊断或测量 ≥140/90
d$dm <- or1(yn(AD$HAD1), as.integer(d$hba1c >= 6.5))                                               # 自报或 HbA1c ≥6.5%
d$dm[is.na(d$dm) & AD$HAD1 == 2 & is.na(d$hba1c)] <- 0L                                            # 自报否且 HbA1c 缺失:按自报记 0,不当缺失
bpmed <- ifelse(AD$HAE2 == 2 & is.na(AD$HAE5A), 0L, yn(AD$HAE5A))                                   # HAE5A 只问自报高血压者;未被问到 = 未服药
d$dyslip <- or1(as.integer(d$tc >= 240), as.integer(d$hdl < 40))
d$cvd_hist <- or1(yn(AD$HAF10), yn(AD$HAC1C), yn(AD$HAC1D))                                        # 心梗、心衰、卒中
d$cancer_hist <- yn(AD$HAC1O)                                                                       # 非皮肤癌("other cancer")
d$arth <- yn(AD$HAC1A); d$obese <- as.integer(d$bmi >= 30)
# eGFR:NHANES III 肌酐先按 Selvin 2007 校准到标准化方法(0.960 × Scr − 0.184)[verify],再用 CKD-EPI 2021(无种族项)
scr <- 0.960 * d$CEP - 0.184; fem <- d$sex == 2; kap <- ifelse(fem, 0.7, 0.9); alp <- ifelse(fem, -0.241, -0.302)
d$egfr <- 142 * pmin(scr / kap, 1)^alp * pmax(scr / kap, 1)^-1.200 * 0.9938^d$age * ifelse(fem, 1.012, 1)
# 代谢综合征(ATP III 五项 ≥3):禁食 ≥8 h 用 TG ≥150 与空腹血糖 ≥100;非禁食用 TG ≥175、HbA1c ≥5.7;糖尿病与降压药计入对应项
fasting <- !is.na(d$fast_h) & d$fast_h >= 8
cr <- cbind(waist = as.integer(ifelse(fem, d$waist > 88, d$waist > 102)),
            tg    = as.integer(ifelse(fasting, d$tg >= 150, d$tg >= 175)),
            hdl   = as.integer(ifelse(fem, d$hdl < 50, d$hdl < 40)),
            bp    = or1(as.integer(d$sbp >= 130), as.integer(d$dbp >= 85), bpmed),
            glu   = or1(d$dm, as.integer(ifelse(fasting, d$SGP >= 100, d$hba1c >= 5.7))))
nmet <- rowSums(cr == 1, na.rm = TRUE); nmis <- rowSums(is.na(cr)); d$mets <- ifelse(nmet >= 3, 1L, ifelse(nmet + nmis < 3, 0L, NA_integer_))
# 胰岛素抵抗:空腹 ≥8 h、未用胰岛素者;HOMA-IR = G1P (mg/dL) × I1P (μU/mL) / 405;哨兵 G1P 88888、I1P 888888
d$fpg <- na_at(LB$G1P, 88888); d$insulin <- na_at(LB$I1P, 888888); ins <- (yn(AD$HAD6) %in% 1L) | (LB$PHPINSU %in% 1)
d$lhoma <- ifelse(fasting & !ins, log(d$fpg * d$insulin / 405), NA); d$linsulin <- ifelse(fasting & !ins, log(d$insulin), NA)
cat(sprintf("HOMA-IR:禁食 ≥8 h、未用胰岛素且两项有效 %d 人(用胰岛素 %d);PHPINSU 取值 %s;HOMA-IR 中位数 %.2f\n", sum(!is.na(d$lhoma)), sum(ins),
    paste(names(table(LB$PHPINSU, useNA = "ifany")), table(LB$PHPINSU, useNA = "ifany"), sep = ":", collapse = " "), exp(median(d$lhoma, na.rm = TRUE))))
cat(sprintf("禁食 ≥8 h:%d 人(%.0f%%);代谢综合征各项阳性率 %s\n", sum(fasting), 100 * mean(fasting), paste(sprintf("%s %.0f%%", colnames(cr), 100 * colMeans(cr, na.rm = TRUE)), collapse = "  ")))

# ---------------- 肌骨与炎症 ----------------
kl <- read.csv(file.path(res, "work", "kl_per_knee.csv")); cat(sprintf("K–L 取值:%s\n", paste(names(table(kl$kl)), table(kl$kl), sep = ":", collapse = "  ")))
klp <- aggregate(kl ~ seqn, kl[kl$kl %in% 0:4, ], max); d$klmax <- klp$kl[match(d$SEQN, klp$seqn)]
d$oa <- as.integer(d$klmax >= 2); d$severe_oa <- as.integer(d$klmax >= 3)
KP <- read.csv(file.path(sp, "knee_pain.csv")); d$pain <- yn(q(KP$HAG22[match(d$SEQN, KP$SEQN)]))
DX <- read.csv(file.path(sp, "dxa_all.csv")); j <- match(d$SEQN, DX$SEQN); d$fn_bmd <- q(DX$BDPFNBMD)[j]; d$tot_bmd <- (q(DX$BDPTOBMC) / q(DX$BDPTOARE))[j]
P <- read.csv(file.path(sp, "physfun.csv")); j <- match(d$SEQN, P$SEQN)
d$gait <- 2.4384 / pmin(q(P$PFPTWLKA), q(P$PFPTWLKB), na.rm = TRUE)[j]; d$gait[!is.finite(d$gait)] <- NA
d$tandem <- na_at(EX$PFPTTAND, 8888)
d$lcrp <- log(d$CRP); d$wbc <- d$WCP; d$glucose <- d$SGP
# 骨折、跌倒、功能受限、贫困(按四分位的基线表用)
qn <- function(x){ x <- suppressWarnings(as.numeric(x)); for(k in 1:6) for(dg in c("8", "9")){ s <- as.numeric(strrep(dg, k)); lo <- x[!is.na(x) & x < s]
  if(any(x == s, na.rm = TRUE) && length(lo) && s >= max(lo) * 1.5) x[x == s] <- NA }; x }
FR <- read.csv(file.path(sp, "fracture.csv")); jf <- match(d$SEQN, FR$SEQN); jp <- match(d$SEQN, P$SEQN)
d$fx_hip <- yn(FR$HAG5A[jf]); fxm <- cbind(d$fx_hip, yn(FR$HAG5B[jf]), yn(FR$HAG5C[jf]))
d$fx_any <- as.integer(rowSums(fxm, na.rm = TRUE) > 0); d$fx_any[rowSums(!is.na(fxm)) == 0] <- NA
d$fall <- as.integer(qn(FR$HAG14)[jf] > 0)
hh <- sapply(c("HAH1", "HAH2", "HAH8", "HAH9"), function(v) { x <- P[[v]][jp]; as.integer(ifelse(x %in% c(8, 9), NA, x) >= 2) })
d$anylimit <- as.integer(rowSums(hh, na.rm = TRUE) > 0); d$anylimit[rowSums(!is.na(hh)) == 0] <- NA
d$poor <- as.integer(d$pir < 1)

# ---------------- 左右膝读出:同一双膝模型的膝级局部头,按主表型同一折内权重与折外偏倚校正 ----------------
# 注意:局部头带有对侧膝信息(交叉注意力),gap_L / gap_R 不是独立的单膝估计
A5 <- c(resnet18 = "resnet18", resnet50 = "resnet50", efficientnetb0 = "efficientnetb0", convnexttiny = "convnexttiny", densenet121 = "densenet121")
FM <- sapply(A5, function(a) read.csv(file.path(res, "predictions", sprintf("m4_%s_folds.csv", a)))$inner_val_mae); W <- apply(FM, 1, function(m) (1/m) / sum(1/m))
for(sd_ in c("L", "R")){
  M <- sapply(A5, function(a){ x <- read.csv(file.path(res, "predictions", sprintf("m4_%s_oof.csv", a))); x[[paste0("pred_", sd_)]][match(d$SEQN, x$SEQN)] }); stopifnot(!anyNA(M))
  raw <- rep(NA_real_, nrow(d)); for(k in 0:4){ i <- d$fold == k; raw[i] <- as.numeric(M[i, , drop = FALSE] %*% W[, k + 1]) }
  g <- rep(NA_real_, nrow(d)); for(k in 0:4){ te <- d$fold == k; f <- lm(I(raw - age) ~ age, data.frame(raw = raw[!te], age = d$age[!te])); g[te] <- raw[te] - d$age[te] - predict(f, data.frame(age = d$age[te])) }
  d[[paste0("gap_", sd_)]] <- g }
# 核对:两膝读出的均值经同一流程应接近主表型
M <- sapply(A5, function(a){ x <- read.csv(file.path(res, "predictions", sprintf("m4_%s_oof.csv", a))); x$pred_age[match(d$SEQN, x$SEQN)] })
raw <- rep(NA_real_, nrow(d)); for(k in 0:4){ i <- d$fold == k; raw[i] <- as.numeric(M[i, , drop = FALSE] %*% W[, k + 1]) }
cat(sprintf("\n重建主表型原始预测 vs ksaa_oof.csv 最大差 %.2e;左/右膝 gap 与主 gap 相关 %.3f / %.3f;左右相关 %.3f\n",
    max(abs(raw - d$pred_raw)), cor(d$gap_L, d$gap), cor(d$gap_R, d$gap), cor(d$gap_L, d$gap_R)))

# ---------------- 核对:缺失、患病率、年龄方向 ----------------
cat("\n二分类变量:\n"); for(v in c("htn", "dm", "dyslip", "cvd_hist", "cancer_hist", "arth", "obese", "mets", "oa", "severe_oa", "pain")) cat("  ", prev(v), "\n")
cat("\n分类变量:\n"); for(v in c("edu3", "smk3", "alc3", "pa3")) cat(sprintf("  %-5s %s  缺失 %d\n", v, paste(names(table(d[[v]])), table(d[[v]]), sep = ":", collapse = "  "), sum(is.na(d[[v]]))))
cat("\n连续变量 (n, 均值, SD, min, max):\n")
for(v in c("sbp", "dbp", "tc", "hdl", "tg", "hba1c", "egfr", "waist", "fib", "lcrp", "wbc", "glucose", "fn_bmd", "tot_bmd", "gait", "tandem", "pa_metsess", "phenoage", "pa_accel_raw", "pa_accel"))
  cat(sprintf("  %-10s n %4d  %8.2f  %7.2f  %8.2f  %8.2f\n", v, sum(!is.na(d[[v]])), mean(d[[v]], na.rm = TRUE), sd(d[[v]], na.rm = TRUE), min(d[[v]], na.rm = TRUE), max(d[[v]], na.rm = TRUE)))
cat("\n年龄方向(每 10 岁 OR,应与常识一致):\n")
for(v in c("htn", "dm", "cvd_hist", "cancer_hist", "oa", "mets")) cat(sprintf("  %-12s %.2f\n", v, exp(10 * coef(glm(d[[v]] ~ d$age, family = binomial))[2])))
cat(sprintf("  不活跃 vs 其余 %.2f;当前饮酒 vs 其余 %.2f;当前吸烟 vs 其余 %.2f\n", exp(10 * coef(glm(I(pa3 == "inactive") ~ age, binomial, d))[2]),
    exp(10 * coef(glm(I(alc3 == "current") ~ age, binomial, d))[2]), exp(10 * coef(glm(I(smk3 == "current") ~ age, binomial, d))[2])))

saveRDS(d, file.path(out, "analysis_data.rds")); cat(sprintf("\n-> %s(%d × %d)\n", file.path(out, "analysis_data.rds"), nrow(d), ncol(d)))
})

# ======== 用本队列(4,875 人)上重算的 PhenoAge(build_analysis_data.R)覆盖 work/oof_pheno.rds,供下游按 SEQN 读取;原文件备份为 oof_pheno_v1cohort4872.rds。 ========
.fresh_rng(); local({
# 用本队列(4,875 人)上重算的 PhenoAge(build_analysis_data.R)覆盖 work/oof_pheno.rds,供下游按 SEQN 读取;原文件备份为 oof_pheno_v1cohort4872.rds。
# 下游只用 SEQN、phenoage、pa_accel(按名字取第一个含 "accel" 的列,须是 pa_accel);其余列按 SEQN 保留。
res <- here::here("results"); sp <- file.path(res, "work")
bk <- file.path(sp, "oof_pheno_v1cohort4872.rds"); old <- readRDS(file.path(sp, "oof_pheno.rds"))
if (!file.exists(bk)) { stopifnot(nrow(old) == 4872); saveRDS(old, bk); cat("备份 ->", bk, "\n") } else old <- readRDS(bk)
v <- readRDS(file.path(res, "analysis", "analysis_data.rds"))
# pa_accel 为 1/99 缩尾值(分析用),pa_accel_raw 为原始值(描述用)
new <- data.frame(SEQN = v$SEQN, phenoage = v$phenoage, pa_accel = v$pa_accel, zpa = as.numeric(scale(v$pa_accel)), pa_accel_raw = v$pa_accel_raw)
keep <- setdiff(names(old), c("phenoage", "pa_accel", "zpa", "phenoage_nc", "znc", "pa_accel_raw"))
new <- merge(new, old[, keep, drop = FALSE], by = "SEQN", all.x = TRUE)
stopifnot(nrow(new) == 4875, grep("accel", names(new), value = TRUE)[1] == "pa_accel")
i <- match(old$SEQN, new$SEQN)
cat(sprintf("新 oof_pheno.rds:%d 人,PhenoAge 非缺失 %d(旧 %d);共有者 phenoage 最大差 %.1e,pa_accel 最大差 %.3f(残差在新样本上重估)\n",
    nrow(new), sum(!is.na(new$phenoage)), sum(!is.na(old$phenoage)), max(abs(new$phenoage[i] - old$phenoage), na.rm = TRUE), max(abs(new$pa_accel[i] - old$pa_accel), na.rm = TRUE)))
saveRDS(new, file.path(sp, "oof_pheno.rds")); cat("->", file.path(sp, "oof_pheno.rds"), "\n")
})
