# EWAS 清洗与扫描。EWAS_STAGE=audit(发布流程所用):只出逐变量计数与清洗后表型,供 phewas_audit.R / phewas.R,不拟合模型。
# 默认 run 阶段:KSAA(gap,岁)~ 变量 + 年龄 + 性别 + 种族,不加权 OLS;连续变量按每 SD(偏态者先取对数),另报标准化 β。
#   分类变量(≥3 级)用整体 F 检验;所有变量另报偏 R²(加入该变量后残差平方和的相对下降),作跨类型可比的效应量。
# 缺失码:按 phewas_dict.R 从官方 codebook 读出的逐变量编码处理:
#   ① codebook 未列出的值 → 缺失(若干 *S 变量的 "Don't know" 在数据中是 9998 而非 codebook 的 9999);
#   ② 标签为 Blank but applicable / Don't know / Refused / Could not be assessed 等 → 缺失;
#   ③ 连续变量中的 7-fill "Less than 1 …" → 0.5 个单位("less than one half hour" → 0.25;"less than the lowest possible" → 区间下限的一半);
#   ④ 连续变量中 "Never" 记 0(频率题),但在"多久以前 / 几岁时"类题中 → 缺失;"Says can't do""Skinfold too large" 等 → 缺失。
# 纳入:DuckDB 数值列;三个文件的公共人口学表头只取 adult 一份;队列内非缺失 ≥ 1000;连续变量 ≥ 5 个不同取值;分类变量每级 ≥ 30(不足的级设缺失,剩余 ≥ 2 级)。
# 剔除:抽样权重 / 设计变量 / 行政与流程字段 / 牙齿逐面与牙周逐位点编码 / SI 单位重复 / 协变量本身 / codebook 核对失败者,另加人工清单(见 EXCL_MANUAL,附理由)。
# 多重比较:全部检验 BH-FDR;另报 Bonferroni 与 Li–Ji 有效检验数;FDR 显著者按 |r| > 0.7 贪心归并为独立信号。
# 输入 results/analysis/analysis_data.rds、results/analysis/phewas/ewas_{vars,codes,dict_check}.csv、nhanes_lab/data/nhanes.duckdb
# 输出 results/analysis/phewas/:audit 阶段 ewas_audit_counts.csv、ewas_clean.rds;run 阶段 ewas_candidates.csv、ewas_{results,domains,top}{_main,}.csv(_main 不含骨折、跌倒与身体功能)、ewas_cor.rds、ewas_clean_top.rds
# 由 run_all.sh 以 EWAS_STAGE=audit 运行;只出候选清单:EWAS_STAGE=list Rscript code/analysis/phewas_clean.R
suppressMessages({library(DBI); library(duckdb)})
res <- here::here("results"); out <- file.path(res, "analysis"); ew <- file.path(out, "phewas")
STAGE <- Sys.getenv("EWAS_STAGE", "run")
fp <- function(p) ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p))
d <- readRDS(file.path(out, "analysis_data.rds")); SDG <- sd(d$gap)
V <- read.csv(file.path(ew, "ewas_vars.csv")); C <- read.csv(file.path(ew, "ewas_codes.csv")); K <- read.csv(file.path(ew, "ewas_dict_check.csv"))
V$section[is.na(V$section)] <- ""; V$subsection[is.na(V$subsection)] <- ""; C$code_label[is.na(C$code_label)] <- ""
V$label <- ifelse(is.na(V$sas_label) | !nzchar(V$sas_label), V$index_label, V$sas_label)

# ---------- 1. 候选与剔除 ----------
V <- V[V$column_type %in% c("DOUBLE", "INTEGER", "TINYINT", "SMALLINT", "BIGINT"), ]
V <- V[!(V$file != "adult" & V$variable %in% V$variable[V$file == "adult"]), ]          # 公共表头只留 adult
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
ex(grepl("SI$", V$variable) & sub("SI$", "", V$variable) %in% V$variable, "SI-unit duplicate")      # 只有 SI 版本的(如 ALP)保留
ex(grepl("check item|interviewer|examiner|technician|recorder|comment|\\bflag\\b|reliab|language of|version|proxy|who (answered|completed)|reason|why (not|did not)|\\bstatus\\b|session|date of|time of day|day of (the )?week|today's date|\\bmonth of\\b|\\byear of (exam|interview)|serial|machine|calibration|\\bcode\\b",
         V$label, ignore.case = TRUE), "administrative field (label)")
bad <- K[with(K, !ok_all & (n_uncovered_values > 0.01 * n_rows | abs(sum_counts - n_rows) > 0.01 * n_rows)), ]
ex(paste(V$file, V$variable) %in% paste(bad$file, bad$variable), "codebook/data mismatch > 1%")
# 人工清单(逐条看过候选标签后补充;理由写明)
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
# 规则误杀、人工救回:城乡、定向力(今日日期 / 星期 / 地址)、访员观察到的行动状况、婚姻、因病改变饮食、无牙、医师评定健康、胰岛素使用、无食物原因
RESCUE <- c("DMPMETRO", "HAA4", "HAA5", "HAA6D", "HAM14", grep("^HAM15", V$variable, value = TRUE), grep("^HAY11", V$variable, value = TRUE),
            "HFA12", "HFF6A", "HFF6B", "HFF6C", "HFF6D", "DEPEDENT", "PEP13A", "PHPINSU")
V$excl_note[V$variable %in% RESCUE & V$excl != ""] <- paste("rescued from rule:", V$excl[V$variable %in% RESCUE & V$excl != ""])
V$excl[V$variable %in% RESCUE] <- ""

# ---------- 2. 队列数据 ----------
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

# ---------- 3. 逐变量清洗 ----------
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
  # 连续:有整数跨度 ≥ 4 的区间,或区间带小数(BMD、钙 mmol/L、腰臀比的跨度都 < 4);或未加标签的逐值枚举 ≥ 5 个
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
    # 分类变量:"Unable to do / Can't do" 是真实结局等级,保留;"Can't grade / could not be assessed / never found pulse" 是测不到,置缺失
    for (r in which(grepl(CAT_NA_RE, lab_single$code_label, ignore.case = TRUE))) {
      n_spc <- n_spc + sum(eq(x, lab_single$lo[r])); x[eq(x, lab_single$lo[r])] <- NA
    }
    lab_single <- lab_single[!grepl(CAT_NA_RE, lab_single$code_label, ignore.case = TRUE), ]
  }
  labs <- setNames(trimws(sub("\\s*\\(.*$", "", lab_single$code_label)), format(lab_single$lo, scientific = FALSE, trim = TRUE))   # 去掉跳题提示 "(HAS26)"
  list(x = x, cont = cont, n_uncovered = n_unc, n_special = n_spc, n_recoded = n_rec, labs = labs)
}
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3

# ---------- 4. 逐变量模型 ----------
cand <- V[V$excl == "", ]
base <- data.frame(gap = d$gap, age = d$age, sex = factor(d$sex), race = factor(d$race))
R <- list(); X <- list(); XR <- list(); FILT <- character(); AUD <- list(); LABS <- list(); CONT <- list()   # AUD:审核阶段的逐变量计数(EWAS_STAGE=audit)
for (i in seq_len(nrow(cand))) {
  v <- cand$variable[i]; f <- cand$file[i]; cc <- C[C$file == f & C$variable == v, ]
  cl <- clean_one(RAW[[v]], cc, cand$label[i]); x <- cl$x
  if (STAGE == "audit") {   # 只出计数,不拟合、不看效应量;Blank = 不适用/未问,Blank but applicable = 适用但缺失
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
    XR[[v]] <- cl$x; LABS[[v]] <- cl$labs; CONT[[v]] <- cl$cont      # 清洗后原值,供 phewas.R
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
    s$xi <- qnorm((rank(s$x, ties.method = "average") - 0.5) / nrow(s))                 # 敏感性:秩逆正态
  } else {
    tb <- table(s$x); keep <- names(tb)[tb >= 30]
    if (length(keep) < 2) { FILT[v] <- "fewer than 2 levels with >= 30"; next }
    s <- s[as.character(s$x) %in% keep, ]; n <- nrow(s)
    if (n < 1000) { FILT[v] <- sprintf("n %d < 1000 after dropping sparse levels", n); next }
    # 参照级:标签为 No / None / Never / No difficulty / Normal 等者优先,否则取最小编码(NHANES III 有序量表多以 1 = 无 / 最轻)
    lv <- as.character(sort(as.numeric(keep))); ll <- cl$labs[lv]
    neg <- which(grepl("^(no|none|never|normal|not)\\b|no difficulty|no findings|no pain", ll, ignore.case = TRUE))
    ref <- if (length(neg)) lv[neg[1]] else lv[1]
    s$x <- relevel(factor(as.character(s$x), levels = lv), ref = ref)
  }
  cov <- paste(c("age", if (length(unique(s$sex)) > 1) "sex", if (length(unique(s$race)) > 1) "race"), collapse = " + ")   # 单性别题(生殖健康)不调整性别
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
    # 多分类:整体 F 检验为主;另记与参照级差别最大的一级及其 95% CI(描述用,beta_y 等列留空以免与每 SD 混读)
    k <- which.max(abs(cf[xr, 1])); ci <- confint(m1)[xr[k], ]
    row$ref_level <- lab_of(levels(s$x)[1]); row$top_level <- lab_of(sub("^x", "", rownames(cf)[xr][k]))
    row$top_beta <- cf[xr[k], 1]; row$top_lo <- ci[1]; row$top_hi <- ci[2]
    row$contrast <- sprintf("%d levels; ref %s; largest %s %+.2f y", nlevels(s$x), row$ref_level, row$top_level, cf[xr[k], 1])
  }
  if (type == "continuous") { mi <- lm(as.formula(paste("gap ~ xi +", cov)), s); row$beta_int <- coef(mi)[["xi"]]; row$p_int <- summary(mi)$coefficients["xi", 4] }
  for (nm in c("top_beta", "top_lo", "top_hi")) if (is.null(row[[nm]])) row[[nm]] <- NA
  R[[v]] <- row
  # 供相关矩阵用的数值版本(行号即队列行号;多分类以"非参照级"二值化,只用于归并与有效检验数)
  z <- rep(NA_real_, nrow(d)); z[as.integer(rownames(s))] <- if (type == "continuous") s$x else as.numeric(s$x != levels(s$x)[1]); X[[v]] <- z
  XR[[v]] <- cl$x                                                                            # 清洗后、未变换的原值
}
if (STAGE == "audit") {
  A <- do.call(rbind, AUD); rownames(A) <- NULL
  A$n_applicable <- A$n_cohort - A$n_blank; A$validity_applicable <- A$n_valid / A$n_applicable; A$validity_cohort <- A$n_valid / A$n_cohort
  write.csv(A, file.path(ew, "ewas_audit_counts.csv"), row.names = FALSE)
  saveRDS(list(SEQN = d$SEQN, x = XR, labs = LABS, cont = CONT), file.path(ew, "ewas_clean.rds"))
  cat(sprintf("audit:%d 个规则内候选的计数已写出 ewas_audit_counts.csv(不含效应量与 P 值)\n", nrow(A))); quit(save = "no")
}
R <- do.call(rbind, R); rownames(R) <- NULL

cand_out <- V[, c("file", "variable", "label", "section", "subsection", "excl", "excl_note")]
cand_out$analysed <- cand_out$variable %in% R$variable & cand_out$excl == ""
cand_out$filter <- ifelse(cand_out$variable %in% names(FILT) & cand_out$excl == "", FILT[cand_out$variable], "")
write.csv(cand_out, file.path(ew, "ewas_candidates.csv"), row.names = FALSE)
cat(sprintf("数值列 %d;规则剔除 %d;进入清洗 %d;清洗后不满足纳入 %d;分析 %d(连续 %d / 二分类 %d / 多分类 %d)\n", nrow(V), sum(V$excl != ""), nrow(cand),
            length(FILT), nrow(R), sum(R$type == "continuous"), sum(R$type == "categorical" & R$levels == 2), sum(R$type == "categorical" & R$levels > 2)))
print(table(V$excl[V$excl != ""]))
if (STAGE == "list") { write.csv(R[, c("file", "variable", "label", "section", "subsection", "type", "levels", "n")], file.path(ew, "ewas_analysable_preview.csv"), row.names = FALSE); quit(save = "no") }

# ---------- 5. 领域 ----------
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
    if (v == "PHPINSU") return("Diabetes history")                          # 采血前问卷"是否在用胰岛素"
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
# 大类:身体功能单列(作为结局分析),骨与关节(含骨代谢生化)为另一大类
GROUP <- c("Musculoskeletal (joints, bone, fracture)" = "Bone & joints", "Physical function & mobility" = "Physical function",
           "Bone & mineral biochemistry" = "Bone & joints", "Anthropometry & body composition" = "Body size & composition",
           "Physical activity" = "Behaviour", "Tobacco & alcohol" = "Behaviour", "Diet" = "Behaviour",
           "Cardiometabolic biomarkers" = "Cardiometabolic & inflammatory", "Inflammation" = "Cardiometabolic & inflammatory",
           "Cardiovascular & blood pressure" = "Cardiometabolic & inflammatory", "Diabetes history" = "Cardiometabolic & inflammatory")
R$group <- ifelse(R$domain %in% names(GROUP), GROUP[R$domain], "Other")

# ---------- 6. 两个扫描范围 ----------
# 主扫描不含作为结局分析的骨折、跌倒与身体功能;全部变量的扫描作敏感性
FXF <- c("HAG5A", "HAG5B", "HAG5C", grep("^HAG(6|10)", R$variable, value = TRUE), "HAG14", "HAG15")
R$outcome_item <- R$domain == "Physical function & mobility" | R$variable %in% FXF
M <- do.call(cbind, X[R$variable])
t0 <- Sys.time(); CM <- suppressWarnings(cor(M, use = "pairwise.complete.obs")); CM[is.na(CM)] <- 0; diag(CM) <- 1
cat(sprintf("\n相关矩阵 %d×%d(%.0f 秒)\n", ncol(M), ncol(M), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat(sprintf("清洗:codebook 未列出而设缺失 %d 个值(%d 个变量);缺失码设缺失 %d;7-fill 重编码 %d\n",
            sum(R$n_uncovered), sum(R$n_uncovered > 0), sum(R$n_special), sum(R$n_recoded)))

show <- function(x, k) {
  x <- head(x, k)
  for (i in seq_len(nrow(x))) with(x[i, ], cat(sprintf("%3d %-9s %-44s %-34s n %4d  %s  pr %.4f  P %s  FDR %s%s\n", rank, variable, substr(label, 1, 44), substr(domain, 1, 34), n,
    if (is.na(beta_y)) sprintf("%-40s", substr(contrast, 1, 40)) else sprintf("%+.2f y (%+.2f, %+.2f) β %+.3f %-12s", beta_y, lo_y, hi_y, std_beta, substr(contrast, 1, 12)),
    partial_r2, fp(p), fp(fdr), if (!is.na(lead) && lead != variable) paste0("  ~", lead) else "")))
}
# 多重比较(BH、Bonferroni、Li–Ji 有效检验数)、|r| > 0.7 贪心归并独立信号、领域与大类汇总;均在该范围内重算
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
  cat(sprintf("检验 %d(连续 %d / 二分类 %d / 多分类 %d);Li–Ji 有效检验数 %.0f;FDR 线 P = %.2e;Meff 线 P = %.1e\n", nrow(R), sum(R$type == "continuous"),
              sum(R$type == "categorical" & R$levels == 2), sum(R$type == "categorical" & R$levels > 2), meff, max(c(0, R$p[R$fdr < 0.05])), 0.05 / meff))
  cat(sprintf("FDR < 0.05:%d 个(独立信号 %d);Meff 校正 < 0.05:%d;Bonferroni < 0.05:%d\n", sum(R$fdr < 0.05), sum(R$is_lead), sum(R$p_meff < 0.05), sum(R$p_bonf < 0.05)))
  cat("\n--- 前 30(按 P;~ 表示与更前面的独立信号 |r| > 0.7)---\n"); show(R, 30)
  cat("\n--- 独立信号前 30 ---\n"); show(R[R$is_lead, ], 30)
  D <- do.call(rbind, lapply(split(R, R$domain), function(x) data.frame(domain = x$domain[1], group = x$group[1], n_tested = nrow(x), n_fdr = sum(x$fdr < 0.05),
    pct_fdr = 100 * mean(x$fdr < 0.05), n_lead = sum(x$is_lead), best_rank = min(x$rank), best_var = x$variable[which.min(x$rank)],
    max_partial_r2 = max(x$partial_r2), median_partial_r2 = median(x$partial_r2))))
  tot_sig <- sum(R$fdr < 0.05); D$enrich_or <- NA; D$enrich_p <- NA
  for (i in seq_len(nrow(D))) {
    a <- D$n_fdr[i]; b <- D$n_tested[i] - a; c <- tot_sig - a; dd <- (nrow(R) - D$n_tested[i]) - c
    ft <- fisher.test(matrix(c(a, b, c, dd), 2)); D$enrich_or[i] <- unname(ft$estimate); D$enrich_p[i] <- ft$p.value
  }
  D <- D[order(D$best_rank), ]
  cat("\n--- 各领域(按最靠前排名;富集为 FDR 显著比例的 Fisher 检验)---\n")
  for (i in seq_len(nrow(D))) with(D[i, ], cat(sprintf("  %-40s 检验 %3d  FDR 显著 %3d (%4.1f%%)  独立信号 %2d  最佳 #%-3d %-9s 最大偏 R² %.4f  富集 OR %.2f P %s\n",
    domain, n_tested, n_fdr, pct_fdr, n_lead, best_rank, best_var, max_partial_r2, enrich_or, fp(enrich_p))))
  LD <- R[R$is_lead, ]; G <- sort(unique(R$group))
  cat("\n--- 大类构成(独立信号按 P 排序的前 k 个;末列为全部检验)与排名检验(该类 P 排名 vs 其余,Wilcoxon)---\n")
  for (g in G) { ing <- R$group == g; wp <- if (sum(ing) > 1) wilcox.test(R$rank[ing], R$rank[!ing])$p.value else NA
    cat(sprintf("  %-32s 前10 %2d 前20 %2d 前50 %2d 前100 %3d | 全部 %3d (%2.0f%%) | 中位排名 %4.0f vs %4.0f P %s\n", g, sum(head(LD$group, 10) == g),
                sum(head(LD$group, 20) == g), sum(head(LD$group, 50) == g), sum(head(LD$group, 100) == g), sum(ing), 100 * mean(ing),
                median(R$rank[ing]), median(R$rank[!ing]), if (is.na(wp)) "—" else fp(wp))) }
  write.csv(R, file.path(ew, paste0("ewas_results", suffix, ".csv")), row.names = FALSE)
  write.csv(D, file.path(ew, paste0("ewas_domains", suffix, ".csv")), row.names = FALSE)
  write.csv(R[R$is_lead, ], file.path(ew, paste0("ewas_top", suffix, ".csv")), row.names = FALSE)
  invisible(R)
}
RM <- finish(R[!R$outcome_item, ], sprintf("主扫描(图 4b):除骨折、跌倒与身体功能外的 %d 个变量", sum(!R$outcome_item)), "_main")
RA <- finish(R, sprintf("敏感性:全部 %d 个变量", nrow(R)), "")
saveRDS(CM, file.path(ew, "ewas_cor.rds"))
# 通过有效检验数校正的独立信号(全扫描与主扫描领先信号的并集):清洗后原值 + 类型信息,供反向建模(变量 ~ KSAA)
TL <- union(RA$variable[RA$is_lead & RA$p_meff < 0.05], RM$variable[RM$is_lead & RM$p_meff < 0.05])
saveRDS(list(SEQN = d$SEQN, X = XR[TL], meta = RA[match(TL, RA$variable), ]), file.path(ew, "ewas_clean_top.rds"))

# ---------- 7. 稳健性(全部变量)----------
cat("\n=========== 稳健性 ===========\n")
cc <- RA[RA$type == "continuous", ]
cat(sprintf("  连续变量 %d 个:秩逆正态变换后方向一致 %d / %d;原 FDR 显著者在 INT 下 P < 0.05:%d / %d\n", nrow(cc),
            sum(sign(cc$beta_y) == sign(cc$beta_int)), nrow(cc), sum(cc$fdr < 0.05 & cc$p_int < 0.05), sum(cc$fdr < 0.05)))
top <- rbind(head(RA, 20), head(RM[!RM$variable %in% head(RA$variable, 20), ], 10))
cat("  取值检查(全部扫描前 20 + 主扫描另外的前 10;连续:清洗后分位数;分类:各级人数)—— 有无残留缺失码:\n")
for (i in seq_len(nrow(top))) {
  v <- top$variable[i]; x <- clean_one(RAW[[v]], C[C$file == top$file[i] & C$variable == v, ], top$label[i])$x
  if (top$type[i] == "continuous") cat(sprintf("    %-9s min %g | p1 %g | 中位 %g | p99 %g | max %g | 缺失码剔除 %d 未列出值 %d\n", v, min(x, na.rm = TRUE),
    quantile(x, .01, na.rm = TRUE), median(x, na.rm = TRUE), quantile(x, .99, na.rm = TRUE), max(x, na.rm = TRUE), top$n_special[i], top$n_uncovered[i]))
  else { tb <- table(x); cat(sprintf("    %-9s %s | 缺失码剔除 %d\n", v, paste(sprintf("%s=%d", names(tb), tb), collapse = " "), top$n_special[i])) }
}
cat("  与 S16 同一测量的对照(本扫描 β / S16 β,岁每 SD):\n")
s16 <- read.csv(file.path(out, "correlates.csv"))
for (pr in list(c("BDPFNBMD", "Femoral neck BMD"), c("BDPTOBMD", "Total hip BMD"), c("BMPBMI", "BMI"), c("GHP", "HbA1c"), c("BMPWAIST", "Waist circumference"), c("WCP", "White blood cell count")))
  if (pr[1] %in% RA$variable) cat(sprintf("    %-9s %+.2f(n %d)/ %s %+.2f(n %d)\n", pr[1], RA$beta_y[RA$variable == pr[1]], RA$n[RA$variable == pr[1]], pr[2],
    s16$beta_y[s16$factor == pr[2]][1], s16$n[s16$factor == pr[2]][1]))
cat("\n-> ewas_candidates.csv, ewas_{results,domains,top}{_main,}.csv, ewas_cor.rds, ewas_clean_top.rds\n")
