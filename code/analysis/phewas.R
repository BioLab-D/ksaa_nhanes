# PheWAS:按审核表(phewas_audit.R)的加权全表型扫描,方向 KSAA → 表型。
# 模型:每个表型一个,暴露 = KSAA 每 SD(全队列 SD),协变量固定为年龄、性别、种族/族裔、教育(edu3)、PIR(不按单变量显著性增减)。
#   连续 → svyglm 高斯(偏态先取对数再标准化);二分类 → quasibinomial(OR);有序 ≥3 级 → svyolr(稀疏级 <30 并入相邻级);名义 → 反向联合检验,只报 P。
#   设计:WTPFEX6 + SDPPSU6/SDPSTRA6(nest);空腹亚样本变量用 WTPFSD6;生殖史题限女性、不调整性别;子集一律 subset(design)。
# 复制:phase 1 / phase 2 各用相位设计(WTPFEXn + SDPPSUn/SDPSTRAn)拟合,"复制" = 两相位同向且各 P < 0.05。
#   P 与 CI 用 t 分布、自由度 = 设计自由度(6 年设计 49,相位设计 23);不用 6 年设计切相位(层内只剩 1 个 PSU)。
# 多重校正:主分析集、补充集、扩展模块各自做 BH-FDR、Bonferroni、Li–Ji 有效检验数。
# 输入 results/analysis/{analysis_data.rds, fxf_vars.csv, ewas/ewas_clean.rds, ewas/ewas_audit_table.csv, ewas/ewas_multilevel_plan.csv}
# 输出 results/analysis/phewas/phewas_results.csv、phewas_y.rds、phewas_summary.txt
suppressMessages({library(survey); library(MASS); library(DBI); library(duckdb)})
options(survey.lonely.psu = "adjust")
res <- here::here("results"); out <- file.path(res, "analysis"); ew <- file.path(out, "phewas")
d <- readRDS(file.path(out, "analysis_data.rds"))
CL <- readRDS(file.path(ew, "ewas_clean.rds")); stopifnot(identical(CL$SEQN, d$SEQN))
A <- read.csv(file.path(ew, "ewas_audit_table.csv"), stringsAsFactors = FALSE)
A <- A[A$in_pool & A$status %in% c("main", "supplementary", "exposure_module"), ]
PL <- read.csv(file.path(ew, "ewas_multilevel_plan.csv"), stringsAsFactors = FALSE); CLS <- setNames(PL$class, PL$variable)
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
dbWriteTable(con, "coh", data.frame(SEQN = d$SEQN), temporary = TRUE)
W <- dbGetQuery(con, "select SEQN, WTPFSD6, WTPFEX1, WTPFSD1, SDPPSU1, SDPSTRA1, WTPFEX2, WTPFSD2, SDPPSU2, SDPSTRA2 from nhanes3_adult join coh using (SEQN)"); dbDisconnect(con, shutdown = TRUE)
W <- W[match(d$SEQN, W$SEQN), ]; for (nm in setdiff(names(W), "SEQN")) { d[[nm]] <- W[[nm]]; d[[nm]][is.na(d[[nm]])] <- 0 }
stopifnot(all(d$WTPFEX1[d$phase == 1] > 0), all(d$WTPFEX2[d$phase == 2] > 0))
SDG <- sd(d$gap); d$z <- d$gap / SDG
d$sex <- factor(d$sex); d$race <- factor(d$race); d$edu3 <- factor(d$edu3)
COV <- c("age", "sex", "race", "edu3", "pir")
mkdes <- function(w, psu, str, rows) svydesign(id = as.formula(paste0("~", psu)), strata = as.formula(paste0("~", str)), weights = as.formula(paste0("~", w)), nest = TRUE, data = d[rows, ])
# 设计:6 年设计、相位设计、空腹亚样本权重;rows 用于对齐 y
DS <- list(all = list(des = mkdes("WTPFEX6", "SDPPSU6", "SDPSTRA6", TRUE), rows = rep(TRUE, nrow(d))),
           p1 = list(des = mkdes("WTPFEX1", "SDPPSU1", "SDPSTRA1", d$phase == 1), rows = d$phase == 1),
           p2 = list(des = mkdes("WTPFEX2", "SDPPSU2", "SDPSTRA2", d$phase == 2), rows = d$phase == 2),
           fall = list(des = mkdes("WTPFSD6", "SDPPSU6", "SDPSTRA6", d$WTPFSD6 > 0), rows = d$WTPFSD6 > 0),
           fp1 = list(des = mkdes("WTPFSD1", "SDPPSU1", "SDPSTRA1", d$WTPFSD1 > 0), rows = d$WTPFSD1 > 0),
           fp2 = list(des = mkdes("WTPFSD2", "SDPPSU2", "SDPSTRA2", d$WTPFSD2 > 0), rows = d$WTPFSD2 > 0))
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3

# ---------- 合成得分:故事回忆即时 / 延迟各 6 段 → 0–6;地址定向 4 题 → 0–4 ----------
lab_of <- function(v) CL$labs[[v]]
score <- function(vs, re) { M <- sapply(vs, function(v) { x <- CL$x[[v]]; lb <- lab_of(v); yes <- as.numeric(names(lb)[grepl(re, lb, ignore.case = TRUE)])
  ifelse(is.na(x), NA, as.numeric(x %in% yes)) }); ifelse(rowSums(is.na(M)) > 0, NA, rowSums(M)) }
COMP <- list(MAPB_story_recall_immediate_total = list(items = paste0("MAPB2", LETTERS[1:6]), re = "^yes", label = "Story recall, immediate (details recalled, 0-6)"),
             MAPD_story_recall_delayed_total = list(items = paste0("MAPD2", LETTERS[1:6]), re = "^yes", label = "Story recall, delayed (details recalled, 0-6)"),
             HAA_orientation_score = list(items = paste0("HAA6", LETTERS[1:4]), re = "^correct", label = "Orientation: address items correct (0-4)"))
for (nm in names(COMP)) { it <- COMP[[nm]]; stopifnot(all(it$items %in% A$variable))
  CL$x[[nm]] <- score(it$items, it$re); CL$cont[[nm]] <- TRUE; CL$labs[[nm]] <- character()
  r <- A[A$variable == it$items[1], ]; r$variable <- nm; r$label <- it$label; r$type <- "continuous"; r$summarise_to <- ""; r$reason <- paste("composite of", length(it$items), "items")
  A <- rbind(A[!(A$variable %in% it$items), ], r) }

# ---------- 其余重复测量 / 多题测验的汇总,定义与表 S10 一致 ----------
# 新构造表型同样过门槛:连续 ≥500 且 ≥5 个不同值;二分类 ≥500 且每类 ≥100;有序 ≥500。
gate <- function(y, type) { y <- y[!is.na(y)]; if (type == "continuous") length(y) >= 500 && length(unique(y)) >= 5
  else if (length(unique(y)) == 2) length(y) >= 500 && min(table(y)) >= 100 else length(y) >= 500 }
add_pheno <- function(nm, y, type, label, tmpl, cls, src, note, labs = character()) {
  stopifnot(length(y) == nrow(d), !(nm %in% A$variable))
  r <- A[A$variable == tmpl, ]; if (nrow(r) == 0) { r <- read.csv(file.path(ew, "ewas_audit_table.csv"), stringsAsFactors = FALSE); r <- r[r$variable == tmpl, ] }
  stopifnot(nrow(r) == 1)
  r$variable <- nm; r$label <- label; r$type <- type; r$health_class <- cls; r$summarise_to <- ""; r$s8_version <- ""; r$reason <- note
  r$status <- if (gate(y, type)) "main" else "supplementary"
  CL$x[[nm]] <<- y; CL$cont[[nm]] <<- type == "continuous"; CL$labs[[nm]] <<- labs
  A <<- rbind(A[!(A$variable %in% src), ], r) }
xv <- function(v) CL$x[[v]]
lev_is <- function(v, re) { ll <- CL$labs[[v]][as.character(CL$x[[v]])]; ifelse(is.na(CL$x[[v]]), NA, as.numeric(grepl(re, ll, ignore.case = TRUE))) }
# 重复测量:步速取两次较快(= 表 S10 gait,逐位核对)、步数取均值、脉搏取家访与体检均值
sp <- 2.4384 / suppressWarnings(pmin(xv("PFPTWLKA"), xv("PFPTWLKB"), na.rm = TRUE)); sp[!is.finite(sp)] <- NA
stopifnot(identical(is.na(sp), is.na(d$gait)), max(abs(sp - d$gait), na.rm = TRUE) < 1e-9)
add_pheno("PFP_gait_speed", sp, "continuous", "Gait speed (faster of two 8-ft walks, m/s)", "PFPTWLKA", "function", c("PFPTWLKA", "PFPTWLKB"),
          "v1.2: repeated trials summarised (faster of two), identical to the Table S8 gait speed")
mn2 <- function(a, b) { m <- rowMeans(cbind(a, b), na.rm = TRUE); m[is.nan(m)] <- NA; m }
add_pheno("PFP_steps_8ft_mean", mn2(xv("PFPNSTPA"), xv("PFPNSTPB")), "continuous", "Steps to walk 8 ft (mean of two trials)", "PFPNSTPA", "function",
          c("PFPNSTPA", "PFPNSTPB"), "v1.2: repeated trials summarised (mean)")
add_pheno("pulse_mean", mn2(xv("HAZA5R"), xv("PEP6DR")), "continuous", "Pulse rate (mean of household and examination readings)", "HAZA5R", "objective",
          c("HAZA5R", "PEP6DR"), "v1.2: repeated measurement summarised (mean of household and MEC readings, as for blood pressure)")
# 连减 5 题 → 答对题数(正确答案 17/14/11/8/5):55 "can't do"、99 "don't know" 记答错,88 与空白为缺失,任一题缺失则得分缺失。
#   清洗表把 55/99 置为缺失,故从 duckdb 取原始码。
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE); dbWriteTable(con, "coh", data.frame(SEQN = d$SEQN), temporary = TRUE)
H18 <- dbGetQuery(con, "select SEQN, HAP18A, HAP18B, HAP18C, HAP18D, HAP18E from nhanes3_adult join coh using (SEQN)"); dbDisconnect(con, shutdown = TRUE)
H18 <- H18[match(d$SEQN, H18$SEQN), ]; KEY <- c(HAP18A = 17, HAP18B = 14, HAP18C = 11, HAP18D = 8, HAP18E = 5)
M18 <- sapply(names(KEY), function(v) { r <- H18[[v]]; ifelse(is.na(r) | r == 88, NA, as.numeric(r == KEY[[v]])) })
for (v in names(KEY)) stopifnot(all(xv(v)[!is.na(xv(v))] == H18[[v]][!is.na(xv(v))]))   # 清洗表的有效值与原始码一致
add_pheno("HAP18_serial_subtraction", ifelse(rowSums(is.na(M18)) > 0, NA, rowSums(M18)), "continuous", "Serial subtraction of $3, answers correct (0-5)", "HAP18A", "objective",
          names(KEY), "v1.2: 5 items scored correct / incorrect ('can't do' and 'don't know' = incorrect); replaces the answers analysed as continuous values")
sc3 <- function(vs) { M <- sapply(vs, lev_is, re = "^correct"); ifelse(rowSums(is.na(M)) > 0, NA, rowSums(M)) }
L03 <- c(`0` = "0 correct", `1` = "1 correct", `2` = "2 correct", `3` = "3 correct")
add_pheno("HAP17_word_registration", sc3(c("HAP17A1", "HAP17A2", "HAP17A3")), "categorical", "Three-word registration, first trial (words correct, 0-3)", "HAP17A2", "objective",
          c("HAP17A1", "HAP17A2", "HAP17A3"), "v1.2: 3 items summed (as for story recall); ordinal (fewer than 5 values)", L03)
add_pheno("HAP19_word_recall", sc3(c("HAP19A", "HAP19B", "HAP19C")), "categorical", "Three-word recall (words correct, 0-3)", "HAP19A", "objective",
          c("HAP19A", "HAP19B", "HAP19C"), "v1.2: 3 items summed (as for story recall); ordinal (fewer than 5 values)", L03)
# 与表 S10 定义统一,原变量就地重编码(每个变量在扫描中只出现一次):
#   困难程度题(HAH1–HAH12、PEP13E1A–E5A)→ 任一困难(含做不到)vs 无;跌倒 HAG14 → 过去 12 个月 ≥1 次;
#   串联站立 PFPTTAND → 不能坚持 10 秒(含未能尝试者);骨折题不变。新编码与表 S10 派生变量逐位核对。
FXV <- read.csv(file.path(out, "fxf_vars.csv")); FXV <- FXV[match(d$SEQN, FXV$SEQN), ]
cmp <- function(a, b) { k <- !is.na(a) & !is.na(b); stopifnot(all(a[k] == b[k])); invisible(TRUE) }
recode_bin <- function(v, y, lab0, lab1, note) { stopifnot(v %in% A$variable, length(y) == nrow(d))
  CL$x[[v]] <<- as.numeric(y); CL$cont[[v]] <<- FALSE; CL$labs[[v]] <<- setNames(c(lab0, lab1), c("0", "1"))
  i <- which(A$variable == v); A$type[i] <<- "categorical"; A$reason[i] <<- note
  A$status[i] <<- if (gate(as.numeric(y), "categorical")) "main" else "supplementary" }
DIFF <- c(paste0("HAH", 1:12), paste0("PEP13E", 1:5, "A"))
for (v in DIFF) { ll <- CL$labs[[v]][as.character(CL$x[[v]])]; stopifnot(any(grepl("^no difficulty", ll, ignore.case = TRUE)))
  y <- 1 - lev_is(v, "^no difficulty")
  if (v %in% c("HAH1", "HAH2", "HAH8", "HAH9")) { cmp(FXV[[paste0("dif_", v)]], y); stopifnot(identical(is.na(FXV[[paste0("dif_", v)]]), is.na(y))) }
  recode_bin(v, y, "No difficulty", "Any difficulty", "v1.2: any difficulty (including unable) vs none, as in Table S8") }
cmp(d$fall, as.numeric(xv("HAG14") >= 1)); stopifnot(identical(is.na(d$fall), is.na(xv("HAG14"))))
recode_bin("HAG14", d$fall, "No", "Yes", "v1.2: at least one fall in the past 12 months, as in Table S8")
cmp(FXV$tand_fail, as.numeric(xv("PFPTTAND") < 10))
recode_bin("PFPTTAND", FXV$tand_fail, "No", "Yes", "v1.2: unable to hold the tandem stand for 10 s, as in Table S8")
cat(sprintf("审核表 v1.2:主分析表型 %d(补充 %d;扩展模块 %d)\n", sum(A$status == "main"), sum(A$status == "supplementary"), sum(A$status == "exposure_module")))

# ---------- 逐变量类型判定与重编码(规则事先写定,理由随结果输出)----------
NA3 <- c(HAF2 = "Never walks uphill or hurries", HAL14D = "Doesn't work", HAL17 = "Never received allergy shots")   # 第三级 = 不适用 → 缺失,二分类
ORDMAP <- list(HAK4 = c(`3` = 0, `1` = 1, `2` = 2),          # never < 1-2 < 3+ times
               PFPSCOOT = c(`2` = 0, `1` = 1, `3` = 2))      # no < yes(uses arms) < unable
ORD_AS_IS <- c("PFPUSTND", "PFPLSER", "PFPRSER", "ANP", "HZP", "MIP", "DEPDQL1", "DEPDQU1")   # 计划表误标 NOM/LAT/YNA,实为有序
BIN_RULE <- list(DEPLPTYP = list(one = "partial", zero = "full", drop = "unable", note = "partial vs full denture; 'unable to evaluate' -> missing"),
                 DEPUPTYP = list(one = "partial", zero = "full", drop = "unable", note = "partial vs full denture; 'unable to evaluate' -> missing"),
                 GUPFDX1R = list(one = "gallstone|cholecystectomy", zero = "^normal", drop = "no conclusion|abnormal", note = "gallbladder disease (stones or cholecystectomy) vs normal; inconclusive/other -> missing"),
                 HAP2 = list(one = "eyeglasses|contact|both", zero = "neither", drop = "^$", note = "any vision correction vs none"))
merge_sparse <- function(y, min_n = 30) {   # 有序:<30 的级并入相邻级(端点并向内,中间并向人数少的一侧)
  notes <- character(); repeat { tb <- table(y); if (length(tb) <= 2 || min(tb) >= min_n) break
    i <- which.min(tb); L <- names(tb); j <- if (i == 1) 2 else if (i == length(tb)) i - 1 else if (tb[i - 1] <= tb[i + 1]) i - 1 else i + 1
    a <- L[min(i, j)]; b <- L[max(i, j)]; new <- paste(a, b, sep = "+"); notes <- c(notes, sprintf("merged [%s]+[%s]", a, b))
    y <- as.character(y); y[y %in% c(a, b)] <- new; lv <- L; lv[lv %in% c(a, b)] <- new; y <- factor(y, unique(lv), ordered = TRUE) }
  list(y = y, notes = notes) }

prep <- function(v, r) {
  x <- CL$x[[v]]; lb <- CL$labs[[v]]; notes <- character(); typ <- r$type; y <- x
  if (typ == "continuous") { if (min(y, na.rm = TRUE) >= 0 && skew(y[!is.na(y)]) > 1) { c0 <- if (min(y, na.rm = TRUE) > 0) 0 else min(y[y > 0], na.rm = TRUE) / 2; y <- log(y + c0); notes <- c(notes, "log") }
    y <- (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE); return(list(y = y, model = "gaussian", notes = notes, levels = "")) }
  ll <- lb[as.character(x)]; ll[is.na(x)] <- NA; k <- length(unique(x[!is.na(x)])); cls <- if (v %in% names(CLS)) CLS[[v]] else "ORD"
  if (v %in% names(BIN_RULE)) { b <- BIN_RULE[[v]]; y <- ifelse(is.na(ll), NA, ifelse(grepl(b$one, ll, ignore.case = TRUE), 1, ifelse(grepl(b$zero, ll, ignore.case = TRUE), 0, NA))); notes <- b$note; k <- 2 }
  else if (v %in% names(NA3)) { y <- ifelse(is.na(ll), NA, ifelse(grepl("^yes", ll, ignore.case = TRUE), 1, ifelse(grepl("^no$", ll, ignore.case = TRUE), 0, NA))); notes <- sprintf("'%s' -> missing", NA3[[v]]); k <- 2 }
  else if (v %in% names(ORDMAP)) { m <- ORDMAP[[v]]; y <- unname(m[as.character(x)]); notes <- "levels re-ordered to a monotone scale"; k <- length(m) }
  else if (cls == "LAT" && any(grepl("right|left", ll, ignore.case = TRUE))) { y <- ifelse(is.na(ll), NA, ifelse(grepl("^no\\b", ll, ignore.case = TRUE), 0, ifelse(grepl("both", ll, ignore.case = TRUE), 2, 1))); notes <- "none < unilateral < bilateral"; k <- 3 }
  else if (cls == "NOM" && !(v %in% ORD_AS_IS) && k > 2) { y <- factor(ll); return(list(y = y, model = "nominal", notes = "nominal: reverse joint test only", levels = paste(sprintf("%s=%d", names(table(y)), as.integer(table(y))), collapse = " | "))) }
  if (k == 2) { u <- sort(unique(y[!is.na(y)])); if (!identical(as.numeric(u), c(0, 1))) {   # 参照级:No/None/Never/Normal 优先,否则最小编码
      l2 <- lb[as.character(u)]; neg <- which(grepl("^(no|none|never|normal|not)\\b|no difficulty|no findings|no pain|^correct|^yes, recalled", l2, ignore.case = TRUE))
      ref <- if (length(neg)) u[neg[1]] else u[1]; y <- as.numeric(y != ref); notes <- c(notes, sprintf("1 = %s", lb[as.character(setdiff(u, ref))])) }
    tb <- table(y); return(list(y = y, model = "binomial", notes = notes, levels = sprintf("0=%d | 1=%d", tb[["0"]], tb[["1"]]))) }
  yo <- factor(y, sort(unique(y[!is.na(y)])), ordered = TRUE); if (v %in% names(lb) || TRUE) levels(yo) <- ifelse(is.na(lb[levels(yo)]) | !nzchar(lb[levels(yo)]), levels(yo), lb[levels(yo)])
  ms <- merge_sparse(yo); yo <- ms$y; notes <- c(notes, ms$notes)
  if (nlevels(yo) == 2) { y <- as.numeric(yo) - 1; tb <- table(y); return(list(y = y, model = "binomial", notes = c(notes, "2 levels after merge"), levels = sprintf("0=%d | 1=%d", tb[["0"]], tb[["1"]]))) }
  list(y = yo, model = "ordinal", notes = notes, levels = paste(sprintf("%s=%d", levels(yo), as.integer(table(yo))), collapse = " | "))
}

fit <- function(model, des, cov) {
  v <- des$variables; cov <- cov[sapply(cov, function(cv) length(unique(v[[cv]])) > 1)]      # 子集内单一水平的协变量(单性别题的 sex)剔除
  df <- degf(des); tp <- function(b, se) 2 * pt(-abs(b / se), df)
  f <- as.formula(paste("y ~ z +", paste(cov, collapse = " + ")))
  err <- ""; note <- ""; used <- model    # used = 实际拟合的模型(有序拟合失败会回退成二分类)
  o <- tryCatch({
    if (model == "gaussian") { m <- svyglm(f, des); s <- summary(m)$coefficients; c(s["z", 1], s["z", 2], tp(s["z", 1], s["z", 2])) }
    else if (model == "binomial") { m <- svyglm(f, des, family = quasibinomial()); s <- summary(m)$coefficients; c(s["z", 1], s["z", 2], tp(s["z", 1], s["z", 2])) }
    else if (model == "ordinal") {
      r <- tryCatch({ m <- svyolr(f, des); b <- coef(m)[["z"]]; se <- sqrt(vcov(m)["z", "z"]); c(b, se, tp(b, se)) }, error = function(e) NULL)
      if (is.null(r)) { des$variables$y <- as.numeric(as.integer(des$variables$y) > 1)
        note <- "svyolr singular -> refitted as binomial: any level above the lowest vs the lowest"; used <- "binomial (ordinal fallback)"
        m <- svyglm(f, des, family = quasibinomial()); s <- summary(m)$coefficients; r <- c(s["z", 1], s["z", 2], tp(s["z", 1], s["z", 2])) }
      r }
    else { m <- svyglm(as.formula(paste("z ~ y +", paste(cov, collapse = " + "))), des); c(NA, NA, regTermTest(m, ~y)$p[1]) }
  }, error = function(e) { err <<- conditionMessage(e); c(NA, NA, NA) })
  list(b = o[1], se = o[2], p = o[3], err = err, note = note, df = df, used = used)
}

R <- list(); YN <- list()
run <- function(key, y, model, keep, cov) { D <- DS[[key]]; n <- sum(keep & D$rows)
  if (n < 50) return(list(fit = list(b = NA, se = NA, p = NA, err = if (n == 0) "not asked in this phase" else sprintf("n = %d < 50", n), note = "", df = NA, used = model), n = n))   # 相位特异题
  des <- D$des; des$variables$y <- y[D$rows]; ds <- subset(des, keep[D$rows]); list(fit = fit(model, ds, cov), n = n) }
for (i in seq_len(nrow(A))) {
  v <- A$variable[i]; r <- A[i, ]; pr <- prep(v, r)
  women <- r$health_class == "reproductive_exposure"
  keep <- !is.na(pr$y) & complete.cases(d[, COV]) & (!women | d$sex == "2")
# 空腹亚样本(两条并集):① 糖耐量组的血糖 / 胰岛素 / C 肽(HbA1c 除外),非晨间者不是空腹值,限定晨间亚样本并用其权重;
#   ② 只在晨间亚样本里有值的变量(如 LDL,LCP)。
  fast <- (grepl("DIABETES TESTING", r$section) && v != "GHP") || (sum(keep) >= 50 && mean(d$WTPFSD6[keep] > 0) > 0.95)
  cov <- if (women) setdiff(COV, "sex") else COV
  a1 <- run(if (fast) "fp1" else "p1", pr$y, pr$model, keep, cov); a2 <- run(if (fast) "fp2" else "p2", pr$y, pr$model, keep, cov)
# 只在一期采集的项目:主分析用该期设计(该期权重 + PSU/层);判定只看两期有效样本量。
  key <- if (a1$n >= 50 && a2$n < 50) (if (fast) "fp1" else "p1") else if (a2$n >= 50 && a1$n < 50) (if (fast) "fp2" else "p2") else (if (fast) "fall" else "all")
  single <- key != (if (fast) "fall" else "all")
  a <- run(key, pr$y, pr$model, keep, cov)
  ex <- pr$model %in% c("binomial", "ordinal"); tr <- function(b) if (ex) exp(b) else b; f <- a$fit; f1 <- a1$fit; f2 <- a2$fit
  R[[v]] <- data.frame(variable = v, phenotype = r$label, set = r$status, health_class = r$health_class, model = pr$model, n = a$n, n_phase1 = a1$n, n_phase2 = a2$n,
    levels = pr$levels,
    weight = if (single) sprintf("%s (single-phase item)", c(p1 = "WTPFEX1", p2 = "WTPFEX2", fp1 = "WTPFSD1", fp2 = "WTPFSD2")[[key]])
             else if (fast) "WTPFSD6 (phase: WTPFSD1/2)" else "WTPFEX6 (phase: WTPFEX1/2)",
    population = if (women) "women" else "all", df = f$df,
    measure = c(gaussian = "beta (SD of phenotype per SD KSAA)", binomial = "OR per SD KSAA", ordinal = "OR per SD KSAA (proportional odds)",
                `binomial (ordinal fallback)` = "OR per SD KSAA (any level above the lowest vs the lowest)", nominal = "joint test only")[[f$used]],
    model_used = f$used,
    est = tr(f$b), lo = tr(f$b - qt(0.975, f$df) * f$se), hi = tr(f$b + qt(0.975, f$df) * f$se), p = f$p,   # CI 与 P 用同一自由度的 t 分布
    est_phase1 = tr(f1$b), p_phase1 = f1$p, est_phase2 = tr(f2$b), p_phase2 = f2$p,
    replicated = if (single) NA else (!is.na(f1$b) && !is.na(f2$b) && sign(f1$b) == sign(f2$b) && f1$p < 0.05 && f2$p < 0.05),   # 单期项目无法做两期复制
    notes = paste(c(pr$notes, f$note, if (single) "single-phase item: main analysis uses that phase's design; no cross-phase replication" else NULL), collapse = "; "), error = f$err, section = r$section, subsection = r$subsection, stringsAsFactors = FALSE)
  YN[[v]] <- if (pr$model == "nominal") NULL else { z <- rep(NA_real_, nrow(d)); z[keep] <- as.numeric(pr$y)[keep]; z }
  if (i %% 50 == 0) cat(sprintf("  %d / %d\n", i, nrow(A)))
}
R <- do.call(rbind, R); rownames(R) <- NULL
FASTV <- sort(R$variable[R$set == "main" & grepl("WTPFSD", R$weight)])
cat("空腹亚样本变量(主集):", paste(FASTV, collapse = ", "), "\n")
stopifnot(setequal(FASTV, c("C1P", "C2P", "G1P", "G2P", "I1P", "I2P", "LCP")))   # 六项糖耐量/胰岛素 + LDL;HbA1c(GHP)不在其中

# ---------- 多重校正:各集内 ----------
liji <- function(M) { M <- M[, colSums(!is.na(M)) > 0, drop = FALSE]; C <- suppressWarnings(cor(M, use = "pairwise.complete.obs")); C[is.na(C)] <- 0; diag(C) <- 1
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values; ev[ev < 0] <- 0; sum(as.numeric(ev >= 1) + (ev - floor(ev))) }
for (s in unique(R$set)) { k <- R$set == s & !is.na(R$p); R$fdr[k] <- p.adjust(R$p[k], "BH"); R$bonf[k] <- p.adjust(R$p[k], "bonferroni")
  vs <- R$variable[k & R$variable %in% names(YN)]; M <- do.call(cbind, YN[vs]); me <- liji(M); R$meff[k] <- me; R$p_liji[k] <- pmin(1, R$p[k] * me) }
# 独立信号:主分析 FDR<5% 者按 P 升序贪心归并(与已选信号 |r| > 0.7 者归入其簇;r 用清洗后数值版,成对完整)
k <- which(R$set == "main" & !is.na(R$fdr) & R$fdr < 0.05 & R$variable %in% names(YN)); k <- k[order(R$p[k])]
R$clump_lead <- NA_character_; leads <- character()
for (i in k) { v <- R$variable[i]; hit <- NA
  for (l in leads) { r <- suppressWarnings(cor(YN[[v]], YN[[l]], use = "pairwise.complete.obs")); if (!is.na(r) && abs(r) > 0.7) { hit <- l; break } }
  if (is.na(hit)) { leads <- c(leads, v); R$clump_lead[i] <- v } else R$clump_lead[i] <- hit }
saveRDS(YN, file.path(ew, "phewas_y.rds"))
R <- R[order(R$set != "main", R$p), ]
write.csv(R, file.path(ew, "phewas_results.csv"), row.names = FALSE)

sink(file.path(ew, "phewas_summary.txt"))
cat("PheWAS(审核表 v1.2;KSAA → 表型;加权;协变量 年龄/性别/种族/教育/PIR)\n\n")
for (s in c("main", "supplementary", "exposure_module")) { k <- R$set == s
  cat(sprintf("[%s] 表型 %d(拟合失败 %d);模型:%s;FDR<5%% %d;Bonferroni %d;Li–Ji meff %.0f → %d;两相位复制 %d\n", s, sum(k), sum(k & is.na(R$p)),
    paste(names(table(R$model[k])), table(R$model[k]), collapse = " "), sum(R$fdr[k] < 0.05, na.rm = TRUE), sum(R$bonf[k] < 0.05, na.rm = TRUE), R$meff[k][1], sum(R$p_liji[k] < 0.05, na.rm = TRUE), sum(R$replicated[k], na.rm = TRUE))) }
cat(sprintf("\n主分析独立信号(FDR<5%% 者按 |r|>0.7 归并):%d 个;领头:%s\n", length(unique(na.omit(R$clump_lead))), paste(unique(na.omit(R$clump_lead)), collapse = " ")))
cat("\n主分析 FDR<5% 按健康类别:\n"); print(table(R$health_class[R$set == "main" & R$fdr < 0.05]))
cat("\n主分析前 30:\n"); k <- which(R$set == "main")[1:30]
for (i in k) cat(sprintf("  %-9s %-46s %-9s n=%-5d %s %s  P=%.1e FDR=%.1e  ph1 %.2f/%.1e ph2 %.2f/%.1e %s\n", R$variable[i], substr(R$phenotype[i], 1, 46), R$model[i], R$n[i],
  ifelse(is.na(R$est[i]), "   -  ", sprintf("%6.2f", R$est[i])), ifelse(is.na(R$lo[i]), "", sprintf("(%.2f-%.2f)", R$lo[i], R$hi[i])), R$p[i], R$fdr[i], R$est_phase1[i], R$p_phase1[i], R$est_phase2[i], R$p_phase2[i], ifelse(R$replicated[i], "REP", "")))
cat("\n拟合失败:\n"); print(R[is.na(R$p), c("variable", "phenotype", "model", "n", "error")], row.names = FALSE)
cat("\n回退/重编码说明:\n"); print(table(R$notes[grepl("singular|merged", R$notes)]))
sink()
cat(readLines(file.path(ew, "phewas_summary.txt")), sep = "\n")
