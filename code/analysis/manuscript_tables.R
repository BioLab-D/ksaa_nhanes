# 稿件用表(HTML 片段):table1(表 1)、tableS1(表 S1,骨干性能)、tableS3(表 S3,按 KSAA 四分位的特征);全部由 results/analysis/*.csv 与 predictions OOF 现算,不手抄;写到 results/analysis/tables/<name>.html。
# build_manuscript_imrad.py 以 <!--HTML:<name>--> 嵌入。稿件以 Office 文件为准,表注中的表号字符串只影响重新生成的 HTML。
res <- here::here("results"); v3 <- file.path(res, "analysis"); td <- file.path(v3, "tables")
dir.create(td, showWarnings = FALSE)
mi <- function(x) gsub("-", "−", x, fixed = TRUE)
f <- function(x, d = 2) ifelse(is.na(x), "", mi(sub("^-(0\\.0+)$", "\\1", formatC(round(x + sign(x) * 1e-9, d), format = "f", digits = d))))   # 常规四舍五入(同 fig_common.R)
cm <- function(x) formatC(round(x), format = "d", big.mark = ",")
sg <- function(x, d) { s <- sprintf(paste0("%+.", d, "f"), x); s[grepl("^[+-]0\\.0+$", s)] <- sprintf(paste0("%.", d, "f"), 0); mi(s) }   # 零不带正负号
th <- function(h) paste0("<thead><tr>", paste0(sprintf("<th%s>%s</th>", ifelse(seq_along(h) == 1, "", ' class="num"'), h), collapse = ""), "</tr></thead>")
row <- function(cells, ind = FALSE, cls = "") paste0(sprintf('<tr class="%s">', cls), paste0(sprintf('<td class="%s">%s</td>', c(if (ind) "ind" else "", rep("num", length(cells) - 1)), cells), collapse = ""), "</tr>")
hdr <- function(txt, k) sprintf('<tr class="hdr"><td colspan="%d">%s</td></tr>', k, txt)
put <- function(name, head, body, notes, thead = NULL) {   # thead 可传入自定义(双层)表头,缺省 th(head)
  html <- paste0('<div class="twrap"><table>', if (is.null(thead)) th(head) else thead, "<tbody>", paste0(body, collapse = ""), "</tbody></table></div>", '<div class="tnote">', paste0("<p>", notes, "</p>", collapse = ""), "</div>")
  html <- gsub("<(?=[0-9=\\s])", "&lt;", html, perl = TRUE)
  html <- gsub("[Kk]nee age gap", "KSAA", html)   # 主暴露显示为 KSAA(CSV 与匹配键不变);上一行转义标签外的 <(BMI <25、HDL <40 等)
  writeLines(html, file.path(td, paste0(name, ".html")), useBytes = TRUE); cat("->", name, "\n")
}
d <- readRDS(file.path(v3, "analysis_data.rds")); qc <- quantile(d$gap, 1:3 / 4)

# ---------- 骨干与集成的预测性能 ----------
A5 <- c("ResNet-18" = "resnet18", "ResNet-50" = "resnet50", "EfficientNet-B0" = "efficientnetb0", "ConvNeXt-Tiny" = "convnexttiny", "DenseNet-121" = "densenet121")
met <- function(p, a) { e <- p - a; c(MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)), CS5 = 100 * mean(abs(e) <= 5), r = cor(p, a), R2 = 1 - sum(e^2) / sum((a - mean(a))^2)) }
fm <- function(p) { m <- sapply(0:4, function(k) met(p[d$fold == k], d$age[d$fold == k]))
  c(sprintf("%s ± %s", f(mean(m["MAE", ]), 2), f(sd(m["MAE", ]), 2)), sprintf("%s ± %s", f(mean(m["RMSE", ]), 2), f(sd(m["RMSE", ]), 2)),
    sprintf("%s ± %s", f(mean(m["CS5", ]), 1), f(sd(m["CS5", ]), 1)), sprintf("%s ± %s", f(mean(m["r", ]), 3), f(sd(m["r", ]), 3))) }
body <- c()
for (nm in names(A5)) { x <- read.csv(file.path(res, "predictions", sprintf("m4_%s_oof.csv", A5[[nm]]))); p <- x$pred_age[match(d$SEQN, x$SEQN)]; stopifnot(!anyNA(p)); body <- c(body, row(c(nm, fm(p)))) }
body <- c(body, row(c("Weighted ensemble<sup>a</sup>", fm(d$pred_raw)), cls = "bold"))
# 指标均为未校正预测(偏倚校正用到实足年龄);表注报告各折 MAE 范围与训练折中位数基线。
me <- sapply(0:4, function(k) met(d$pred_raw[d$fold == k], d$age[d$fold == k])["MAE"])
base <- rep(NA_real_, nrow(d)); for (k in 0:4) base[d$fold == k] <- median(d$age[d$fold != k])
bm <- sapply(0:4, function(k) { i <- d$fold == k; mean(abs(base[i] - d$age[i])) })
put("tableS1", c("Model", "MAE, years", "RMSE, years", "CS (α = 5 years), %", "Pearson r"), body, c(
  sprintf("Mean ± SD across the five held-out folds of participant-level cross-fitting (n = %s). All metrics are for the uncorrected prediction, because the age-bias correction uses chronological age (Methods). The ensemble MAE ranged from %s to %s years across folds; assigning every participant the median age of the training folds gave an MAE of %s years.", cm(nrow(d)), f(min(me), 2), f(max(me), 2), f(mean(bm), 2)),
  "<sup>a</sup>Backbone predictions weighted by the inverse internal-validation MAE within the training folds; the primary prediction."))

# ---------- 按 KSAA 四分位的基线特征 ----------
T1 <- read.csv(file.path(v3, "quartile_characteristics.csv"), stringsAsFactors = FALSE); T1[is.na(T1)] <- ""
nq <- T1[T1$lab == "Participants, n", ]
# 不含分节小标题与死亡块;quartile_characteristics.csv 本身不动(表 1 的自检 gq() 要用其中的死亡行)。
S2DROP <- c("All-cause deaths", "Cardiovascular deaths", "Cancer deaths")
T1 <- T1[T1$type != "hdr" & !(T1$lab %in% S2DROP), ]
lab3 <- function(x) { x <- sub(", y$", ", years", x); x <- sub("Poverty-income", "Poverty–income", x, fixed = TRUE); sub("hip, wrist, or spine fracture", "hip, wrist or vertebral fracture", x, fixed = TRUE) }
body <- c()
for (i in seq_len(nrow(T1))) { r <- T1[i, ]
  body <- c(body, if (r$type == "hdr") hdr(r$lab, 6) else row(c(lab3(gsub("<sup>[^<]*</sup>", "", r$lab)), r$Q1, r$Q2, r$Q3, r$Q4, r$p), r$ind == 1)) }   # lab3() 统一单位与行名写法(quartile_characteristics.csv 写法不同)
body <- body[!grepl("Participants, n", body, fixed = TRUE)]
put("tableS3", c("Characteristic", sprintf("Q1 (n = %s)", nq$Q1), sprintf("Q2 (n = %s)", nq$Q2), sprintf("Q3 (n = %s)", nq$Q3), sprintf("Q4 (n = %s)", nq$Q4), "P"), body, c(
  # 列名用 P:多分类行为 χ² 检验,不是趋势检验;年龄呈 U 形,表注说明
  sprintf("Data are mean (SD) or n (%%), unweighted; definitions as in Table 1 and the Supplementary Methods, and whole-cohort values in Table 1. Percentages were calculated among participants with available data. Current alcohol use is ≥12 alcoholic drinks in the past 12 months and inactive is no leisure-time physical activity, as in Table 1. Quartile cut points %s, %s and %s years (Q4, knee older than expected for age). P values are from linear regression (continuous characteristics) or logistic regression (binary characteristics) on the quartile number, and from a χ² test for characteristics with more than two categories (race and ethnicity, education, smoking, alcohol use and leisure-time physical activity), for which no trend is implied. Mean age is highest in the lowest and in the highest quartile, so for age the regression on quartile number does not describe the U-shaped pattern.", f(qc[1]), f(qc[2]), f(qc[3]))))

# ---------- 表 1(fxf_vars.csv 由 fracture_function_vars.R 生成)----------
V <- read.csv(file.path(v3, "fxf_vars.csv")); d1 <- d; for (v in setdiff(names(V), "SEQN")) d1[[v]] <- V[[v]][match(d1$SEQN, V$SEQN)]
msd <- function(x, dig) { x <- x[!is.na(x)]; sprintf("%s (%s)", f(mean(x), dig), f(sd(x), dig)) }
pb <- function(x) { x <- x[!is.na(x)]; sprintf("%s (%s)", cm(sum(x == 1)), f(100 * mean(x == 1), 1)) }
mrow <- function(lab, x, dig, ind = FALSE) row(c(lab, cm(sum(!is.na(x))), msd(x, dig)), ind)
brow <- function(lab, x, ind = FALSE) row(c(lab, cm(sum(!is.na(x))), pb(x)), ind)
crow <- function(lab, x, lev, labs) { n <- sum(!is.na(x)); c(row(c(lab, cm(n), "")), unname(mapply(function(l, lb) row(c(lb, "", sprintf("%s (%s)", cm(sum(x == l, na.rm = TRUE)), f(100 * sum(x == l, na.rm = TRUE) / n, 1))), TRUE), lev, labs))) }
T1Q <- read.csv(file.path(v3, "quartile_characteristics.csv")); gq <- function(lab) { v <- T1Q$overall[T1Q$lab == lab]; stopifnot(length(v) == 1); v }
chk <- list(c("Age, years", msd(d1$age, 1)), c("Sex, female", pb(as.integer(d1$sex == 2))), c("Body mass index, kg/m²", msd(d1$bmi, 1)), c("Diabetes<sup>b</sup>", pb(d1$dm)),
  c("Total cholesterol, mmol/L", msd(d1$tc / 38.67, 2)), c("HDL cholesterol, mmol/L", msd(d1$hdl / 38.67, 2)), c("eGFR, mL/min/1.73 m²<sup>b</sup>", msd(d1$egfr, 1)), c("Phenotypic Age acceleration, years<sup>c</sup>", msd(d1$pa_accel, 1)),
  c("Total hip BMD, g/cm²<sup>d</sup>", msd(d1$tot_bmd, 3)), c("Radiographic knee OA<sup>e</sup>", pb(as.integer(d1$klmax >= 2))), c("Knee pain", pb(d1$pain)),
  c("Gait speed, m/s", msd(d1$gait, 2)), c("≥1 fall in the past 12 months", pb(d1$fall)), c("Knee age gap, mean (SD), y", msd(d1$gap, 2)), c("Cardiovascular deaths", pb(d1$e_cvd)))
for (z in chk) if (gq(z[1]) != z[2]) stop("表 1 与表 S2 总体列不一致: ", z[1], " ", gq(z[1]), " vs ", z[2])
rk <- summary(survival::survfit(survival::Surv(t / 12, 1 - e) ~ 1, data = d1))$table; qq <- quantile(survival::survfit(survival::Surv(t / 12, 1 - e) ~ 1, data = d1), c(.25, .5, .75))$quantile
# 表 1:有缺失的项目在脚注列出;教育用受教育年数(HFA8R,ewas_clean.rds 的清洗值);种族只展示非西班牙裔白人
CLN <- readRDS(file.path(v3, "phewas", "ewas_clean.rds")); stopifnot(identical(CLN$SEQN, d1$SEQN)); d1$edu_y <- CLN$x[["HFA8R"]]
NMISS <- list(); mrow <- function(lab, x, dig, ind = FALSE) { if (sum(!is.na(x)) < nrow(d1)) NMISS[[lab]] <<- sum(!is.na(x)); row(c(lab, msd(x, dig)), ind) }
brow <- function(lab, x, ind = FALSE) { if (sum(!is.na(x)) < nrow(d1)) NMISS[[lab]] <<- sum(!is.na(x)); row(c(lab, pb(x)), ind) }
body <- c(hdr("Demographics and socioeconomic status", 2),
  mrow("Age, years", d1$age, 1), brow("Female sex", as.integer(d1$sex == 2)),
  brow("Non-Hispanic White", as.integer(d1$race == 1)),
  mrow("Education, years", d1$edu_y, 1),
  brow("Family income below the poverty line", d1$poor),
  hdr("Lifestyle and cardiometabolic factors (Model 2 and 3 covariates)", 2),
  brow("Former or current smoking", as.integer(d1$smk3 %in% c("former", "current"))),
  brow("Current alcohol use", as.integer(d1$alc3 == "current")),
  brow("Physically inactive in leisure time", as.integer(d1$pa3 == "inactive")),
  mrow("Body mass index, kg/m²", d1$bmi, 1), brow("Hypertension", d1$htn), brow("Diabetes", d1$dm),
  mrow("Total cholesterol, mmol/L", d1$tc / 38.67, 2), mrow("HDL cholesterol, mmol/L", d1$hdl / 38.67, 2), mrow("eGFR, mL/min/1.73 m²", d1$egfr, 1),   # SI 单位(mg/dL ÷ 38.67)
  brow("History of cardiovascular disease", d1$cvd_hist), brow("History of cancer", d1$cancer_hist),
  hdr("Knee skeletal age and comparison measures", 2),
  row(c("KSAA, years, range (quartile cut points)", sprintf("%s to %s (%s, %s, %s)", f(min(d1$gap), 1), f(max(d1$gap), 1), f(quantile(d1$gap, .25), 2), f(quantile(d1$gap, .5), 2), f(quantile(d1$gap, .75), 2)))),
  mrow("Phenotypic Age acceleration, years", d1$pa_accel_raw, 1),   # 描述性用原始值(分析用 1/99 缩尾值)
  brow("Radiographic knee OA, K–L grade ≥2 (phase 2)", as.integer(d1$klmax >= 2)), brow("Knee pain", d1$pain),
  mrow("Total hip BMD, g/cm²", d1$tot_bmd, 3),
  hdr("Fracture, falls and physical function", 2),
  brow("Prior hip fracture (self-reported)", d1$fx_hip), brow("Any prior hip, wrist or vertebral fracture (self-reported)", d1$fx_any),
  brow("Fall in the past 12 months", d1$fall), mrow("Gait speed, m/s", d1$gait, 2),
  hdr("Mortality through 31 December 2019", 2),
  row(c("Follow-up time, years, median (IQR)", sprintf("%s (%s–%s)", f(median(d1$t / 12), 1), f(quantile(d1$t / 12, .25), 1), f(quantile(d1$t / 12, .75), 1)))),   # 实际随访(PERMTH_EXM/12)
  brow("All-cause deaths", d1$e), brow("Cardiovascular deaths", d1$e_cvd, TRUE), brow("Cancer deaths", d1$e_cancer, TRUE))
miss_note <- sprintf("Numbers with data where fewer than %s: %s.", cm(nrow(d1)), paste(sprintf("%s, %s", sub(",.*$", "", sub(" \\(.*$", "", names(NMISS))), cm(unlist(NMISS))), collapse = "; "))
put("table1", c("Characteristic", sprintf("Overall (n = %s)", cm(nrow(d1)))), body, c(
  "Data are mean (SD), n (%) or median (IQR), unweighted; percentages are of participants with data. Numbers with available data and survey-weighted estimates are given in Table S1, and characteristics by quartile of KSAA in Table S3.",
  sprintf("KSAA: residual of the out-of-fold predicted knee age on chronological age (SD %s years). Phenotypic Age acceleration is shown unwinsorized (winsorized at the 1st and 99th percentiles in all analyses). K–L grade: higher grade of the two knees, survey phase 2 only. Covariate definitions are given in the Supplementary Methods.", f(sd(d1$gap), 2)),
  "Follow-up time: from the examination to death or 31 December 2019. Cardiovascular deaths: heart disease or cerebrovascular disease as the underlying cause.",
  "Abbreviations: BMD, bone mineral density; eGFR, estimated glomerular filtration rate; HDL, high-density lipoprotein; IQR, interquartile range; K–L, Kellgren–Lawrence; KSAA, knee skeletal age acceleration; OA, osteoarthritis."))
