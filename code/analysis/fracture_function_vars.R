# 表 1、表 S10(weighted.R / weighted_common.R)与 phewas.R 所需的骨折、跌倒及身体功能结局变量
#  变量定义:骨折 HAG5A–C(8/9 缺失)、低能量 HAG10A–F、跌倒 HAG14(analysis_data 已去 88/99)、活动困难 HAH ≥2、串联站立未能坚持 10 秒;计数有断言核对
# 输入 results/analysis/analysis_data.rds、results/work/{fracture,physfun}.csv;输出 results/analysis/fxf_vars.csv
# 由 run_all.sh 运行
res <- here::here("results"); out <- file.path(res, "analysis"); sp <- file.path(res, "work")
d <- readRDS(file.path(out, "analysis_data.rds"))
FR <- read.csv(file.path(sp, "fracture.csv")); PF <- read.csv(file.path(sp, "physfun.csv")); stopifnot(!anyDuplicated(FR$SEQN), !anyDuplicated(PF$SEQN))
g <- function(D, v) suppressWarnings(as.numeric(D[[v]][match(d$SEQN, D$SEQN)]))
yn <- function(x) ifelse(x %in% c(8, 9), NA, as.integer(x == 1))
same <- function(a, b) all((!is.na(a) & !is.na(b) & a == b) | (is.na(a) & is.na(b)))
hip <- yn(g(FR, "HAG5A")); d$fx_wri <- yn(g(FR, "HAG5B")); d$fx_spi <- yn(g(FR, "HAG5C"))
stopifnot(same(hip, d$fx_hip))
fx <- cbind(hip, d$fx_wri, d$fx_spi); any_ <- as.integer(rowSums(fx, na.rm = TRUE) > 0); any_[rowSums(!is.na(fx)) == 0] <- NA; stopifnot(same(any_, d$fx_any))
# 低能量骨折:HAG10A–F 为各处骨折的受伤机制(1 = 从站立高度或更低处跌倒),只对 50 岁及以后的骨折(HAG9)询问;8/9 = 不知道/拒答。
#   任一处为 1 → 阳性;否则若有机制 8/9、或 50 岁后骨折而机制缺失 → 不可判定(NA);其余 → 阴性。
#   HAG10D–F 在库中为 BOOLEAN,转数值后即原码 1。
num <- function(x) if (is.logical(x)) as.numeric(x) else suppressWarnings(as.numeric(x))
idx <- match(d$SEQN, FR$SEQN)
MEC <- sapply(paste0("HAG10", LETTERS[1:6]), function(v) num(FR[[v]])[idx])
AGE <- sapply(paste0("HAG9", LETTERS[1:6], "R"), function(v) { a <- num(FR[[v]])[idx]; a[a >= 500] <- NA; a })   # 555 / 888 等特殊码不是年龄
n_low <- rowSums(MEC == 1 & !is.na(MEC))
n_unk <- rowSums((MEC %in% c(8, 9) & !is.na(MEC)) | (is.na(MEC) & !is.na(AGE) & AGE >= 50))
d$fx_frag <- ifelse(is.na(d$fx_any), NA_integer_, ifelse(d$fx_any == 0, 0L, ifelse(n_low > 0, 1L, ifelse(n_unk > 0, NA_integer_, 0L))))
cat(sprintf("  低能量骨折:阳性 %d;机制不可判定转为缺失 %d 人(明确 8/9 %d;50 岁及以后骨折但机制缺失 %d)\n",
            sum(d$fx_frag == 1, na.rm = TRUE), sum(is.na(d$fx_frag)) - sum(is.na(d$fx_any)),
            sum(d$fx_any == 1 & rowSums(MEC %in% c(8, 9) & !is.na(MEC)) > 0, na.rm = TRUE),
            sum(d$fx_any == 1 & rowSums(is.na(MEC) & !is.na(AGE) & AGE >= 50) > 0, na.rm = TRUE)))
for (v in c("HAH1", "HAH2", "HAH8", "HAH9")) { x <- g(PF, v); x[x %in% c(8, 9)] <- NA; d[[paste0("dif_", v)]] <- ifelse(x >= 2, 1L, ifelse(x == 1, 0L, NA)) }
d$tand_fail <- as.integer(d$tandem < 10)
cnt <- function(x) c(sum(x == 1, na.rm = TRUE), sum(!is.na(x)))
EXP <- list(fx_hip = c(143, 4875), fx_wri = c(422, 4875), fx_spi = c(103, 4874), fx_frag = c(245, 4838), fx_any = c(614, 4875), fall = c(1113, 4852),   # fx_frag 分母 4838:机制不可判定者为缺失
            dif_HAH1 = c(1401, 4720), dif_HAH2 = c(1341, 4716), dif_HAH8 = c(344, 4855), dif_HAH9 = c(1026, 4853))
# 活动困难四项的参照计数为 BMI 非缺失者(少 16 人),故在该子集中比对
for (v in names(EXP)) { w <- if (grepl("^dif_", v)) !is.na(d$bmi) else rep(TRUE, nrow(d)); k <- cnt(d[[v]][w]); kr <- cnt(d[[v]])
  cat(sprintf("  %-9s 原始 %5d / %5d;比对口径 %5d / %5d(旧表 4:%d / %d)\n", v, kr[1], kr[2], k[1], k[2], EXP[[v]][1], EXP[[v]][2])); stopifnot(all(k == EXP[[v]])) }
cat(sprintf("  tand_fail %5d / %5d;gait n %d\n", cnt(d$tand_fail)[1], cnt(d$tand_fail)[2], sum(!is.na(d$gait))))
write.csv(d[, c("SEQN", "fx_wri", "fx_spi", "fx_frag", "tand_fail", "dif_HAH1", "dif_HAH2", "dif_HAH8", "dif_HAH9")], file.path(out, "fxf_vars.csv"), row.names = FALSE)
