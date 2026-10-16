# 加权主分析的公共设置(CDC 膝 / 手 X 光数据说明要求影像相关分析使用体检样本权重;加权为主分析,不加权为敏感性分析)。
# 做法:先在**全部 NHANES III 成人体检者**(WTPFEX6 > 0,8,805 + 8,900 人)上建立完整调查设计,再以研究子样本为 domain 取 subset;
#   仅第二期的分析(K–L 判读)用第二期权重 WTPFEX2 与 SDPPSU2 / SDPSTRA2,在全部第二期体检者上建设计后取 subset。
# 由 weighted.R(模型)与 weighted_boot.R(重复权重的 bootstrap)source。
suppressMessages({library(survival); library(survey); library(splines); library(duckdb)})
options(survey.lonely.psu = "adjust")
RES <- here::here("results"); V3 <- file.path(RES, "analysis"); SPW <- file.path(RES, "work")
W <- file.path(V3, "weighted"); dir.create(W, showWarnings = FALSE)
d <- readRDS(file.path(V3, "analysis_data.rds")); SDG <- sd(d$gap); d$z <- d$gap / SDG
QL <- paste0("Q", 1:4); d$q <- factor(d$q, QL)
d$qmed <- tapply(d$gap, d$q, median)[as.character(d$q)]
FX <- read.csv(file.path(V3, "fxf_vars.csv")); d <- merge(d, FX, by = "SEQN", all.x = TRUE)
d$cause <- factor(ifelse(d$e == 0, "censor", ifelse(d$e_cvd == 1, "cvd", ifelse(d$e_cancer == 1, "cancer", "other"))), c("censor", "cvd", "cancer", "other"))

con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
A <- dbGetQuery(con, "SELECT SEQN, SDPPHASE AS phase_all, WTPFEX6 AS W6, SDPPSU6 AS P6, SDPSTRA6 AS S6, WTPFEX2 AS W2, SDPPSU2 AS P2, SDPSTRA2 AS S2,
                      HSAGEIR, HSSEX, DMARETHN, HAB1, HAH1, HAH9, HFA8R FROM nhanes3_design_all WHERE WTPFEX6 > 0")
MORT <- dbGetQuery(con, "SELECT SEQN, MORTSTAT, PERMTH_EXM FROM nhanes3_design_all")
dbDisconnect(con, shutdown = TRUE)
drop <- c("WTPFEX6", "WTPFEX2", "SDPPSU6", "SDPSTRA6", "SDPPSU2", "SDPSTRA2")
F <- merge(A, d[, setdiff(names(d), drop)], by = "SEQN", all.x = TRUE)
F$inc <- F$SEQN %in% d$SEQN
stopifnot(sum(F$inc) == nrow(d), nrow(A) == 17705)
chk <- merge(d[, c("SEQN", "WTPFEX6", "SDPPSU6", "SDPSTRA6")], A[, c("SEQN", "W6", "P6", "S6")], by = "SEQN")
stopifnot(nrow(chk) == nrow(d), max(abs(chk$WTPFEX6 - chk$W6)) < 1e-6, all(chk$SDPPSU6 == chk$P6), all(chk$SDPSTRA6 == chk$S6))

M3V <- c("age", "sex", "race", "edu3", "smk3", "alc3", "pa3", "bmi", "htn", "dm", "tc", "hdl", "egfr", "cvd_hist", "cancer_hist")
COV <- c("age", "factor(sex)", "factor(race)", "edu3", "smk3", "alc3", "pa3", "bmi", "htn", "dm", "tc", "hdl", "egfr", "cvd_hist", "cancer_hist")
MODS <- c(m1 = "age + factor(sex)", m2 = "age + factor(sex) + factor(race) + edu3 + smk3 + alc3 + pa3 + bmi", m3 = paste(COV, collapse = " + "))
EV <- c(allcause = "e", cvd = "e_cvd", cancer = "e_cancer")
F$cc <- F$inc & complete.cases(F[, M3V])
stopifnot(sum(F$cc) == 4347)

D6 <- svydesign(id = ~P6, strata = ~S6, weights = ~W6, nest = TRUE, data = F)
F2 <- F[F$phase_all == 2 & !is.na(F$W2) & F$W2 > 0, ]
D2 <- svydesign(id = ~P2, strata = ~S2, weights = ~W2, nest = TRUE, data = F2)
dom6 <- function(x) subset(D6, x)                       # x: 与 F 等长的逻辑向量
dom2 <- function(x) subset(D2, x[match(F2$SEQN, F$SEQN)])

ci3 <- function(b, se, k = 1) c(hr = exp(k * b), lo = exp(k * (b - 1.96 * se)), hi = exp(k * (b + 1.96 * se)))
wcox <- function(des, rhs, term = "z", ev = "e", k = 1) {   # 加权 Cox;返回该项 HR / CI / P 与 domain 的人数、事件数
  m <- svycoxph(as.formula(sprintf("Surv(t, %s) ~ %s", ev, rhs)), design = des)
  b <- coef(m)[[term]]; se <- sqrt(vcov(m)[term, term]); dd <- model.frame(des)
  c(ci3(b, se, k), p = 2 * pnorm(-abs(b / se)), n = m$n, events = m$nevent)
}
f2 <- function(x, k = 2) formatC(round(x + sign(x) * 1e-9, k), format = "f", digits = k)
fci <- function(h, l, u) sprintf("%s (%s–%s)", f2(h), f2(l), f2(u))
