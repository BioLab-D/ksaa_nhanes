# 专家影像判读特征与近端股骨 DXA(偏相关 r 用于 3.4 节):与 KSAA 的偏相关,及各特征自身对全因死亡的 HR(Model 3 调整)。
#   判读特征:K–L、骨赘、硬化、软骨钙化、关节置换、影像 OA、膝痛;样本 = phase 2 判读完整(K–L、骨赘、硬化非缺失)且 Model 3 协变量完整者。
#   偏相关 = 两者各对 Model 3 协变量回归后的残差相关(二分类即点二列偏相关),CI 用 Fisher z(SE = 1/√(n − 3 − k),k 为协变量列数);
#   连续特征每 SD(在该特征的分析样本内标准化),二分类为有 vs 无;P 为 Wald,Bonferroni 按 5 项判读特征(FEATP)。
#   特征构造:每膝 4 区室之和,取双膝较重者;8/9 置缺失;膝痛同 analysis_data 的 pain(HAG22)。
# 输入 results/analysis/analysis_data.rds、results/work/xray_components.csv
# 输出 results/analysis/kl_dxa_partial_cor.csv
# 由 run_all.sh 运行
suppressMessages(library(survival))
res <- here::here("results"); v3 <- file.path(res, "analysis"); spw <- file.path(res, "work")
d <- readRDS(file.path(v3, "analysis_data.rds")); SDG <- sd(d$gap); d$z <- d$gap / SDG
X <- read.csv(file.path(spw, "xray_components.csv")); cl <- function(v) ifelse(v %in% c(8, 9), NA, v)
for (v in setdiff(names(X), "SEQN")) X[[v]] <- cl(X[[v]])
X$osteo <- pmax(rowSums(X[, c("XRPOMFR", "XRPOMTR", "XRPOLFR", "XRPOLTR")]), rowSums(X[, c("XRPOMFL", "XRPOMTL", "XRPOLFL", "XRPOLTL")]), na.rm = TRUE)
X$scler <- pmax(rowSums(X[, c("XRPSMFR", "XRPSMTR", "XRPSLFR", "XRPSLTR")]), rowSums(X[, c("XRPSMFL", "XRPSMTL", "XRPSLFL", "XRPSLTL")]), na.rm = TRUE)
X$kl <- pmax(X$XRPKLR, X$XRPKLL, na.rm = TRUE)
X$chondro <- as.integer(pmax(X$XRPCHOR, X$XRPCHOL, na.rm = TRUE) == 1)
X$tka <- as.integer(pmax(X$XRPJRR, X$XRPJRL, na.rm = TRUE) == 1)
d <- merge(d, X[, c("SEQN", "kl", "osteo", "scler", "chondro", "tka")], by = "SEQN")
d$roa <- as.integer(d$kl >= 2)
M3V <- c("age", "sex", "race", "edu3", "smk3", "alc3", "pa3", "bmi", "htn", "dm", "tc", "hdl", "egfr", "cvd_hist", "cancer_hist")
COV <- "age + factor(sex) + factor(race) + edu3 + smk3 + alc3 + pa3 + bmi + htn + dm + tc + hdl + egfr + cvd_hist + cancer_hist"
stopifnot(all(M3V %in% names(d)))
s <- d[!is.na(d$kl) & !is.na(d$osteo) & !is.na(d$scler) & complete.cases(d[, M3V]), ]
stopifnot(all(s$kl == s$klmax))   # 与 analysis_data 的 klmax 一致
cat(sprintf("样本:第 2 期判读完整且 Model 3 协变量完整 n = %d,死亡 %d\n", nrow(s), sum(s$e)))
FEAT <- c(kl = "K–L grade", osteo = "Osteophyte score", scler = "Subchondral sclerosis score", chondro = "Chondrocalcinosis",
          tka = "Joint replacement", roa = "Radiographic osteoarthritis (K–L ≥2)", pain = "Knee pain")
BIN <- c("chondro", "tka", "roa", "pain")
# Bonferroni 只计五项判读特征(FEATP);影像 OA(K–L 的二分类形式)与膝痛不计入。
FEATP <- c("kl", "osteo", "scler", "chondro", "tka")
pcor <- function(v, t) { k <- ncol(model.matrix(as.formula(paste("~", COV)), t)) - 1
  r <- cor(resid(lm(as.formula(paste(v, "~", COV)), t)), resid(lm(as.formula(paste("z ~", COV)), t))); se <- 1 / sqrt(nrow(t) - 3 - k)
  c(r = r, lo = tanh(atanh(r) - 1.96 * se), hi = tanh(atanh(r) + 1.96 * se)) }
R <- NULL
for (v in names(FEAT)) { t <- s[!is.na(s[[v]]), ]
  t$q <- if (v %in% BIN) t[[v]] else as.numeric(scale(t[[v]]))
  m <- coxph(as.formula(paste("Surv(t, e) ~ q +", COV)), t); b <- coef(m)[["q"]]; se <- sqrt(vcov(m)["q", "q"]); p <- 2 * pnorm(-abs(b / se))
  pr <- pcor(v, t)
  R <- rbind(R, data.frame(key = v, feature = FEAT[[v]], binary = v %in% BIN, n = nrow(t), deaths = sum(t$e), prev = if (v %in% BIN) sum(t[[v]]) else NA,
                           r = pr[["r"]], r_lo = pr[["lo"]], r_hi = pr[["hi"]], hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = p, p_bonf = if (v %in% FEATP) min(1, length(FEATP) * p) else NA)) }
mk <- coxph(as.formula(paste("Surv(t, e) ~ z +", COV)), s); b <- coef(mk)[["z"]]; se <- sqrt(vcov(mk)["z", "z"])
R <- rbind(R, data.frame(key = "ksaa", feature = "KSAA (per SD)", binary = FALSE, n = nrow(s), deaths = sum(s$e), prev = NA, r = NA, r_lo = NA, r_hi = NA,
                         hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), p_bonf = NA))
rownames(R) <- NULL
print(format(R[, -1], digits = 3))

# ---- DXA 块:股骨颈 BMD、全髋 BMD、全髋骨矿含量与骨面积(说明 BMD 偏低来自矿物少而非骨小) ----
#   全 8 哨兵按 q() 规则置缺失;全髋 BMD = BDPTOBMC / BDPTOARE。样本 = 四项均非缺失且 Model 3 协变量完整(两期)。
DX <- read.csv(file.path(spw, "dxa_all.csv"))
q <- function(x) { x <- suppressWarnings(as.numeric(x))
  for (k in 1:6) { s8 <- as.numeric(strrep("8", k)); lo <- x[!is.na(x) & x < s8]
    if (any(x == s8, na.rm = TRUE) && length(lo) && s8 >= max(lo) * 1.5) x[x == s8] <- NA }
  x }
for (v in setdiff(names(DX), "SEQN")) DX[[v]] <- q(DX[[v]])
DX$BDPTOBMD <- DX$BDPTOBMC / DX$BDPTOARE
g <- merge(readRDS(file.path(v3, "analysis_data.rds")), DX[, c("SEQN", "BDPFNBMD", "BDPTOBMD", "BDPTOBMC", "BDPTOARE")], by = "SEQN")
g$z <- g$gap / SDG
DXV <- c(BDPFNBMD = "Femoral neck BMD", BDPTOBMD = "Total hip BMD", BDPTOBMC = "Total hip bone mineral content", BDPTOARE = "Total hip bone area")
gs <- g[complete.cases(g[, c(names(DXV), M3V)]), ]
cat(sprintf("\nDXA 块样本:四项 DXA 均非缺失且 Model 3 协变量完整 n = %d,死亡 %d\n", nrow(gs), sum(gs$e)))
for (v in names(DXV)) { t <- gs; t$q <- as.numeric(scale(t[[v]]))
  m <- coxph(as.formula(paste("Surv(t, e) ~ q +", COV)), t); b <- coef(m)[["q"]]; se <- sqrt(vcov(m)["q", "q"]); pr <- pcor(v, t)
  R <- rbind(R, data.frame(key = v, feature = DXV[[v]], binary = FALSE, n = nrow(t), deaths = sum(t$e), prev = NA, r = pr[["r"]], r_lo = pr[["lo"]], r_hi = pr[["hi"]],
                           hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), p_bonf = NA)) }
mk <- coxph(as.formula(paste("Surv(t, e) ~ z +", COV)), gs); b <- coef(mk)[["z"]]; se <- sqrt(vcov(mk)["z", "z"])
R <- rbind(R, data.frame(key = "ksaa_dxa", feature = "KSAA (per SD)", binary = FALSE, n = nrow(gs), deaths = sum(gs$e), prev = NA, r = NA, r_lo = NA, r_hi = NA,
                         hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), p_bonf = NA))
R$block <- ifelse(R$key %in% c(names(DXV), "ksaa_dxa"), "dxa", "reading")
print(format(R[R$block == "dxa", c("feature", "n", "deaths", "r", "r_lo", "r_hi", "hr", "lo", "hi", "p")], digits = 3))
write.csv(R, file.path(v3, "kl_dxa_partial_cor.csv"), row.names = FALSE)
