# 联合分层:KSAA 四分位 × 对照指标(PhenoAge 加速、全髋 BMD、K–L)交叉分格,格子填相对参照格的 HR(生存结局用 HR 而非 OR);
#   标准化 10 年绝对风险(Cox g-formula)一并写入 CSV。分档用四分位(每格约 280 人);K–L 合并为 0–1 / ≥2(高级别格子人数太少)。
#   心血管、癌症版为死因别 Cox(其他死因删失),1 − S 不是累积发生率,故 risk 记 NA。
# 输入 results/analysis/analysis_data.rds;输出 results/analysis/{joint_cells.csv, joint_within.csv}
source(here::here("code", "analysis", "fig_common.R"))
d <- readRDS(file.path(V3, "analysis_data.rds"))
TIME <- 120; K <- 4                       # 月 = 10 年;四分位
QL <- paste0("Q", 1:K)

qcut <- function(x, k = K) cut(x, quantile(x, 0:k / k, na.rm = TRUE), include.lowest = TRUE, labels = paste0("Q", 1:k))

# Cox g-formula:把全样本依次设成每个格子,平均预测生存,1 - surv 即该格的标准化 10 年死亡风险
std_cell <- function(dat, adj = "age + factor(sex)", time = TIME, ev = "e") {
  fit <- coxph(as.formula(sprintf("Surv(t, %s) ~ cell + %s", ev, adj)), dat)
  levs <- levels(dat$cell)
  risk <- if (ev != "e") setNames(rep(NA_real_, length(levs)), levs) else sapply(levs, function(l) { nd <- dat; nd$cell <- factor(l, levs)
    1 - rowMeans(summary(survfit(fit, newdata = nd), times = time)$surv) })
  cf <- summary(fit)$coefficients
  hr <- hl <- hh <- setNames(rep(1, length(levs)), levs)   # 参照格三者均为 1,作图时标为 ref
  for (l in levs[-1]) { rn <- paste0("cell", l)
    if (rn %in% rownames(cf)) { b <- cf[rn, "coef"]; se <- cf[rn, "se(coef)"]
      hr[l] <- exp(b); hl[l] <- exp(b - 1.96 * se); hh[l] <- exp(b + 1.96 * se) } }
  list(risk = risk, hr = hr, hr_lo = hl, hr_hi = hh)
}

# y = KSAA 分档,x = 对照指标(已定向为越右越差);参照格为两者均最低的 Q1 × 首档
panel_data <- function(dat, xlab_levels, ref_x, tag, ev = "e", ref_y = "Q1") {
  dat <- dat[!is.na(dat$xq) & !is.na(dat$gap), ]
  dat$yq <- qcut(dat$gap)
  dat$cell <- factor(paste(dat$yq, dat$xq, sep = "|"))
  dat$cell <- relevel(dat$cell, ref = paste(ref_y, ref_x, sep = "|"))
  o <- std_cell(dat, ev = ev)
  g <- do.call(rbind, lapply(names(o$risk), function(l) {
    p <- strsplit(l, "|", fixed = TRUE)[[1]]
    data.frame(panel = tag, y = p[1], x = p[2], risk = o$risk[[l]],
               hr = o$hr[[l]], hr_lo = o$hr_lo[[l]], hr_hi = o$hr_hi[[l]],
               ref = l == levels(dat$cell)[1], n = sum(dat$cell == l), deaths = sum(dat$cell == l & dat[[ev]] == 1))
  }))
  g$y <- factor(g$y, QL); g$x <- factor(g$x, xlab_levels)
  sd_gap <- sd(dat$gap)
  within <- do.call(rbind, lapply(xlab_levels, function(l) {
    s <- dat[dat$xq == l, ]
    m <- summary(coxph(as.formula(sprintf("Surv(t, %s) ~ I(gap / sd_gap) + age + factor(sex)", ev)), s))$coefficients[1, ]
    data.frame(panel = tag, stratum = l, n = nrow(s), deaths = sum(s[[ev]] == 1),
               hr = exp(m["coef"]), lo = exp(m["coef"] - 1.96 * m["se(coef)"]), hi = exp(m["coef"] + 1.96 * m["se(coef)"]))
  }))
  list(grid = g, within = within, n = nrow(dat))
}

dM <- d[!is.na(d$pa_accel), ];  dM$xq <- qcut(dM$pa_accel)          # 主图:PhenoAge 加速(越高越差)
dB <- d[!is.na(d$tot_bmd), ];   dB$xq <- qcut(-dB$tot_bmd)          # 补充:全髋 BMD(取负,Q4 = 最低 BMD)
KL <- c("0–1", "≥2")                                       # 补充:K–L 两档(phase 2)
dC <- d[!is.na(d$klmax), ];     dC$xq <- factor(ifelse(dC$klmax <= 1, KL[1], KL[2]), KL)
M <- panel_data(dM, QL, "Q1", "main"); B <- panel_data(dB, QL, "Q1", "bmd"); C <- panel_data(dC, KL, KL[1], "kl")
MV <- panel_data(dM, QL, "Q1", "main_cvd", ev = "e_cvd")
MK <- panel_data(dM, QL, "Q1", "main_cancer", ev = "e_cancer")  # 癌症死亡版(备选)

G <- rbind(M$grid, B$grid, C$grid, MV$grid, MK$grid); W <- rbind(M$within, B$within, C$within, MV$within, MK$within)
write.csv(G, file.path(V3, "joint_cells.csv"), row.names = FALSE)
write.csv(W, file.path(V3, "joint_within.csv"), row.names = FALSE)

ov <- function(dat) 100 * mean(dat$xq[qcut(dat$gap) == "Q4"] == "Q4", na.rm = TRUE)
cat(sprintf("\n重叠度(KSAA 最高四分位中同时为对照最差四分位的百分比;完全独立 = 25%%)\n  PhenoAge 加速 %.1f%%\n  全髋 BMD %.1f%%\n", ov(dM), ov(dB)))
cat("\n各档之内 KSAA 每 SD 的 HR(年龄、性别调整)\n")
for (i in seq_len(nrow(W))) with(W[i, ], cat(sprintf("  %-5s %-6s n=%-5d 死亡 %-5d HR %s\n", panel, stratum, n, deaths, ci(hr, lo, hi))))
cat(sprintf("\n样本:main n=%d  bmd n=%d  kl n=%d\n", M$n, B$n, C$n))
