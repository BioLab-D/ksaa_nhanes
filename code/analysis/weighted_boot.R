# 加权主分析中需要重复权重的部分:图 2 标化绝对风险与曲线、图 3B 与表 S6 的区分度(Harrell C、R²_D)
# 重复权重:在全部体检者的设计上用 Rao–Wu 重缩放 bootstrap(survey::as.svrepdesign, type = "subbootstrap")生成,再取研究样本的行;
#   每个重复用对应的重复权重重新拟合全部模型,95% CI 取分位数。点估计用原始权重。第 2 期(K–L)用第 2 期设计的重复权重。
# 由 run_all.sh 运行
source(here::here("code", "analysis", "weighted_common.R"))
suppressMessages(library(parallel)); NC <- max(1, detectCores() - 1)
set.seed(20260929); R_RISK <- 500; R_C <- 200
RW6 <- weights(as.svrepdesign(D6, type = "subbootstrap", replicates = R_RISK), type = "analysis")
RW2 <- weights(as.svrepdesign(D2, type = "subbootstrap", replicates = R_C), type = "analysis")

# ---------------- 图 2:标化累积发生率(g-formula,权重进入 Cox 拟合与人群平均)----------------
co <- F[F$inc, ]; wco <- co$W6; rco <- RW6[F$inc, ]
stdrisk <- function(dat, w, grid) {
  k <- w > 0; dat <- dat[k, ]; w <- w[k]; dat$w <- w   # 重复权重为 0 的行(该 PSU 未被抽中)不参与拟合,也不贡献人群平均
  step_at <- function(bh, g) { i <- findInterval(g, bh$time); c(0, bh$hazard)[i + 1] }
  cs <- function(ev) coxph(as.formula(sprintf("Surv(t, cause == '%s') ~ q + age + factor(sex)", ev)), dat, weights = w)
  FIT <- list(all = coxph(Surv(t, e) ~ q + age + factor(sex), dat, weights = w), cvd = cs("cvd"), cancer = cs("cancer"), other = cs("other"))
  H0 <- lapply(FIT, function(m) step_at(basehaz(m, centered = FALSE), grid)); wn <- w / sum(w)
  out <- lapply(QL, function(l) { nd <- dat; nd$q <- factor(l, QL)
    r <- lapply(FIT, function(m) exp(predict(m, newdata = nd, type = "lp", reference = "zero")))
    dH <- lapply(H0[c("cvd", "cancer", "other")], function(h) diff(c(0, h)))
    L <- lapply(names(dH), function(k) outer(r[[k]], dH[[k]])); names(L) <- names(dH); dL <- L$cvd + L$cancer + L$other
    S <- exp(-t(apply(dL, 1, cumsum))); pd <- cbind(1, S[, -ncol(S)]) * (1 - exp(-dL))
    cif <- function(k) cumsum(colSums(wn * pd * ifelse(dL > 0, L[[k]] / dL, 0)))
    data.frame(q = l, month = grid, allcause = 1 - colSums(wn * exp(-outer(r$all, H0$all))), cvd = cif("cvd"), cancer = cif("cancer")) })
  do.call(rbind, out) }
CUR <- stdrisk(co, wco, 0:360); write.csv(CUR, file.path(W, "figure2_std_curves.csv"), row.names = FALSE)
est10 <- CUR[CUR$month == 120, ]
bt <- simplify2array(mclapply(seq_len(R_RISK), function(b) { x <- stdrisk(co, rco[, b], 0:120); x <- x[x$month == 120, ]; unlist(x[x$q == "Q4", c("allcause", "cvd", "cancer")] - x[x$q == "Q1", c("allcause", "cvd", "cancer")]) }, mc.cores = NC))
# Model 3 全因(完整病例)
c3 <- F[F$cc, ]; w3 <- c3$W6; r3 <- RW6[F$cc, ]
std3 <- function(w) { k <- w > 0; dat <- c3[k, ]; w <- w[k]; dat$w <- w; fit <- coxph(as.formula(paste("Surv(t, e) ~ q +", MODS[["m3"]])), dat, weights = w)
  bh <- basehaz(fit, centered = FALSE); H <- c(0, bh$hazard)[findInterval(120, bh$time) + 1]
  sapply(QL, function(l) { nd <- dat; nd$q <- factor(l, QL); sum(w * (1 - exp(-exp(predict(fit, newdata = nd, type = "lp", reference = "zero")) * H))) / sum(w) }) }
e3 <- std3(w3); bt3 <- unlist(mclapply(seq_len(R_RISK), function(b) { x <- std3(r3[, b]); x[["Q4"]] - x[["Q1"]] }, mc.cores = NC))
AR <- rbind(do.call(rbind, lapply(c("allcause", "cvd", "cancer"), function(o) data.frame(model = "M1 (full cohort)", outcome = o,
  Q1 = est10[[o]][1], Q2 = est10[[o]][2], Q3 = est10[[o]][3], Q4 = est10[[o]][4], rd = est10[[o]][4] - est10[[o]][1], rd_lo = quantile(bt[o, ], .025), rd_hi = quantile(bt[o, ], .975)))),
  data.frame(model = "M3 (complete cases)", outcome = "allcause", Q1 = e3[[1]], Q2 = e3[[2]], Q3 = e3[[3]], Q4 = e3[[4]], rd = e3[[4]] - e3[[1]], rd_lo = quantile(bt3, .025), rd_hi = quantile(bt3, .975)))
write.csv(AR, file.path(W, "figure2_absrisk.csv"), row.names = FALSE)
cat("=========== 图 2 标化 10 年风险(加权)===========\n")
for (i in seq_len(nrow(AR))) with(AR[i, ], cat(sprintf("  %-20s %-8s Q1 %.1f Q2 %.1f Q3 %.1f Q4 %.1f | Q4−Q1 %.1f (%.1f 到 %.1f)\n", model, outcome, 100 * Q1, 100 * Q2, 100 * Q3, 100 * Q4, 100 * rd, 100 * rd_lo, 100 * rd_hi)))
c30 <- CUR[CUR$month == 360, ]; cat(sprintf("  30 年:癌症 %s;心血管 %s;全因 %s\n", paste(sprintf("%.1f", 100 * c30$cancer), collapse = "/"), paste(sprintf("%.1f", 100 * c30$cvd), collapse = "/"), paste(sprintf("%.1f", 100 * c30$allcause), collapse = "/")))
PT <- sapply(EV, function(ev) wcox(dom6(F$inc), "qmed + age + factor(sex)", term = "qmed", ev = ev)[["p"]])
write.csv(data.frame(outcome = names(EV), p_trend = PT), file.path(W, "figure2_ptrend.csv"), row.names = FALSE)
cat(sprintf("  P for trend(全队列,年龄 + 性别,加权):%s\n", paste(sprintf("%s %.2g", names(EV), PT), collapse = "; ")))

# ---------------- 区分度:加权 Harrell C 与加权 R²_D ----------------
KAP <- sqrt(8 / pi)
wC <- function(tt, ee, lp, w) concordance(Surv(tt, ee) ~ lp, weights = w, reverse = TRUE)$concordance
wR2 <- function(tt, ee, lp, w) { o <- order(lp); cw <- cumsum(w[o]); p <- numeric(length(lp)); p[o] <- (cw - w[o] / 2) / sum(w)
  s <- qnorm(p) / KAP; D <- coef(coxph(Surv(tt, ee) ~ s, weights = w)); r <- unname(D / KAP)^2; r / (pi^2 / 6 + r) }
# 检验实现:不加权时与 survival::royston 的 R.D 一致
sp <- F$inc & !is.na(F$pa_accel); s0 <- F[sp, ]; m0 <- coxph(Surv(t, e) ~ age + factor(sex) + pa_accel, s0)
cat(sprintf("\n自检:不加权 R²_D 本实现 %.4f,royston %.4f;C 本实现 %.4f,concordance %.4f\n", wR2(s0$t, s0$e, predict(m0), rep(1, nrow(s0))), royston(m0)[["R.D"]],
  wC(s0$t, s0$e, predict(m0), rep(1, nrow(s0))), concordance(m0)$concordance))
disc <- function(dat, w, ev, base, add) { k <- w > 0; dat <- dat[k, ]; w <- w[k]; dat$w <- w / mean(w)
  m0 <- coxph(as.formula(sprintf("Surv(t, %s) ~ %s", ev, base)), dat, weights = w); m1 <- coxph(as.formula(sprintf("Surv(t, %s) ~ %s + %s", ev, base, add)), dat, weights = w)
  l0 <- predict(m0); l1 <- predict(m1)
  c(C0 = wC(dat$t, dat[[ev]], l0, dat$w), C1 = wC(dat$t, dat[[ev]], l1, dat$w), R0 = wR2(dat$t, dat[[ev]], l0, dat$w), R1 = wR2(dat$t, dat[[ev]], l1, dat$w)) }
summ <- function(pt, bs, lab) { d <- rbind(dC = bs["C1", ] - bs["C0", ], relC = 100 * (bs["C1", ] / bs["C0", ] - 1), relR = 100 * (bs["R1", ] / bs["R0", ] - 1))
  data.frame(block = lab, C0 = pt[["C0"]], C1 = pt[["C1"]], R0 = pt[["R0"]], R1 = pt[["R1"]], dC = pt[["C1"]] - pt[["C0"]], dC_lo = quantile(d["dC", ], .025), dC_hi = quantile(d["dC", ], .975),
    relC = 100 * (pt[["C1"]] / pt[["C0"]] - 1), relC_lo = quantile(d["relC", ], .025), relC_hi = quantile(d["relC", ], .975),
    relR2 = 100 * (pt[["R1"]] / pt[["R0"]] - 1), relR2_lo = quantile(d["relR", ], .025), relR2_hi = quantile(d["relR", ], .975),
    C0_lo = quantile(bs["C0", ], .025), C0_hi = quantile(bs["C0", ], .975), C1_lo = quantile(bs["C1", ], .025), C1_hi = quantile(bs["C1", ], .975),
    R0_lo = quantile(bs["R0", ], .025), R0_hi = quantile(bs["R0", ], .975), R1_lo = quantile(bs["R1", ], .025), R1_hi = quantile(bs["R1", ], .975)) }   # 各模型自身的 CI(图 3B 画区间用)
cat("\n=========== 图 3B(加权;年龄 + 性别 + PhenoAge 加速,再加 KSAA)===========\n")
rs <- RW6[sp, 1:R_C]; CB <- NULL
for (o in names(EV)) { ev <- EV[[o]]; pt <- disc(s0, s0$W6, ev, "age + factor(sex) + pa_accel", "gap")
  bs <- simplify2array(mclapply(seq_len(R_C), function(b) disc(s0, rs[, b], ev, "age + factor(sex) + pa_accel", "gap"), mc.cores = NC))
  CB <- rbind(CB, cbind(outcome = o, summ(pt, bs, "fig3b"))) }
write.csv(CB, file.path(W, "figure3b_cindex.csv"), row.names = FALSE)
for (i in seq_len(nrow(CB))) with(CB[i, ], cat(sprintf("  %-8s C %.3f -> %.3f  ΔC %.4f (%.4f, %.4f)  相对 %.2f%% (%.2f, %.2f) | R²_D %.3f -> %.3f 相对 %.1f%% (%.1f, %.1f)\n",
  outcome, C0, C1, dC, dC_lo, dC_hi, relC, relC_lo, relC_hi, R0, R1, relR2, relR2_lo, relR2_hi)))

cat("\n=========== 表 S6(加权;全因死亡,嵌套模型)===========\n")
M3 <- MODS[["m3"]]; S6 <- NULL
sb <- F$cc & !is.na(F$pa_accel); x <- F[sb, ]; rb <- RW6[sb, 1:R_C]
st <- function(dat, w, base, add) { dat$w <- w; coxph(as.formula(sprintf("Surv(t, e) ~ %s + %s", base, add)), dat, weights = w) }
pt1 <- disc(x, x$W6, "e", M3, "pa_accel"); bs1 <- simplify2array(mclapply(seq_len(R_C), function(b) disc(x, rb[, b], "e", M3, "pa_accel"), mc.cores = NC))
pt2 <- disc(x, x$W6, "e", paste(M3, "+ pa_accel"), "gap"); bs2 <- simplify2array(mclapply(seq_len(R_C), function(b) disc(x, rb[, b], "e", paste(M3, "+ pa_accel"), "gap"), mc.cores = NC))
pw <- function(des, rhs, term) { m <- svycoxph(as.formula(paste("Surv(t, e) ~", rhs)), design = des); regTermTest(m, as.formula(paste("~", term)), df = Inf)$p }   # 卡方参照(与 CI 一致):第 2 期设计仅 23 个自由度,F 检验分母自由度只剩约 1
S6 <- rbind(S6, cbind(summ(pt1, bs1, "Beyond PhenoAge: + Phenotypic Age acceleration"), n = nrow(x), events = sum(x$e), p_wald = pw(dom6(sb), paste(M3, "+ pa_accel"), "pa_accel")),
               cbind(summ(pt2, bs2, "Beyond PhenoAge: + KSAA"), n = nrow(x), events = sum(x$e), p_wald = pw(dom6(sb), paste(M3, "+ pa_accel + gap"), "gap")))
s2 <- F2$cc & F2$phase_all == 2 & !is.na(F2$klmax); y <- F2[s2, ]; r2 <- RW2[s2, 1:R_C]
pt3 <- disc(y, y$W2, "e", M3, "klmax"); bs3 <- simplify2array(mclapply(seq_len(R_C), function(b) disc(y, r2[, b], "e", M3, "klmax"), mc.cores = NC))
pt4 <- disc(y, y$W2, "e", paste(M3, "+ klmax"), "gap"); bs4 <- simplify2array(mclapply(seq_len(R_C), function(b) disc(y, r2[, b], "e", paste(M3, "+ klmax"), "gap"), mc.cores = NC))
D2s <- subset(D2, s2)
S6 <- rbind(S6, cbind(summ(pt3, bs3, "Beyond K–L (phase 2): + K–L grade"), n = nrow(y), events = sum(y$e), p_wald = pw(D2s, paste(M3, "+ klmax"), "klmax")),
               cbind(summ(pt4, bs4, "Beyond K–L (phase 2): + KSAA"), n = nrow(y), events = sum(y$e), p_wald = pw(D2s, paste(M3, "+ klmax + gap"), "gap")))
write.csv(S6, file.path(W, "tableS6.csv"), row.names = FALSE)
for (i in seq_len(nrow(S6))) with(S6[i, ], cat(sprintf("  %-48s %d/%d  C %.3f -> %.3f  ΔC %.3f (%.3f, %.3f)  Wald P %.2g\n", block, events, n, C0, C1, dC, dC_lo, dC_hi, p_wald)))
cat("\n-> results/analysis/weighted/\n")
