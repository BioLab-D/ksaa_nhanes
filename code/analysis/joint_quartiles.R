# Joint stratification: KSAA quartile x comparison measure (PhenoAge acceleration, total hip BMD, K–L) cells, each filled with the HR vs the reference cell (HR rather than OR for survival);
#   standardized 10-year absolute risks (Cox g-formula) are also written to the CSV. Quartiles (about 280 per cell); K–L collapsed to 0–1 / ≥2 (too few in higher grades).
#   Cardiovascular and cancer versions use cause-specific Cox (other deaths censored); 1 − S is not a cumulative incidence, so risk is NA.
# Input results/analysis/analysis_data.rds; output results/analysis/{joint_cells.csv, joint_within.csv}
source(here::here("code", "analysis", "fig_common.R"))
d <- readRDS(file.path(V3, "analysis_data.rds"))
TIME <- 120; K <- 4                       # months = 10 years; quartiles
QL <- paste0("Q", 1:K)

qcut <- function(x, k = K) cut(x, quantile(x, 0:k / k, na.rm = TRUE), include.lowest = TRUE, labels = paste0("Q", 1:k))

# Cox g-formula: set the whole sample to each cell in turn and average predicted survival; 1 - surv is that cell's standardized 10-year mortality risk
std_cell <- function(dat, adj = "age + factor(sex)", time = TIME, ev = "e") {
  fit <- coxph(as.formula(sprintf("Surv(t, %s) ~ cell + %s", ev, adj)), dat)
  levs <- levels(dat$cell)
  risk <- if (ev != "e") setNames(rep(NA_real_, length(levs)), levs) else sapply(levs, function(l) { nd <- dat; nd$cell <- factor(l, levs)
    1 - rowMeans(summary(survfit(fit, newdata = nd), times = time)$surv) })
  cf <- summary(fit)$coefficients
  hr <- hl <- hh <- setNames(rep(1, length(levs)), levs)   # all three are 1 in the reference cell, labelled ref in plots
  for (l in levs[-1]) { rn <- paste0("cell", l)
    if (rn %in% rownames(cf)) { b <- cf[rn, "coef"]; se <- cf[rn, "se(coef)"]
      hr[l] <- exp(b); hl[l] <- exp(b - 1.96 * se); hh[l] <- exp(b + 1.96 * se) } }
  list(risk = risk, hr = hr, hr_lo = hl, hr_hi = hh)
}

# y = KSAA quartile, x = comparison measure (oriented so that right is worse); reference cell is Q1 x first level, lowest on both
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

dM <- d[!is.na(d$pa_accel), ];  dM$xq <- qcut(dM$pa_accel)          # main panel: PhenoAge acceleration (higher is worse)
dB <- d[!is.na(d$tot_bmd), ];   dB$xq <- qcut(-dB$tot_bmd)          # supplementary: total hip BMD (negated, Q4 = lowest BMD)
KL <- c("0–1", "≥2")                                       # supplementary: K–L in two levels (phase 2)
dC <- d[!is.na(d$klmax), ];     dC$xq <- factor(ifelse(dC$klmax <= 1, KL[1], KL[2]), KL)
M <- panel_data(dM, QL, "Q1", "main"); B <- panel_data(dB, QL, "Q1", "bmd"); C <- panel_data(dC, KL, KL[1], "kl")
MV <- panel_data(dM, QL, "Q1", "main_cvd", ev = "e_cvd")
MK <- panel_data(dM, QL, "Q1", "main_cancer", ev = "e_cancer")  # cancer death version (alternative)

G <- rbind(M$grid, B$grid, C$grid, MV$grid, MK$grid); W <- rbind(M$within, B$within, C$within, MV$within, MK$within)
write.csv(G, file.path(V3, "joint_cells.csv"), row.names = FALSE)
write.csv(W, file.path(V3, "joint_within.csv"), row.names = FALSE)

ov <- function(dat) 100 * mean(dat$xq[qcut(dat$gap) == "Q4"] == "Q4", na.rm = TRUE)
cat(sprintf("\nOverlap (percentage of the top KSAA quartile also in the worst quartile of the comparison measure; independence = 25%%)\n  PhenoAge acceleration %.1f%%\n  Total hip BMD %.1f%%\n", ov(dM), ov(dB)))
cat("\nHR per SD of KSAA within each level (age and sex adjusted)\n")
for (i in seq_len(nrow(W))) with(W[i, ], cat(sprintf("  %-5s %-6s n=%-5d deaths %-5d HR %s\n", panel, stratum, n, deaths, ci(hr, lo, hi))))
cat(sprintf("\nSamples: main n=%d  bmd n=%d  kl n=%d\n", M$n, B$n, C$n))
