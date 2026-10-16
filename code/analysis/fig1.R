# Figure 1:研究流程(A)、交叉拟合与年龄分布(B)、年龄预测性能(C)、全因死亡(D)、与 PhenoAge 加速的联合分层(E)、全表型扫描(F)。
# 输入 results/analysis/*(死亡相关为调查加权结果 weighted/*,先跑 weighted.R / weighted_boot.R)、results/predictions/m4_<骨干>_oof.csv、code/analysis/phewas_labels.R、示例膝片 ROI
# 输出 results/figures/Figure1/Figure1.{png,pdf,tiff}(183 × 166.5 mm)及 panels/ 下的子图 PDF;在项目根目录下运行
suppressPackageStartupMessages({library(grid); library(ggplot2); library(survival); library(ggrepel)})
source("code/analysis/phewas_labels.R")
out <- "results/figures/Figure1"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- readRDS("results/analysis/analysis_data.rds")
stopifnot(nrow(d) == 4875, !anyDuplicated(d$SEQN), sum(d$e) == 4399)

# ---------------- 数字(全部由数据计算并与稿件核对) ----------------
PF <- do.call(rbind, lapply(split(d, d$fold), function(z) data.frame(MAE = mean(abs(z$pred_raw - z$age)), r = cor(z$age, z$pred_raw))))
stopifnot(round(mean(PF$MAE), 2) == 4.37, round(mean(PF$r), 3) == 0.736)
# 五骨干的折外原始预测(与集成 pred_raw 同口径:未做年龄偏倚校正)
MODS <- c(ensemble = "Weighted ensemble", resnet18 = "ResNet-18", resnet50 = "ResNet-50", densenet121 = "DenseNet-121",
          efficientnetb0 = "EfficientNet-B0", convnexttiny = "ConvNeXt-Tiny")
PR <- sapply(names(MODS)[-1], function(a) { x <- read.csv(sprintf("results/predictions/m4_%s_oof.csv", a)); x$pred_age[match(d$SEQN, x$SEQN)] })
PR <- cbind(ensemble = d$pred_raw, PR); stopifnot(!anyNA(PR))
CS5 <- sapply(split(seq_len(nrow(d)), d$fold), function(i) 100 * mean(abs(PR[i, "ensemble"] - d$age[i]) <= 5))
stopifnot(round(mean(CS5), 1) == 64.5, round(sd(CS5), 1) == 1.7)  # CS (α = 5 y) 64.5 ± 1.7
# 主分析为调查加权
M1 <- read.csv("results/analysis/weighted/model1_full_cohort.csv"); M1 <- M1[M1$outcome == "allcause", ]
CI7 <- read.csv("results/analysis/phenoage_combined_cindex.csv"); stopifnot(all(CI7$n == 4522), round(CI7$dC_vs_pa[1], 3) == 0.004)
PW <- read.csv("results/analysis/phewas/phewas_results.csv"); PW <- PW[PW$set == "main" & !is.na(PW$p), ]
stopifnot(nrow(PW) == 482, sum(PW$fdr < 0.05) == 109)
f2 <- function(x, k = 2) sub("^-", "−", formatC(x, format = "f", digits = k))
cm <- function(x) formatC(x, format = "d", big.mark = ",")

# ---------------- 画布与基本元素(坐标单位 mm,y 自上而下) ----------------
W <- 183; H <- 166.5
C <- list(ink = "#1F2D3D", muted = "#65717E", rule = "#CFDCE8", hdr = "#DCEAF6", body = "#FFFFFF", blue = "#2F6FAE", light = "#9EC3E6",
          teal = "#23897F", gold = "#C07A34", grey = "#98A2AC", red = "#CF4A3C",
          mort = c("#FBE4E1", "#B8372C"), frac = c("#DDF0E8", "#1F7A5E"), func = c("#E0EBF7", "#2A5E9A"), mark = c("#ECE6F5", "#674A92"))
QCOL <- c(Q1 = "#9DB4C8", Q2 = "#4E8FD0", Q3 = "#1B3F7A", Q4 = "#C8372D")  # 与 Figure 2 相同
tx <- function(x, y, label, size = 5.6, col = C$ink, bold = FALSE, just = "left", rot = 0, italic = FALSE)
  grid.text(label, x = unit(x, "mm"), y = unit(H - y, "mm"), just = just, rot = rot,
            gp = gpar(fontfamily = "Arial", fontsize = size, col = col, lineheight = 1.05,
                      fontface = if (bold && italic) "bold.italic" else if (bold) "bold" else if (italic) "italic" else "plain"))
rbox <- function(x, y, w, h, fill = "white", border = C$rule, r = 1.2, lwd = .5)
  grid.roundrect(x = unit(x, "mm"), y = unit(H - y, "mm"), width = unit(w, "mm"), height = unit(h, "mm"),
                 just = c("left", "top"), r = unit(r, "mm"), gp = gpar(fill = fill, col = border, lwd = lwd))
rect <- function(x, y, w, h, fill, col = NA, lwd = .4)
  grid.rect(x = unit(x, "mm"), y = unit(H - y, "mm"), width = unit(w, "mm"), height = unit(h, "mm"), just = c("left", "top"), gp = gpar(fill = fill, col = col, lwd = lwd))
seg <- function(x, y, x2, y2, col = C$muted, lwd = .6, lty = 1, arrow = FALSE, head = 1.2)
  grid.lines(x = unit(c(x, x2), "mm"), y = unit(H - c(y, y2), "mm"), gp = gpar(col = col, fill = col, lwd = lwd, lty = lty),
             arrow = if (arrow) grid::arrow(length = unit(head, "mm"), type = "closed") else NULL)
dot <- function(x, y, r = .6, col = C$blue, fill = col, lwd = .6)
  grid.circle(unit(x, "mm"), unit(H - y, "mm"), r = unit(r, "mm"), gp = gpar(fill = fill, col = col, lwd = lwd))
path <- function(x, y, col = C$blue, lwd = 1, fill = NA) {
  if (is.na(fill)) grid.lines(unit(x, "mm"), unit(H - y, "mm"), gp = gpar(col = col, lwd = lwd, lineend = "round", linejoin = "round"))
  else grid.polygon(unit(x, "mm"), unit(H - y, "mm"), gp = gpar(col = col, fill = fill, lwd = lwd))
}
bigarrow <- function(x, y, col = "#7FA9D1", s = 1) grid.polygon(unit(x + c(0, 2.2, 2.2, 4, 2.2, 2.2, 0) * s, "mm"),
  unit(H - (y + c(-.9, -.9, -2, 0, 2, .9, .9) * s), "mm"), gp = gpar(fill = col, col = NA))
card <- function(x, y, w, h, title, hdr = C$hdr, tcol = C$ink, letter = NULL, hh = 6, tsize = 6.6, lsize = 7.2, strip = TRUE) {
  rbox(x, y, w, h, "white", C$rule)
  if (strip) grid.roundrect(x = unit(x, "mm"), y = unit(H - y, "mm"), width = unit(w, "mm"), height = unit(hh, "mm"), just = c("left", "top"),
                 r = unit(1.2, "mm"), gp = gpar(fill = hdr, col = NA))
  if (strip) rect(x, y + hh - 1.3, w, 1.3, hdr)  # 标题栏下沿取直
  if (is.null(letter)) tx(x + w / 2, y + hh / 2, title, tsize, tcol, TRUE, "centre")
  else { tx(x + 2, y + hh / 2, letter, 9, C$ink, TRUE); tx(x + 7.5, y + hh / 2, title, lsize, C$ink, TRUE) }
}
person <- function(x, y, col, s = 1) {
  dot(x, y, .8 * s, col); path(x + c(-1.1, -1.1, 1.1, 1.1) * s, y + c(4.6, 1.3, 1.3, 4.6) * s, col, .4, fill = col)
  seg(x - .55 * s, y + 4.4 * s, x - .55 * s, y + 7 * s, col, 1.6 * s); seg(x + .55 * s, y + 4.4 * s, x + .55 * s, y + 7 * s, col, 1.6 * s)
}
features <- function(x, y, col = C$blue, s = 1) {
  for (i in 2:0) rbox(x + i * .7 * s, y - i * .6 * s, 4 * s, 5.2 * s, "#DCE8F3", col, .2, .45)
  for (i in 0:2) for (j in 0:3) dot(x + (.8 + i * 1.1) * s, y + (.8 + j * 1.15) * s, .17 * s, col)
}
subcard <- function(x, y, w, h, title, hh = 5, hdr = "#E8F1FA") {  # 组卡内的子卡(A 面板)
  rbox(x, y, w, h, "white", C$rule, 1, .45)
  grid.roundrect(x = unit(x, "mm"), y = unit(H - y, "mm"), width = unit(w, "mm"), height = unit(hh, "mm"), just = c("left", "top"), r = unit(1, "mm"), gp = gpar(fill = hdr, col = NA))
  rect(x, y + hh - 1.1, w, 1.1, hdr); tx(x + w / 2, y + hh / 2, title, 6, C$ink, TRUE, "centre")
}
tick <- function(x, y, col = "white") path(x + c(-1, -.2, 1.2), y + c(0, .8, -.9), col, .9)
plot_at <- function(p, x, y, w, h) print(p, vp = viewport(x = unit(x, "mm"), y = unit(H - y, "mm"), width = unit(w, "mm"), height = unit(h, "mm"), just = c("left", "top")))
# 上下两张 ggplot 纵轴线对齐:绘图区两侧各列宽取两图较大者(须在已打开的设备里调用)
align_y <- function(p1, p2) { g1 <- ggplotGrob(p1); g2 <- ggplotGrob(p2); stopifnot(ncol(g1) == ncol(g2)); pl <- g1$layout[g1$layout$name == "panel", ]
  side <- setdiff(seq_len(ncol(g1)), pl$l:pl$r); mw <- grid::unit.pmax(g1$widths[side], g2$widths[side]); g1$widths[side] <- mw; g2$widths[side] <- mw; list(g1, g2) }
grob_at <- function(g, x, y, w, h) { pushViewport(viewport(x = unit(x, "mm"), y = unit(H - y, "mm"), width = unit(w, "mm"), height = unit(h, "mm"), just = c("left", "top"))); grid.draw(g); popViewport() }
# 画 ggplot 并使绘图区下沿(横轴线)恰落在 yaxis;返回视口高(mm)
plot_axis_at <- function(p, x, ytop, w, yaxis) { g <- ggplotGrob(p); pr <- g$layout[g$layout$name == "panel", ]
  hB <- sum(convertHeight(g$heights[(pr$b + 1):length(g$heights)], "mm", valueOnly = TRUE)); h <- yaxis - ytop + hB
  pushViewport(viewport(x = unit(x, "mm"), y = unit(H - ytop, "mm"), width = unit(w, "mm"), height = unit(h, "mm"), just = c("left", "top"))); grid.draw(g); popViewport(); invisible(h) }
th <- theme_classic(base_size = 6, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .22, colour = C$muted), axis.ticks = element_line(linewidth = .2, colour = C$muted),
        axis.ticks.length = unit(.6, "mm"), axis.text = element_text(colour = C$muted, size = 5.4), axis.title = element_text(size = 6, colour = C$ink),
        plot.margin = margin(1, 2, 1, 1), plot.background = element_blank(), panel.background = element_blank(), legend.position = "none")

# ---------------- 子图 ----------------
# B 右:实足年龄分布(年龄在纵轴,5 岁直方条 + 核密度曲线)。NHANES III 90 岁顶码:≥90 岁单列一条,画在 90–95 并标 "90+"。
# 核密度:点位 a + 0.5,带宽 2.5 岁,在入组下限 60 岁处反射校正(否则曲线在 60 岁处虚假收尖);乘 5 换算成每 5 岁人数,与直方条同刻度
AG <- data.frame(lo = seq(60, 90, 5)); AG$Freq <- sapply(AG$lo, function(l) sum(d$age >= l & (d$age < l + 5 | l == 90)))
stopifnot(sum(AG$Freq) == nrow(d), range(d$age) == c(60, 90), AG$Freq[7] == sum(d$age >= 90), AG$Freq[1] == sum(d$age <= 64))
XM <- ceiling(max(AG$Freq) / 200) * 200; stopifnot(round(mean(d$age), 1) == 71.8, round(sd(d$age), 1) == 8.1)
KD <- data.frame(age = seq(60, 95, .1)); KD$n5 <- 5 * sapply(KD$age, function(x) sum(dnorm(x, d$age + .5, 2.5) + dnorm(x, 120 - (d$age + .5), 2.5)))
AGE_W <- 35.5; AGE_H <- NA  # 子图总宽(mm);总高在 draw() 里按与 Pooled OOF 列对齐算出
E_H <- NA   # E 子图总高(mm),在 draw() 里按横轴对齐算出
AGE_COL <- "#1A6B4A"  # 直方条半透明,免得遮住曲线
p_age <- ggplot(AG) + geom_rect(aes(xmin = 0, xmax = Freq, ymin = lo + .12, ymax = lo + 4.88), fill = AGE_COL, alpha = .5) +
  geom_path(data = KD, aes(n5, age), colour = AGE_COL, linewidth = .55, lineend = "round") +
  annotate("text", x = XM, y = 93.3, label = sprintf("n = %s\nMean (SD), %.1f (%.1f) years", cm(nrow(d)), mean(d$age), sd(d$age)), hjust = 1, vjust = 0,
           size = 4.4 / .pt, family = "Arial", colour = C$muted, lineheight = 1) +
  scale_x_continuous(breaks = seq(0, XM, 500), labels = cm, limits = c(0, XM), expand = expansion(mult = c(0, .02))) +
  scale_y_reverse(breaks = seq(60, 90, 10), labels = c(60, 70, 80, "90+"), limits = c(95, 60), expand = c(0, 0)) +  # 绘图区恰为 60–95
  labs(x = "Number of participants", y = "Chronological age\n(years)") + th + theme(axis.title.y = element_text(lineheight = .95))  # 纵轴标题分两行
# C 上:散点 + 拟合线;r 与 MAE 为五个留出折的均值 ± SD
fit <- coef(lm(pred_raw ~ age, d))
SHOW_IDENTITY <- FALSE  # TRUE 时加画恒等线
p_sc <- ggplot(d, aes(age, pred_raw)) + geom_point(size = .18, alpha = .16, colour = C$blue, shape = 16) +
  (if (SHOW_IDENTITY) geom_abline(slope = 1, intercept = 0, colour = "#8D98A3", linetype = 2, linewidth = .3)) +
  geom_abline(intercept = fit[1], slope = fit[2], colour = C$teal, linewidth = .5) +
  annotate("text", x = 60.4, y = 95.2, hjust = 0, vjust = 1, size = 5.4 / .pt, family = "Arial", colour = C$ink, lineheight = 1.1,
           label = sprintf("Pearson r = %s ± %s\nMAE = %s ± %s years", f2(mean(PF$r), 3), f2(sd(PF$r), 3), f2(mean(PF$MAE)), f2(sd(PF$MAE)))) +
  annotate("segment", x = 79.6, xend = 81.4, y = c(54.5, 50.4)[seq_len(1 + SHOW_IDENTITY)], yend = c(54.5, 50.4)[seq_len(1 + SHOW_IDENTITY)],
           colour = c(C$teal, "#8D98A3")[seq_len(1 + SHOW_IDENTITY)], linetype = c(1, 2)[seq_len(1 + SHOW_IDENTITY)], linewidth = c(.5, .3)[seq_len(1 + SHOW_IDENTITY)]) +
  annotate("text", x = 81.9, y = c(54.5, 50.4)[seq_len(1 + SHOW_IDENTITY)], label = c("Pooled linear fit", "Identity line")[seq_len(1 + SHOW_IDENTITY)], hjust = 0,
           size = 5 / .pt, family = "Arial", colour = C$ink) +
  scale_x_continuous(breaks = seq(60, 90, 10), labels = c(60, 70, 80, "90+")) + scale_y_continuous(breaks = seq(50, 90, 10)) +
  coord_cartesian(xlim = c(59.5, 90.5), ylim = c(48, 96), expand = FALSE) +
  labs(x = "Chronological age (years)", y = "Predicted knee age\n(years)") + th + theme(axis.title.y = element_text(lineheight = .95))
# C 下:累积得分 CS(x) = 绝对误差 ≤ x 岁者的百分比(合并折外预测);x = 5 处标集成的 CS5(折均值 ± SD)
XS <- seq(0, 16, .1)
CSD <- do.call(rbind, lapply(names(MODS), function(m) { e <- abs(PR[, m] - d$age); data.frame(model = MODS[[m]], x = XS, cs = 100 * sapply(XS, function(x) mean(e <= x))) }))
CSD$model <- factor(CSD$model, unname(MODS)); CSD <- CSD[order(CSD$model != MODS[["ensemble"]], decreasing = TRUE), ]
MCOL <- setNames(c(C$ink, "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"), MODS)
MLTY <- setNames(c("solid", "42", "22", "42", "22", "42"), MODS); MLWD <- setNames(c(.6, rep(.32, 5)), MODS)
cs5 <- 100 * mean(abs(PR[, "ensemble"] - d$age) <= 5)
# 图例顺序同方法(集成在前)
BRK <- unname(MODS[c("ensemble", "resnet18", "resnet50", "efficientnetb0", "convnexttiny", "densenet121")])
p_cs <- ggplot(CSD, aes(x, cs, colour = model, linetype = model, linewidth = model)) +
  annotate("segment", x = 5, xend = 5, y = 0, yend = cs5, colour = C$grey, linetype = "11", linewidth = .3) +
  geom_line(lineend = "round") +
  annotate("point", x = 5, y = cs5, size = .7, colour = C$ink) +
  annotate("text", x = 5.5, y = cs5 - 5, hjust = 0, vjust = 1, size = 5.2 / .pt, family = "Arial", colour = C$ink, lineheight = 1.05,
           label = sprintf("CS (α = 5 years)\n%s ± %s%%", f2(mean(CS5), 1), f2(sd(CS5), 1))) +
  scale_colour_manual(values = MCOL, breaks = BRK, name = NULL) + scale_linetype_manual(values = MLTY, breaks = BRK, name = NULL) +
  scale_linewidth_manual(values = MLWD, breaks = BRK, name = NULL) +
  scale_x_continuous(breaks = seq(0, 16, 4), limits = c(0, 16), expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100), expand = expansion(mult = c(0, .03))) +
  labs(x = "Error level x (years)", y = "Cumulative score (%)") + th +
  theme(legend.position = c(.99, .03), legend.justification = c(1, 0), legend.background = element_blank(), legend.key = element_blank(),
        legend.key.height = unit(2.05, "mm"), legend.key.width = unit(4.6, "mm"), legend.text = element_text(size = 5, colour = C$ink),
        legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0))
# D:与 Figure 2 同一组年龄、性别标化曲线(weighted.R 写出)、同一四分位配色
STD <- read.csv("results/analysis/weighted/figure2_std_curves.csv"); STD$time <- STD$month / 12; STD$q <- factor(STD$q, names(QCOL))
# 10 年标化死亡(Q1 / Q4),取同一组曲线的 t = 10 年
A10 <- STD[STD$month == 120, ]; A10 <- setNames(A10$allcause, as.character(A10$q))
stopifnot(round(100 * A10[c("Q1", "Q4")], 1) == c(31.7, 41.3))
p_km <- ggplot(STD, aes(time, allcause, colour = q)) + geom_step(linewidth = .45) + scale_colour_manual(values = QCOL, name = "KSAA quartile") +
  scale_x_continuous(breaks = seq(0, 30, 10), limits = c(0, 30), expand = expansion(mult = c(0, .02))) +
  scale_y_continuous(breaks = seq(0, 1, .25), limits = c(0, 1), labels = function(v) sprintf("%.0f%%", 100 * v), expand = c(0, 0)) +
  annotate("segment", x = 10, xend = 10, y = 0, yend = A10[["Q4"]], linetype = "11", linewidth = .3, colour = C$muted) +
  annotate("point", x = 10, y = A10[c("Q1", "Q4")], colour = QCOL[c("Q1", "Q4")], size = .9) +
  # 数值标在两点旁,避免压到曲线
  annotate("text", x = c(9.5, 10.5), y = c(A10[["Q4"]] + .015, A10[["Q1"]] - .015), hjust = c(1, 0), vjust = c(0, 1),
           label = sprintf("%.1f%%", 100 * A10[c("Q4", "Q1")]), size = 5.2 / .pt, family = "Arial", colour = c(QCOL[["Q4"]], "#6E8AA3")) +  # Q1 浅色作文字对比度不够,取同色相加深
  annotate("text", x = .6, y = .97, hjust = 0, vjust = 1, size = 5.4 / .pt, family = "Arial", colour = C$ink, lineheight = 1.05,
           label = sprintf("Age- and sex-adjusted\nHR per SD, %s (95%% CI %s–%s)", f2(M1$hr_sd), f2(M1$lo_sd), f2(M1$hi_sd))) +
  labs(x = "Years since examination", y = "Standardized cumulative mortality") + th + guides(colour = guide_legend(ncol = 2)) +
  theme(legend.position = c(.98, .03), legend.justification = c(1, 0), legend.background = element_blank(), legend.key.size = unit(2, "mm"),
        legend.title = element_text(size = 5.4), legend.text = element_text(size = 5.2), legend.key.height = unit(2.1, "mm"))
# E:全因死亡按 PhenoAge 加速 × KSAA 四分位交叉分格的 HR(年龄、性别调整;参照格 Q1×Q1);配色同 Figure 3
c3 <- read.csv("results/analysis/weighted/figure3_cells.csv"); c3 <- c3[c3$outcome == "allcause", ]
JC1 <- data.frame(panel = "main", y = c3$y, x = c3$x, hr = c3$hr, ref = c3$y == "Q1" & c3$x == "Q1", n = c3$n)
QL4 <- paste0("Q", 1:4); JC1$x <- factor(JC1$x, QL4); JC1$y <- factor(JC1$y, QL4); JC1$dark <- abs(log(JC1$hr)) > 0.8
stopifnot(f2(JC1$hr[JC1$x == "Q4" & JC1$y == "Q4"]) == "3.93", JC1$ref[JC1$x == "Q1" & JC1$y == "Q1"])  # 参照格 Q1×Q1
p_joint <- ggplot(JC1, aes(x, y, fill = log(hr))) + geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = ifelse(ref, "Ref", f2(hr)), colour = dark), size = 5.2 / .pt, family = "Arial", fontface = "bold") +
  scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#eb6834", midpoint = 0, limits = c(-log(5), log(5))) +
  scale_colour_manual(values = c(`FALSE` = C$ink, `TRUE` = "white"), guide = "none") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = QL4, expand = c(0, 0)) +
  labs(x = "Phenotypic Age acceleration quartile", y = "KSAA quartile", subtitle = "All-cause mortality, HR vs Q1 \u00d7 Q1") + th +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(colour = C$ink, size = 5.4),
        plot.subtitle = element_text(size = 5.4, colour = C$muted, hjust = .5, margin = margin(b = 1.5)))
# F 曼哈顿:与 Figure 4 同一批 482 个表型、同一类别配色;类别按字母序,灰白交替底;红点线 FDR 5%
CORD <- sort(unname(CLS_NAME[unique(PW$health_class)]))
PW$cat <- factor(CLS_NAME[PW$health_class], CORD); PW <- PW[order(PW$cat, PW$variable), ]; PW$x <- NA_real_; cur <- 0
for (cc in CORD) { i <- which(PW$cat == cc); PW$x[i] <- cur + seq_along(i); cur <- cur + length(i) + 6 }
BAND <- do.call(rbind, lapply(seq_along(CORD), function(k) { z <- PW$x[PW$cat == CORD[k]]
  data.frame(xmin = min(z) - 3, xmax = max(z) + 3, mid = mean(range(z)), cat = CORD[k], shade = k %% 2 == 1) }))
PW$y <- -log10(PW$p); PW$dir <- factor(ifelse((grepl("^beta", PW$measure) & PW$est > 0) | (!grepl("^beta", PW$measure) & PW$est > 1), "pos", "neg"), c("pos", "neg")); lead <- PW[!is.na(PW$clump_lead) & PW$clump_lead == PW$variable, ]; lead <- head(lead[order(lead$p), ], 3)
LBL <- rbind(data.frame(x = lead$x, y = lead$y, lab = c("Trochanter BMD \u2193", "Abnormal locomotion \u2191", "Walks slowly \u2191")),
             data.frame(x = PW$x[PW$variable == "CRP"], y = PW$y[PW$variable == "CRP"], lab = "CRP"))
stopifnot(identical(lead$variable, c("BDPTRBMD", "PEP1", "HAY11D")), identical(as.character(lead$dir), c("neg", "pos", "pos")))
LABEL_CRP <- FALSE  # TRUE 时标注参照项 CRP
p_man <- ggplot(PW, aes(x, y)) +
  geom_rect(data = BAND[BAND$shade, ], aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#EDEDED") +
  geom_hline(yintercept = -log10(max(PW$p[PW$fdr < 0.05])), colour = "red", linetype = "11", linewidth = .35) +
  geom_point(aes(colour = health_class, shape = dir), size = .42) + scale_shape_manual(values = c(pos = 16, neg = 17)) +
  geom_text_repel(data = LBL[LBL$lab != "CRP", ], aes(label = lab), size = 5 / .pt, family = "Arial", colour = C$ink, segment.size = .2, segment.colour = C$muted,
                  min.segment.length = 0, box.padding = .25, nudge_x = 40, direction = "y", hjust = 0, seed = 1) +
  (if (LABEL_CRP) geom_text_repel(data = LBL[LBL$lab == "CRP", ], aes(label = lab), size = 5 / .pt, family = "Arial", colour = C$ink, segment.size = .2, segment.colour = C$muted,
                  min.segment.length = 0, nudge_y = 5.5, nudge_x = 15, seed = 1)) +
  annotate("point", x = min(BAND$xmin) + 12, y = c(17.6, 15.9), shape = c(16, 17), size = .8, colour = C$muted) +
  annotate("text", x = min(BAND$xmin) + 23, y = c(17.6, 15.9), label = c("Positive", "Negative"), hjust = 0, size = 4.8 / .pt, family = "Arial", colour = C$muted) +
  scale_colour_manual(values = CLS_COL) +
  scale_x_continuous(breaks = BAND$mid, labels = NULL, expand = c(0, 0), limits = c(min(BAND$xmin), max(BAND$xmax))) +
  scale_y_continuous(limits = c(0, 19), expand = c(0, 0)) +
  labs(x = "Phenotype", y = expression(-log[10]*italic(P))) + th +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ---------------- 组版 ----------------
draw <- function() {
  grid.newpage()
  # ===== A:人群 + 膝片 → Knee-Age model → KSAA =====
  Y0 <- 1.5; HA <- 49; YM <- 31.5
  card(1.5, Y0, 63, HA, "Study population and knee images", letter = "A", hh = 6.2)
  subcard(3.3, Y0 + 8, 28.5, HA - 9.7, "Population"); subcard(34.2, Y0 + 8, 28.5, HA - 9.7, "Knee images")
  tx(17.55, 17.7, "NHANES III", 6.2, bold = TRUE, just = "centre"); tx(17.55, 20.8, "1988–1994", 6, bold = TRUE, just = "centre")
  tx(17.55, 23.9, bquote(bolditalic(n) ~ bold(.(paste("=", cm(nrow(d)))))), 6, just = "centre"); tx(17.55, 27, "Aged ≥60 years", 5.6, just = "centre")
  for (r in 1:2) { yy <- c(29.4, 36.6)[r]; cc <- c("#A9C6E2", "#3F74A9")[r]
    for (xx in c(8.2, 12.6, 22.5, 26.9)) person(xx, yy, cc, .75); for (xx in c(16.4, 17.55, 18.7)) dot(xx, yy + 2.6, .28, cc) }
  tx(17.55, 44.6, "Community-dwelling", 4.8, C$muted, just = "centre"); tx(17.55, 46.8, "US adults", 4.8, C$muted, just = "centre")
  # 左膝按原始方向显示(镜像仅为模型输入的预处理)
  img <- lapply(c("L", "R"), function(side) { z <- png::readPNG(paste0("data/images/example_roi/10018K1_", side, ".png"))
    q <- quantile(z, c(.005, .995)); z[] <- pmax(0, pmin(1, (z - q[1]) / (q[2] - q[1]))); z })
  for (i in 1:2) grid.raster(img[[i]], x = unit(35.25 + (i - 1) * 13.5, "mm"), y = unit(H - 22.8, "mm"), width = unit(12.9, "mm"), height = unit(12.9, "mm"), just = c("left", "top"))
  tx(41.7, 38.3, "Left knee", 5.2, just = "centre"); tx(55.2, 38.3, "Right knee", 5.2, just = "centre")
  card(66, Y0, 59.5, HA, "Knee-Age model", hh = 6.2, tsize = 7.2)
  tx(95.75, Y0 + 9.6, "Multi-view deep learning model", 5.8, just = "centre"); tx(95.75, Y0 + 12.6, "(5 backbones, ensemble)", 5.6, just = "centre")
  bn <- c("ResNet-18", "ResNet-50", "EfficientNet-B0", "ConvNeXt-Tiny", "DenseNet-121")  # 顺序同方法
  for (i in 1:5) { rbox(68, 18.1 + (i - 1) * 5.1, 16.5, 4, "#E3EDF7", NA, .6); tx(76.25, 20.1 + (i - 1) * 5.1, bn[i], 5, just = "centre") }
  path(c(85.2, 86.1, 86.1, 85.2), c(18.4, 18.4, 42.2, 42.2), C$muted, .5); seg(86.1, 30.3, 88.3, 30.3, C$muted, .5, arrow = TRUE, head = .9)
  features(90.6, 22.2, s = 1.15); features(90.6, 33.6, s = 1.15); tx(89.9, 25.2, "L", 5.2, bold = TRUE, just = "right"); tx(89.9, 36.6, "R", 5.2, bold = TRUE, just = "right")
  seg(97.6, 25.2, 100.4, 28.4, C$muted, .5, arrow = TRUE, head = .9); seg(97.6, 36.4, 100.4, 32.4, C$muted, .5, arrow = TRUE, head = .9)
  for (i in 0:2) for (j in 0:2) seg(102.2, 25.3 + i * 5, 109.6, 25.3 + j * 5, "#9FCBC5", .4)
  for (i in 0:2) { dot(102.2, 25.3 + i * 5, .95, C$teal, "white"); dot(109.6, 25.3 + i * 5, .95, C$teal, "#CDE7E3") }
  tx(105.9, 20.8, "Transformer", 5.4, just = "centre")
  seg(111.6, 30.3, 114.2, 30.3, C$muted, .5, arrow = TRUE, head = .9)
  for (i in 1:9) { h <- 7.6 * exp(-((i - 5) / 2.1)^2); rect(115 + (i - 1) * 1.02, 33.6 - h, .78, h, C$gold) }
  seg(114.6, 33.6, 124.3, 33.6, C$muted, .5)
  tx(119.5, 37.6, "Mean of two\nknee ages", 4.8, C$muted, just = "centre")
  card(127, Y0, 54.5, HA, "KSAA", hdr = C$frac[1], tcol = C$teal, hh = 6.2, tsize = 7.6)
  tx(154.25, Y0 + 10, "Age-bias-corrected residual", 5.8, C$teal, just = "centre")
  seg(137, 42.5, 137, 16.5, C$ink, .6); seg(137, 42.5, 178, 42.5, C$ink, .6)
  tx(132.2, 29.5, "Predicted knee age\n(years)", 5, just = "centre", rot = 90); tx(157.5, 45.4, "Chronological age (years)", 5, just = "centre")
  ly <- function(x) 38.7 - (x - 138.5) * 16.5 / 32
  seg(138.5, ly(138.5), 170.5, ly(170.5), "#8D98A3", .7, 2)
  seg(161, ly(161) - .2, 161, 20.9, C$teal, .8, arrow = TRUE, head = 1); dot(161, 19.7, 1, C$teal)
  tx(163.6, 17.6, "+ KSAA", 5.6, C$teal, TRUE); tx(163.6, 19.9, "(older-appearing)", 4.8, C$teal)
  seg(164.2, ly(164.2) + .3, 166.1, 26.3, C$muted, .3)
  tx(166.6, 26.6, "Expected", 4.5, C$muted); tx(166.6, 28.7, "predicted knee age", 4.5, C$muted); tx(166.6, 30.8, "(α + β × age)", 4.5, C$muted)
  seg(148.5, ly(148.5) + .2, 148.5, 38.1, C$gold, .8, arrow = TRUE, head = 1); dot(148.5, 39.3, 1, C$gold)
  tx(151, 37.4, "− KSAA", 5.6, C$gold, TRUE); tx(151, 39.7, "(younger-appearing)", 4.8, C$gold)
  for (xx in c(31.3, 63.55, 124.55)) bigarrow(xx, YM, s = .85)  # 大箭头压在卡片边线上,须在所有卡片之后画
  # ===== B 交叉拟合 + 年龄分布 / C 年龄预测(等高) =====
  Y1 <- 52; HB <- 63.5
  card(1.5, Y1, 98.5, HB, "Cross-fitting", letter = "B", hh = 6.2)
  # 年龄分布绘图区上、下沿与汇总列最上、最下一格对齐
  tx(4.5, Y1 + 9.4, "Outer five-fold cross-fitting", 6, bold = TRUE)
  GX <- 11.5; GP <- 7.6; GW <- 6.9; GY <- Y1 + 16.2; RP <- 4.3; RH <- 3.5; GE <- GY + 4 * RP + RH; GM <- (GY + GE) / 2; PX <- 55.6; PW_ <- 5.2
  tx(4.8, GY - 2.3, "Round", 5); for (j in 1:5) tx(GX + (j - 1) * GP + GW / 2, GY - 2.3, paste("Fold", j), 5, just = "centre")
  for (i in 1:5) { yc <- GY + (i - 1) * RP + RH / 2; tx(6.9, yc, as.character(i), 5.2, just = "centre")
    for (j in 1:5) { rbox(GX + (j - 1) * GP, GY + (i - 1) * RP, GW, RH, if (i == j) C$teal else "#E3EBF2", NA, .5)
      if (i == j) tick(GX + (j - 1) * GP + GW / 2 - .1, yc + .05) }
    rbox(PX, GY + (i - 1) * RP, PW_, RH, C$teal, NA, .5) }
  path(c(49.7, 50.5, 50.5, 49.7), c(GY - .3, GY - .3, GE + .3, GE + .3), C$muted, .6); bigarrow(51.4, GM, s = .8)
  tx(PX + PW_ / 2, GY - 3.2, "Pooled OOF\npredictions", 4.8, C$teal, TRUE, just = "centre"); tx(PX + PW_ / 2, GE + 2, sprintf("(n = %s)", cm(nrow(d))), 4.8, just = "centre")
  tx(64, Y1 + 9.4, "Chronological age distribution", 6, bold = TRUE)
  gA <- ggplotGrob(p_age); pr <- gA$layout[gA$layout$name == "panel", ]  # 量出绘图区上下方高度,使绘图区恰落在 GY … GE
  hA <- sum(convertHeight(gA$heights[seq_len(pr$t - 1)], "mm", valueOnly = TRUE)); hB <- sum(convertHeight(gA$heights[(pr$b + 1):length(gA$heights)], "mm", valueOnly = TRUE))
  AGE_H <<- hA + (GE - GY) + hB
  pushViewport(viewport(x = unit(63.5, "mm"), y = unit(H - (GY - hA), "mm"), width = unit(AGE_W, "mm"), height = unit(AGE_H, "mm"), just = c("left", "top"))); grid.draw(gA); popViewport()
  LY <- GE + 4.4
  rbox(5, LY, 4.6, 3.4, "#E3EBF2", NA, .4); tx(10.8, LY + 1.7, "Development (80%)", 5.2)
  rbox(29.5, LY, 4.6, 3.4, C$teal, NA, .4); tx(35.3, LY + 1.7, "Held-out prediction (20%)", 5.2)
  grid.roundrect(x = unit(4.5, "mm"), y = unit(H - (LY + 5.6), "mm"), width = unit(92.5, "mm"), height = unit(13.2, "mm"), just = c("left", "top"), r = unit(1.2, "mm"),
                 gp = gpar(fill = NA, col = "#AEBDCC", lwd = .6, lty = "22"))
  tx(6.8, LY + 8.6, "Within each development set", 6, bold = TRUE)
  rbox(7.3, LY + 11, 63.9, 4.6, "#C9DCEC", NA, .5); rbox(71.5, LY + 11, 7.1, 4.6, "#E6C39D", NA, .5)
  tx(39.25, LY + 13.3, "90% training", 5.4, just = "centre"); tx(80.2, LY + 13.3, "10% internal\nvalidation", 5.2, C$gold)
  card(101.5, Y1, 80, HB, "Performance of age prediction", letter = "C", hh = 6.2)
  gC <- align_y(p_sc, p_cs)  # 上下两图纵轴线对齐
  grob_at(gC[[1]], 102.5, Y1 + 7, 78, 27.8); grob_at(gC[[2]], 102.5, Y1 + 35, 78, 28)  # 图例在图内
  # ===== D / E / F =====
  Y2 <- Y1 + HB + 1.5; HE <- 48
  card(1.5, Y2, 59, HE, "All-cause mortality by KSAA quartile", C$mort[1], C$mort[2], letter = "D", hh = 6, lsize = 6.2)
  YAX <- Y2 + 39.4  # D、E、F 三图横轴线对齐
  plot_axis_at(p_km, 2.5, Y2 + 7, 57, YAX)
  card(62, Y2, 59, HE, "Joint stratification with Phenotypic Age", C$func[1], C$func[2], letter = "E", hh = 6, lsize = 6.2)
  E_H <<- plot_axis_at(p_joint, 63, Y2 + 6.8, 57, YAX)
  card(122.5, Y2, 59, HE, "Phenome-wide correlates of KSAA", C$frac[1], C$frac[2], letter = "F", hh = 6, lsize = 6.2)
  # 阈值线图例省略(红点线 = FDR 5%,在图注说明)
  plot_axis_at(p_man, 123.5, Y2 + 7.4, 57, YAX)
  KEY <- list(c("disease", 124.5, Y2 + 44.1), c("healthcare", 141.5, Y2 + 44.1), c("symptom", 157, Y2 + 44.1),
              c("objective", 124.5, Y2 + 46.5), c("function", 141.5, Y2 + 46.5), c("self_rated", 157, Y2 + 46.5))
  for (z in KEY) { dot(as.numeric(z[2]), as.numeric(z[3]), .6, CLS_COL[[z[1]]]); tx(as.numeric(z[2]) + 1.4, as.numeric(z[3]), CLS_NAME[[z[1]]], 4.5) }
}
for (ext in c("png", "tiff")) {
  f <- file.path(out, paste0("Figure1.", ext))
  switch(ext, png = ragg::agg_png(f, width = W, height = H, units = "mm", res = 350),
         tiff = ragg::agg_tiff(f, width = W, height = H, units = "mm", res = 600, compression = "lzw"))
  draw(); invisible(dev.off())
}
# PDF 最后写:cairo 加载失败时关掉残留设备,改用 quartz(同 fig_common.R 的 save_fig)
write_pdf <- function(f, w, h, fun) {
  if (file.exists(f)) file.remove(f); had <- file.exists("Rplots.pdf"); on.exit(if (!had && file.exists("Rplots.pdf")) file.remove("Rplots.pdf"))  # cairo 失败时会留下 Rplots.pdf
  ok <- tryCatch({ grDevices::cairo_pdf(f, width = w / 25.4, height = h / 25.4, family = "Arial"); fun(); invisible(dev.off()); file.exists(f) && file.size(f) > 0 },
                 error = function(e) { while (dev.cur() > 1) dev.off(); FALSE })
  if (!ok) { if (file.exists(f)) file.remove(f)
    grDevices::quartz(file = f, type = "pdf", width = w / 25.4, height = h / 25.4, family = "Arial"); fun(); invisible(dev.off()) }
  if (file.exists(f)) (if (ok) "cairo" else "quartz") else "未写出"
}
cat("PDF:", write_pdf(file.path(out, "Figure1.pdf"), W, H, draw), "\n")
# 子图另存单独的矢量 PDF,尺寸与组版内一致(mm)
dir.create(file.path(out, "panels"), showWarnings = FALSE)
PANELS <- list(b_age_distribution = list(p_age, AGE_W, AGE_H), c_scatter = list(function() grid.draw(align_y(p_sc, p_cs)[[1]]), 78, 27.8), c_cumulative_score = list(function() grid.draw(align_y(p_sc, p_cs)[[2]]), 78, 28), e_joint = list(p_joint, 57, E_H))
for (nm in names(PANELS)) cat(sprintf("  panels/%s.pdf: %s\n", nm, write_pdf(file.path(out, "panels", paste0(nm, ".pdf")), PANELS[[nm]][[2]], PANELS[[nm]][[3]], function() { o <- PANELS[[nm]][[1]]; if (is.function(o)) o() else print(o) })))
cat(sprintf("Figure 1 v3 -> %s; MAE %.2f ± %.2f, r %.3f ± %.3f, CS5 %.1f ± %.1f; HR/SD %.2f; dC all-cause %.4f; PheWAS %d phenotypes, %d at FDR<5%%\n",
            out, mean(PF$MAE), sd(PF$MAE), mean(PF$r), sd(PF$r), mean(CS5), sd(CS5), M1$hr_sd, CI7$dC_vs_pa[1], nrow(PW), sum(PW$fdr < 0.05)))
