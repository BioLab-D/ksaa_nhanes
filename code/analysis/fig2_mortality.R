# Figure 2(KSAA 四分位与全因 / 心血管 / 癌症死亡;每行一个结局:左为标化 10 年绝对风险柱图,右为标化累积发生曲线)及补充图 S3(亚组森林图)。
# 均为调查加权结果,由 weighted.R / weighted_boot.R 预先算好(标化曲线、10 年绝对风险、P for trend、亚组 HR),本脚本只作图。
# 输入 results/analysis/{analysis_data.rds, weighted/{figure2_std_curves, figure2_absrisk, figure2_ptrend, figureS3}.csv}
# 输出 figures/Figure2、FigureS3
source(here::here("code", "analysis", "fig_common.R"))
d <- readRDS(file.path(V3, "analysis_data.rds")); d$ty <- d$t / 12
QL <- paste0("Q", 1:4); d$q <- factor(d$q, QL)
QPAL <- c(Q1 = "#9DB4C8", Q2 = "#4E8FD0", Q3 = "#1B3F7A", Q4 = "#C8372D")
qc <- quantile(d$gap, 0:4 / 4); fq <- function(x) f2(x, 1)
QLAB <- c(sprintf("Q1 (< %s years)", fq(qc[2])), sprintf("Q2 (%s to %s years)", fq(qc[2]), fq(qc[3])), sprintf("Q3 (%s to %s years)", fq(qc[3]), fq(qc[4])), sprintf("Q4 (> %s years)", fq(qc[4])))
TMAX <- 30                                              # 横轴上限(年)
FS <- 6                                                  # 成图字号(pt)
thf <- th + theme(text = element_text(size = FS), axis.text = element_text(size = FS - .5, colour = INK2), axis.title = element_text(size = FS, colour = INK),
                  plot.title = element_text(size = 8, face = "bold"), plot.margin = margin(2, 4, 1, 2))

# ---------------- 1. 标化曲线 ----------------
WD <- file.path(V3, "weighted")  # 主分析为调查加权(weighted.R / weighted_boot.R → results/analysis/weighted/)
STD <- read.csv(file.path(WD, "figure2_std_curves.csv")); STD$q <- factor(STD$q, QL); STD$time <- STD$month / 12; s10 <- STD[STD$month == 120, ]
cat(sprintf("标化 10 年:全因 %s;心血管 %s;癌症 %s\n", paste(sprintf("%.1f", 100 * s10$allcause), collapse = "/"),
            paste(sprintf("%.1f", 100 * s10$cvd), collapse = "/"), paste(sprintf("%.1f", 100 * s10$cancer), collapse = "/")))

OUT <- c(allcause = "A", cvd = "B", cancer = "C")
YL <- c(allcause = "Cumulative all-cause mortality", cvd = "Cumulative incidence of cardiovascular death", cancer = "Cumulative incidence of cancer death")
YLS <- setNames(paste("Standardized", sub("^C", "c", YL)), names(YL))
YMAX <- c(allcause = 1, cvd = 0.5, cancer = 0.25); YMAXS <- c(allcause = 1, cvd = 0.5, cancer = 0.25)
# 版式:无网格、开口坐标轴、图内图例
FS <- 5.2  # 全图字号(面板字母 8 pt 粗体除外)
FS_SUPP <- 6.5  # 补充图字号(docx 内约缩到 0.87 倍)
mk_thk <- function(FS) theme_classic(base_family = "Arial", base_size = FS) +
  theme(axis.line = element_line(colour = INK, linewidth = .35), axis.ticks = element_line(colour = INK, linewidth = .35), axis.ticks.length = unit(2, "pt"),
        axis.text = element_text(colour = INK, size = FS), axis.title = element_text(colour = INK, size = FS, lineheight = .95),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)), plot.title.position = "plot",
        plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin(2, 6, 0, 2), legend.position = "none")
thk <- mk_thk(FS)
set_fs <- function(v) { FS <<- v; FR <<- v; thk <<- mk_thk(v) }   # 临时切换全局字号(ftab / curve_grid 在调用时读 FS、FR、thk)
sr0 <- summary(survfit(Surv(ty, e) ~ q, d), times = seq(0, 30, 10), extend = TRUE)
RK <- data.frame(time = sr0$time, n = sr0$n.risk, q = factor(sub("^q=", "", sr0$strata), QL))
# x 范围用 coord_cartesian 定在 0–30.8,在险人数的彩色圆点画在坐标区外(clip off)
XS <- scale_x_continuous(breaks = seq(0, 30, 10), expand = c(0, 0))
XC <- function(ylim) coord_cartesian(xlim = c(0, TMAX + .8), ylim = ylim, expand = FALSE, clip = "off")
FR <- FS                                                                 # 在险人数字号
cap <- guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both"))
# 纵轴写明为累积量:心血管 / 癌症为竞争风险下的累积发生率
YT <- c(allcause = "Cumulative all-cause mortality (%)", cvd = "Cumulative incidence of cardiovascular death (%)", cancer = "Cumulative incidence of cancer death (%)")   # 标化见图注
curve_plot <- function(o, legend = FALSE, S = NULL, ytitle = NULL, note = NULL, lab = NULL, marks = TRUE, xaxis = FALSE, inlegend = FALSE) {
  std <- is.null(S)
  if (std) { S <- STD[, c("time", "q", o)]; names(S)[3] <- "y" }
  S$y <- 100 * S$y; ym <- 100 * YMAXS[[o]]
  pc <- ggplot(S, aes(time, y, colour = q)) + geom_step(linewidth = .5, key_glyph = "point") +
    scale_colour_manual(values = QPAL, labels = QLAB, name = "KSAA quartile") + XS +
    scale_y_continuous(breaks = seq(0, ym, by = c(allcause = 25, cvd = 10, cancer = 5)[[o]]), expand = c(0, 0)) + XC(c(-ym * .03, ym)) +
    labs(title = if (is.null(lab)) OUT[[o]] else lab, x = if (xaxis) "Years since examination" else NULL,
         y = if (is.null(ytitle)) YT[[o]] else ytitle) + thk + cap
  if (!xaxis) pc <- pc + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())  # xaxis = FALSE 时横轴由下方在险人数表给出
  pc <- pc + annotate("segment", x = 0, xend = TMAX, y = -ym * .03, yend = -ym * .03, colour = INK, linewidth = .35)
  if (!is.null(note)) pc <- pc + annotate("text", x = TMAX, y = ym * .02, hjust = 1, vjust = 0, label = note, size = FS / .pt, family = "Arial", colour = INK)
  # 标化曲线在 10 年处加虚线与 Q1 / Q4 的数值
  if (marks && std && o %in% c("allcause", "cvd")) {
    s10 <- STD[STD$month == 120, ]; v <- 100 * setNames(s10[[o]], as.character(s10$q))[c("Q1", "Q4")]
    pc <- pc + annotate("segment", x = 10, xend = 10, y = 0, yend = v[["Q4"]], linetype = "11", linewidth = .3, colour = MUTED) +
      annotate("point", x = 10, y = v, colour = QPAL[c("Q1", "Q4")], size = .9) +
      annotate("text", x = c(9.5, 10.5), y = c(v[["Q4"]] + .015 * ym, v[["Q1"]] - .015 * ym), hjust = c(1, 0), vjust = c(0, 1),
               label = sprintf("%.1f%%", v[c("Q4", "Q1")]), size = FS / .pt, family = "Arial", colour = c(QPAL[["Q4"]], "#6E8AA3"))  # Q1 浅色作文字对比度不够,加深
  }
  if (inlegend)  # 图例在每格左上角、单列
    return(pc + guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.2))) +
      theme(legend.position = "inside", legend.position.inside = c(.015, .99), legend.justification = c(0, 1),
            legend.title = element_blank(), legend.text = element_text(size = FS), legend.background = element_blank(),
            legend.key = element_blank(), legend.key.size = unit(5, "pt"), legend.key.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0)))
  pc + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
    theme(legend.title = element_text(size = FS), legend.text = element_text(size = FS),
          legend.key.size = unit(6, "pt"), legend.key.spacing.x = unit(5, "pt"), legend.margin = margin(0, 0, 4, 0))
}
risk_plot <- function() ggplot(RK, aes(time, q)) + geom_text(aes(label = cm(n)), size = FR / .pt, family = "Arial", colour = INK) +
  geom_point(aes(x = -2.3, colour = q), size = 1) + scale_colour_manual(values = QPAL, guide = "none") + XS + XC(c(.4, 4.6)) +
  scale_y_discrete(limits = rev(QL)) +
  labs(x = "Years since examination", y = NULL, subtitle = "Number at risk") + thk + guides(x = guide_axis(cap = "both")) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        plot.subtitle = element_text(size = FR, colour = INK, margin = margin(t = 2, b = 1), hjust = 0), plot.margin = margin(0, 6, 2, 2))
fmtp <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %s", f2(p, if (p < 0.1) 3 else 2))  # 0.0538 不能显示成 0.05
EVC <- c(allcause = "e", cvd = "e_cvd", cancer = "e_cancer")
PW <- read.csv(file.path(WD, "figure2_ptrend.csv")); PTR <- setNames(PW$p_trend, PW$outcome)[names(EVC)]
cat(sprintf("P for trend(全队列,年龄 + 性别;心血管 / 癌症为死因别):%s\n", paste(signif(PTR, 3), collapse = " / ")))
ZM <- plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))  # 拼图外边距归零(默认 5.5 pt 会层层叠加)

# ---------------- 3. 亚组 ----------------
# 性别、年龄、种族、BMI、影像学 OA、股骨颈 BMD 的 WHO 分类;加权亚组 HR 与 Cochran Q 异质性 P 由 weighted.R 写入 weighted/figureS3.csv
EV <- c(allcause = "e", cvd = "e_cvd", cancer = "e_cancer")
SW <- read.csv(file.path(WD, "figureS3.csv"), stringsAsFactors = FALSE); SW$dth <- NA
REFH <- SW[SW$fam == "All participants", c("outcome", "n", "dth", "hr", "lo", "hi")]; SGA <- SW[SW$fam != "All participants", c("outcome", "fam", "lev", "n", "dth", "hr", "lo", "hi", "p_het")]
# 全髋 BMD 四分位不画(事后分组、切点无临床含义、合并切分后各组性别构成差异大);骨密度分层只留股骨颈 WHO 分类
FAMD <- c("Sex", "Age", "Race and ethnicity", "BMI", "Radiographic OA (phase 2)", "Femoral neck BMD (WHO)")
FLAB <- c(BMI = "BMI, kg/m\u00b2")  # 只改显示名
# 首行 "All participants":全队列、同一组调整变量的 HR,作为各亚组的参照
mk <- function(fams, o) { r <- REFH[REFH$outcome == o, ]
  top <- data.frame(lab = "All participants", hdr = FALSE, n = cm(r$n), dth = cm(r$dth), est = r$hr, lo = r$lo, hi = r$hi, txt = ci(r$hr, r$lo, r$hi), p = "", all = TRUE)
  rbind(top, do.call(rbind, lapply(fams, function(fm) { g <- SGA[SGA$outcome == o & SGA$fam == fm, ]
  rbind(data.frame(lab = if (fm %in% names(FLAB)) FLAB[[fm]] else fm, hdr = TRUE, n = "", dth = "", est = NA, lo = NA, hi = NA, txt = "", p = if (g$p_het[1] > 0.99) "> 0.99" else if (g$p_het[1] < 0.001) "< 0.001" else f2(g$p_het[1], if (g$p_het[1] < 0.05) 3 else 2), all = FALSE),
        data.frame(lab = g$lev, hdr = FALSE, n = cm(g$n), dth = cm(g$dth), est = g$hr, lo = g$lo, hi = g$hi, txt = ci(g$hr, g$lo, g$hi), p = "", all = FALSE)) }))) }
# P 为各层分别拟合后的 Cochran Q 异质性检验,不是交互项检验
XL <- c(.9, 1.62); BR <- c(1, 1.2, 1.4, 1.6); SQ <- "#4E7FC4"  # 只画 HR = 1 一条参考线
# 全因、心血管两块共用左侧亚组名与人数列;各块内位置以毫米给出,两块森林区物理宽度相同(30 mm)
BW <- c(d = 96, e = 62); GAP <- 6
ftab <- function(D, span, letter, labels) {   # span:结局名,作为块上方的跨列标题
  FT <- FS
  w <- if (labels) BW[["d"]] else BW[["e"]]; o <- if (labels) 34 else 0          # 左侧标签区宽度(mm)
  fr <- c(n = 32.5, hr = o + 19, f0 = o + 22, f1 = o + 52, p = o + 56.9) / w  # 各列位置(占本块宽度的比例)
  nr <- nrow(D); D$row <- nr - seq_len(nr) + 1; top <- nr + 1.05
  L0 <- log10(XL[1]); L1 <- log10(XL[2]); S <- (L1 - L0) / (fr[["f1"]] - fr[["f0"]]); st <- L0 - fr[["f0"]] * S
  pos <- function(f) 10^(st + f * S)
  G <- D[!D$hdr, ]; H <- D[D$hdr, ]
  D$grp <- cumsum(D$hdr); PV <- do.call(rbind, lapply(split(D[D$grp > 0, ], D$grp[D$grp > 0]), function(z) data.frame(y = mean(z$row[!z$hdr]), p = z$p[z$hdr])))  # P 放在各因素水平行的正中(全队列行 grp = 0,无 P)
  rule_end <- pos(fr[["hr"]])
  hd <- data.frame(x = c(if (labels) c(pos(0), pos(fr[["n"]])), pos(fr[["hr"]]), pos(fr[["p"]])),
                   lab = c(if (labels) c("Subgroup", "Number"), "HR (95% CI)", "P for\nheterogeneity"), hj = c(if (labels) c(0, 1), 1, .5))
  sp <- c(pos(fr[["hr"]] - 19 / w), pos(min(1, fr[["p"]] + 6 / w)))                    # 跨列标题范围:HR 列左缘到 P 列右缘
  g <- ggplot() +
    annotate("segment", x = pos(0), xend = rule_end, y = nr + .62, yend = nr + .62, colour = INK, linewidth = .45) +
    geom_segment(data = rbind(H, D[D$all, names(H)]), aes(x = pos(0), xend = rule_end, y = row - .5, yend = row - .5), colour = "#9A9A9A", linewidth = .25) +
    annotate("segment", x = 1, xend = 1, y = .5, yend = nr + .3, colour = INK, linewidth = .4) +
    geom_errorbar(data = G, aes(xmin = lo, xmax = hi, y = row), orientation = "y", width = .32, colour = INK, linewidth = .3) +  # 两端带短竖线
    geom_point(data = G[!G$all, ], aes(est, row), shape = 15, size = 1.5, colour = SQ) +
    geom_point(data = G[G$all, ], aes(est, row), shape = 18, size = 2.3, colour = INK) +   # 全队列:菱形
    geom_text(data = G, aes(x = pos(fr[["hr"]]), y = row, label = txt), hjust = 1, size = FT / .pt, family = "Arial", colour = INK) +
    geom_text(data = PV, aes(x = pos(fr[["p"]]), y = y, label = p), hjust = .5, size = FT / .pt, family = "Arial", colour = INK) +
    geom_text(data = hd, aes(x = x, y = top, label = lab, hjust = hj), vjust = 0, lineheight = .95, size = FT / .pt, family = "Arial", colour = INK) +
    annotate("text", x = 10^mean(log10(sp)), y = top + 1.95, label = span, vjust = 0, size = FT / .pt, family = "Arial", fontface = "bold", colour = INK) +
    annotate("segment", x = sp[1], xend = sp[2], y = top + 1.7, yend = top + 1.7, colour = INK, linewidth = .3)
  if (labels) g <- g +
    geom_text(data = D, aes(x = ifelse(hdr | all, pos(0), pos(.015)), y = row, label = lab), hjust = 0, size = FT / .pt, family = "Arial", colour = INK) +
    geom_text(data = D, aes(x = pos(fr[["n"]]), y = row, label = n), hjust = 1, size = FT / .pt, family = "Arial", colour = INK)
  g + annotate("segment", x = XL[1], xend = XL[2], y = .3, yend = .3, colour = INK, linewidth = .35) +
    annotate("text", x = 10^((L0 + L1) / 2), y = -1.35, label = "HR per SD of KSAA", size = FT / .pt, family = "Arial", colour = INK) +
    scale_x_log10(limits = c(pos(0), pos(1)), breaks = BR, labels = c("1.0", "1.2", "1.4", "1.6"), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + coord_cartesian(ylim = c(.3, top + 2.45), clip = "off") +
    labs(title = letter, x = NULL, y = NULL) + thk +
    theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_line(colour = INK, linewidth = .35),
          axis.text.x = element_text(size = FT), plot.margin = margin(2, 0, 15, 0))
}
DD <- mk(FAMD, "allcause"); DE <- mk(FAMD, "cvd")
set_fs(FS_SUPP)
pdD <- ftab(DD, "All-cause mortality", " ", labels = TRUE)   # 字母 d 在左侧空列里
pdE <- ftab(DE, "Cardiovascular mortality", " ", labels = FALSE)
set_fs(5.2)

# ---------------- 组图(183 mm)----------------
# ---------------- 2b. 标化 10 年绝对风险(柱图) ----------------
AB <- subset(read.csv(file.path(WD, "figure2_absrisk.csv")), model == "M1 (full cohort)")  # 全队列(与曲线同口径)
s10b <- STD[STD$month == 120, ]
for (o in names(OUT)) stopifnot(all.equal(unlist(AB[AB$outcome == o, paste0("Q", 1:4)], use.names = FALSE),
                                          s10b[[o]][match(QL, as.character(s10b$q))], tolerance = 1e-8))   # 柱与曲线同源
OUTB <- c(allcause = "D", cvd = "E", cancer = "F")
YMAXB <- c(allcause = 60, cvd = 30, cancer = 15); YSTEP <- c(allcause = 20, cvd = 10, cancer = 5)
YTB <- c(allcause = "10-year all-cause mortality (%)", cvd = "10-year cardiovascular death (%)", cancer = "10-year cancer death (%)")
bar_plot <- function(o, lab = NULL, xtitle = TRUE, rdlab = TRUE) {
  a <- AB[AB$outcome == o, ]; ym <- YMAXB[[o]]
  B <- data.frame(q = factor(QL, QL), y = 100 * unlist(a[, paste0("Q", 1:4)], use.names = FALSE))
  rd <- sprintf("Q4 \u2212 Q1  %s pp (%s to %s)", f2(100 * a$rd, 1), f2(100 * a$rd_lo, 1), f2(100 * a$rd_hi, 1))
  p0 <- ggplot(B, aes(q, y, fill = q)) + geom_col(width = .58) +
    geom_text(aes(label = sprintf("%.1f", y)), vjust = -.45, size = (FS - .5) / .pt, family = "Arial", colour = INK) +
    scale_fill_manual(values = QPAL, guide = "none") +
    scale_y_continuous(breaks = seq(0, ym, by = YSTEP[[o]]), expand = c(0, 0)) +
    coord_cartesian(xlim = c(.4, 4.6), ylim = c(0, ym), expand = FALSE, clip = "off") +  # 两端各留约半个柱宽
    labs(title = if (is.null(lab)) OUTB[[o]] else lab, x = if (xtitle) "KSAA quartile" else NULL, y = YTB[[o]]) + thk +
    guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) + theme(plot.margin = margin(6, 6, 2, 2)) +
    # 离散轴上 guide_axis(cap) 不生效,故关掉主题横轴线,改画一条铺满面板、落在 y = 0 的线段,与纵轴在原点相交
    theme(axis.line.x = element_blank()) + annotate("segment", x = .4, xend = 4.6, y = 0, yend = 0, colour = INK, linewidth = .35, lineend = "square")
  if (rdlab) p0 <- p0 + annotate("text", x = 4.5, y = ym * .99, hjust = 1, vjust = 1, label = rd, size = (FS - .5) / .pt, family = "Arial", colour = INK)
  p0
}
# 每行一个结局,左柱图右曲线;面板字母按阅读顺序 A–F
YTR <- c(allcause = "Cumulative mortality (%)", cvd = "Cumulative incidence (%)", cancer = "Cumulative incidence (%)")  # 格子太矮放不下长标题;结局名由同行柱图的纵轴标题给出
fig2_grid <- function() {
  os <- c("allcause", "cvd", "cancer"); LB <- LETTERS[1:6]
  cells <- list()
  for (i in seq_along(os)) {
    o <- os[[i]]
    cells <- c(cells, list(bar_plot(o, lab = LB[2 * i - 1], rdlab = FALSE),
                           curve_plot(o, lab = LB[2 * i], ytitle = YTR[[o]], marks = FALSE, xaxis = TRUE, inlegend = TRUE,
                                      note = sub("^P ", "P for trend ", fmtp(PTR[[o]])))))
  }
  wrap_plots(cells, ncol = 2, widths = c(1, 1.1)) + ZM
}
top <- fig2_grid()
NOM <- theme(plot.margin = margin(0, 0, 0, 0))  # 去掉外层包裹的默认 5.5 pt 边距
warm <- function(x, w, h) ggsave(tempfile(fileext = ".png"), x, width = w, height = h, dpi = 300, bg = "white")  # 首次渲染文字尺寸会量错,先热身
p <- wrap_elements(full = top) + NOM + plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
warm(p, 5.2, 5.2)
save_fig(p, "Figure2", 5.2, 5.2)

# 补充图:亚组森林图(全因 | 心血管两块共用左侧亚组名与人数列)
bottom <- pdD + plot_spacer() + pdE + plot_layout(widths = c(BW[["d"]], GAP, BW[["e"]])) + ZM
psg <- wrap_elements(full = bottom) + NOM + plot_annotation(theme = theme(plot.margin = margin(4, 12, 4, 4)))  # 右侧表头留边距,免被裁
warm(psg, 7.2, 3.2)
save_fig(psg, "FigureS3", 7.2, 3.2)
for (o in names(EV)) { x <- unique(SGA[SGA$outcome == o, c("fam", "p_het")]); cat(sprintf("亚组 %-8s 全队列 HR %.2f;P 异质性 %s\n", o, REFH$hr[REFH$outcome == o], paste(sprintf("%s %.2f", x$fam, x$p_het), collapse = "; "))) }
print(SGA[, c("outcome", "fam", "lev", "n", "dth", "hr", "lo", "hi")], digits = 3, row.names = FALSE)
