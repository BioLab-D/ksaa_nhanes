# Figure 2 (KSAA quartiles and all-cause / cardiovascular / cancer mortality; one outcome per row: standardized 10-year absolute risk bars on the left, standardized cumulative incidence curves on the right) and Fig. S3 (subgroup forest plot).
# All results are survey-weighted and precomputed by weighted.R / weighted_boot.R (standardized curves, 10-year absolute risks, P for trend, subgroup HRs); this script only plots.
# Input: results/analysis/{analysis_data.rds, weighted/{figure2_std_curves, figure2_absrisk, figure2_ptrend, figureS3}.csv}
# Output: figures/Figure2, FigureS3
source(here::here("code", "analysis", "fig_common.R"))
d <- readRDS(file.path(V3, "analysis_data.rds")); d$ty <- d$t / 12
QL <- paste0("Q", 1:4); d$q <- factor(d$q, QL)
QPAL <- c(Q1 = "#9DB4C8", Q2 = "#4E8FD0", Q3 = "#1B3F7A", Q4 = "#C8372D")
qc <- quantile(d$gap, 0:4 / 4); fq <- function(x) f2(x, 1)
QLAB <- c(sprintf("Q1 (< %s years)", fq(qc[2])), sprintf("Q2 (%s to %s years)", fq(qc[2]), fq(qc[3])), sprintf("Q3 (%s to %s years)", fq(qc[3]), fq(qc[4])), sprintf("Q4 (> %s years)", fq(qc[4])))
TMAX <- 30                                              # x-axis limit (years)
FS <- 6                                                  # font size (pt)
thf <- th + theme(text = element_text(size = FS), axis.text = element_text(size = FS - .5, colour = INK2), axis.title = element_text(size = FS, colour = INK),
                  plot.title = element_text(size = 8, face = "bold"), plot.margin = margin(2, 4, 1, 2))

# ---------------- 1. Standardized curves ----------------
WD <- file.path(V3, "weighted")  # primary analyses are survey-weighted (weighted.R / weighted_boot.R -> results/analysis/weighted/)
STD <- read.csv(file.path(WD, "figure2_std_curves.csv")); STD$q <- factor(STD$q, QL); STD$time <- STD$month / 12; s10 <- STD[STD$month == 120, ]
cat(sprintf("standardized 10-year risk: all-cause %s; cardiovascular %s; cancer %s\n", paste(sprintf("%.1f", 100 * s10$allcause), collapse = "/"),
            paste(sprintf("%.1f", 100 * s10$cvd), collapse = "/"), paste(sprintf("%.1f", 100 * s10$cancer), collapse = "/")))

OUT <- c(allcause = "A", cvd = "B", cancer = "C")
YL <- c(allcause = "Cumulative all-cause mortality", cvd = "Cumulative incidence of cardiovascular death", cancer = "Cumulative incidence of cancer death")
YLS <- setNames(paste("Standardized", sub("^C", "c", YL)), names(YL))
YMAX <- c(allcause = 1, cvd = 0.5, cancer = 0.25); YMAXS <- c(allcause = 1, cvd = 0.5, cancer = 0.25)
# Layout: no grid, open axes, legend inside the panel
FS <- 5.2  # font size of the whole figure (except 8-pt bold panel letters)
FS_SUPP <- 6.5  # font size of the supplementary figure (scaled to about 0.87 in the docx)
mk_thk <- function(FS) theme_classic(base_family = "Arial", base_size = FS) +
  theme(axis.line = element_line(colour = INK, linewidth = .35), axis.ticks = element_line(colour = INK, linewidth = .35), axis.ticks.length = unit(2, "pt"),
        axis.text = element_text(colour = INK, size = FS), axis.title = element_text(colour = INK, size = FS, lineheight = .95),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)), plot.title.position = "plot",
        plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin(2, 6, 0, 2), legend.position = "none")
thk <- mk_thk(FS)
set_fs <- function(v) { FS <<- v; FR <<- v; thk <<- mk_thk(v) }   # temporarily switch the global font size (ftab / curve_grid read FS, FR, thk when called)
sr0 <- summary(survfit(Surv(ty, e) ~ q, d), times = seq(0, 30, 10), extend = TRUE)
RK <- data.frame(time = sr0$time, n = sr0$n.risk, q = factor(sub("^q=", "", sr0$strata), QL))
# x range fixed at 0–30.8 with coord_cartesian; the colored dots of the number-at-risk table are drawn outside the panel (clip off)
XS <- scale_x_continuous(breaks = seq(0, 30, 10), expand = c(0, 0))
XC <- function(ylim) coord_cartesian(xlim = c(0, TMAX + .8), ylim = ylim, expand = FALSE, clip = "off")
FR <- FS                                                                 # font size of numbers at risk
cap <- guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both"))
# y-axis titles state cumulative quantities: cardiovascular / cancer are cumulative incidence under competing risks
YT <- c(allcause = "Cumulative all-cause mortality (%)", cvd = "Cumulative incidence of cardiovascular death (%)", cancer = "Cumulative incidence of cancer death (%)")   # standardization is described in the legend
curve_plot <- function(o, legend = FALSE, S = NULL, ytitle = NULL, note = NULL, lab = NULL, marks = TRUE, xaxis = FALSE, inlegend = FALSE) {
  std <- is.null(S)
  if (std) { S <- STD[, c("time", "q", o)]; names(S)[3] <- "y" }
  S$y <- 100 * S$y; ym <- 100 * YMAXS[[o]]
  pc <- ggplot(S, aes(time, y, colour = q)) + geom_step(linewidth = .5, key_glyph = "point") +
    scale_colour_manual(values = QPAL, labels = QLAB, name = "KSAA quartile") + XS +
    scale_y_continuous(breaks = seq(0, ym, by = c(allcause = 25, cvd = 10, cancer = 5)[[o]]), expand = c(0, 0)) + XC(c(-ym * .03, ym)) +
    labs(title = if (is.null(lab)) OUT[[o]] else lab, x = if (xaxis) "Years since examination" else NULL,
         y = if (is.null(ytitle)) YT[[o]] else ytitle) + thk + cap
  if (!xaxis) pc <- pc + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())  # with xaxis = FALSE the x axis is given by the number-at-risk table below
  pc <- pc + annotate("segment", x = 0, xend = TMAX, y = -ym * .03, yend = -ym * .03, colour = INK, linewidth = .35)
  if (!is.null(note)) pc <- pc + annotate("text", x = TMAX, y = ym * .02, hjust = 1, vjust = 0, label = note, size = FS / .pt, family = "Arial", colour = INK)
  # dashed line at 10 years with the Q1 / Q4 values
  if (marks && std && o %in% c("allcause", "cvd")) {
    s10 <- STD[STD$month == 120, ]; v <- 100 * setNames(s10[[o]], as.character(s10$q))[c("Q1", "Q4")]
    pc <- pc + annotate("segment", x = 10, xend = 10, y = 0, yend = v[["Q4"]], linetype = "11", linewidth = .3, colour = MUTED) +
      annotate("point", x = 10, y = v, colour = QPAL[c("Q1", "Q4")], size = .9) +
      annotate("text", x = c(9.5, 10.5), y = c(v[["Q4"]] + .015 * ym, v[["Q1"]] - .015 * ym), hjust = c(1, 0), vjust = c(0, 1),
               label = sprintf("%.1f%%", v[c("Q4", "Q1")]), size = FS / .pt, family = "Arial", colour = c(QPAL[["Q4"]], "#6E8AA3"))  # the light Q1 color lacks contrast as text, so darken it
  }
  if (inlegend)  # legend in the top-left of each panel, one column
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
fmtp <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %s", f2(p, if (p < 0.1) 3 else 2))  # 0.0538 must not be shown as 0.05
EVC <- c(allcause = "e", cvd = "e_cvd", cancer = "e_cancer")
PW <- read.csv(file.path(WD, "figure2_ptrend.csv")); PTR <- setNames(PW$p_trend, PW$outcome)[names(EVC)]
cat(sprintf("P for trend (full cohort, age + sex; cause-specific for cardiovascular / cancer): %s\n", paste(signif(PTR, 3), collapse = " / ")))
ZM <- plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))  # zero outer margins of the combined plot (the default 5.5 pt accumulates)

# ---------------- 3. Subgroups ----------------
# sex, age, race, BMI, radiographic OA, WHO category of femoral neck BMD; weighted subgroup HRs and Cochran Q heterogeneity P values are written by weighted.R to weighted/figureS3.csv
EV <- c(allcause = "e", cvd = "e_cvd", cancer = "e_cancer")
SW <- read.csv(file.path(WD, "figureS3.csv"), stringsAsFactors = FALSE); SW$dth <- NA
REFH <- SW[SW$fam == "All participants", c("outcome", "n", "dth", "hr", "lo", "hi")]; SGA <- SW[SW$fam != "All participants", c("outcome", "fam", "lev", "n", "dth", "hr", "lo", "hi", "p_het")]
# total hip BMD quartiles are not shown (post hoc groups, cut points without clinical meaning, very different sex mix across pooled groups); BMD strata use only the femoral neck WHO categories
FAMD <- c("Sex", "Age", "Race and ethnicity", "BMI", "Radiographic OA (phase 2)", "Femoral neck BMD (WHO)")
FLAB <- c(BMI = "BMI, kg/m\u00b2")  # display name only
# first row "All participants": full-cohort HR with the same adjustment set, as the reference for the subgroups
mk <- function(fams, o) { r <- REFH[REFH$outcome == o, ]
  top <- data.frame(lab = "All participants", hdr = FALSE, n = cm(r$n), dth = cm(r$dth), est = r$hr, lo = r$lo, hi = r$hi, txt = ci(r$hr, r$lo, r$hi), p = "", all = TRUE)
  rbind(top, do.call(rbind, lapply(fams, function(fm) { g <- SGA[SGA$outcome == o & SGA$fam == fm, ]
  rbind(data.frame(lab = if (fm %in% names(FLAB)) FLAB[[fm]] else fm, hdr = TRUE, n = "", dth = "", est = NA, lo = NA, hi = NA, txt = "", p = if (g$p_het[1] > 0.99) "> 0.99" else if (g$p_het[1] < 0.001) "< 0.001" else f2(g$p_het[1], if (g$p_het[1] < 0.05) 3 else 2), all = FALSE),
        data.frame(lab = g$lev, hdr = FALSE, n = cm(g$n), dth = cm(g$dth), est = g$hr, lo = g$lo, hi = g$hi, txt = ci(g$hr, g$lo, g$hi), p = "", all = FALSE)) }))) }
# P is a Cochran Q heterogeneity test across separately fitted strata, not an interaction test
XL <- c(.9, 1.62); BR <- c(1, 1.2, 1.4, 1.6); SQ <- "#4E7FC4"  # only one reference line at HR = 1
# the all-cause and cardiovascular blocks share the left subgroup-name and n columns; positions within each block are in mm, and both forest areas are 30 mm wide
BW <- c(d = 96, e = 62); GAP <- 6
ftab <- function(D, span, letter, labels) {   # span: outcome name, used as the spanning header above the block
  FT <- FS
  w <- if (labels) BW[["d"]] else BW[["e"]]; o <- if (labels) 34 else 0          # width of the left label area (mm)
  fr <- c(n = 32.5, hr = o + 19, f0 = o + 22, f1 = o + 52, p = o + 56.9) / w  # column positions (fraction of the block width)
  nr <- nrow(D); D$row <- nr - seq_len(nr) + 1; top <- nr + 1.05
  L0 <- log10(XL[1]); L1 <- log10(XL[2]); S <- (L1 - L0) / (fr[["f1"]] - fr[["f0"]]); st <- L0 - fr[["f0"]] * S
  pos <- function(f) 10^(st + f * S)
  G <- D[!D$hdr, ]; H <- D[D$hdr, ]
  D$grp <- cumsum(D$hdr); PV <- do.call(rbind, lapply(split(D[D$grp > 0, ], D$grp[D$grp > 0]), function(z) data.frame(y = mean(z$row[!z$hdr]), p = z$p[z$hdr])))  # P centered on the level rows of each factor (full-cohort row has grp = 0 and no P)
  rule_end <- pos(fr[["hr"]])
  hd <- data.frame(x = c(if (labels) c(pos(0), pos(fr[["n"]])), pos(fr[["hr"]]), pos(fr[["p"]])),
                   lab = c(if (labels) c("Subgroup", "Number"), "HR (95% CI)", "P for\nheterogeneity"), hj = c(if (labels) c(0, 1), 1, .5))
  sp <- c(pos(fr[["hr"]] - 19 / w), pos(min(1, fr[["p"]] + 6 / w)))                    # spanning header from the left edge of the HR column to the right edge of the P column
  g <- ggplot() +
    annotate("segment", x = pos(0), xend = rule_end, y = nr + .62, yend = nr + .62, colour = INK, linewidth = .45) +
    geom_segment(data = rbind(H, D[D$all, names(H)]), aes(x = pos(0), xend = rule_end, y = row - .5, yend = row - .5), colour = "#9A9A9A", linewidth = .25) +
    annotate("segment", x = 1, xend = 1, y = .5, yend = nr + .3, colour = INK, linewidth = .4) +
    geom_errorbar(data = G, aes(xmin = lo, xmax = hi, y = row), orientation = "y", width = .32, colour = INK, linewidth = .3) +  # short caps at both ends
    geom_point(data = G[!G$all, ], aes(est, row), shape = 15, size = 1.5, colour = SQ) +
    geom_point(data = G[G$all, ], aes(est, row), shape = 18, size = 2.3, colour = INK) +   # full cohort: diamond
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
pdD <- ftab(DD, "All-cause mortality", " ", labels = TRUE)   # letter d sits in the empty left column
pdE <- ftab(DE, "Cardiovascular mortality", " ", labels = FALSE)
set_fs(5.2)

# ---------------- Assemble the figure (183 mm)----------------
# ---------------- 2b. Standardized 10-year absolute risk (bar chart) ----------------
AB <- subset(read.csv(file.path(WD, "figure2_absrisk.csv")), model == "M1 (full cohort)")  # full cohort (same basis as the curves)
s10b <- STD[STD$month == 120, ]
for (o in names(OUT)) stopifnot(all.equal(unlist(AB[AB$outcome == o, paste0("Q", 1:4)], use.names = FALSE),
                                          s10b[[o]][match(QL, as.character(s10b$q))], tolerance = 1e-8))   # bars and curves come from the same estimates
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
    coord_cartesian(xlim = c(.4, 4.6), ylim = c(0, ym), expand = FALSE, clip = "off") +  # about half a bar width of space at each end
    labs(title = if (is.null(lab)) OUTB[[o]] else lab, x = if (xtitle) "KSAA quartile" else NULL, y = YTB[[o]]) + thk +
    guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")) + theme(plot.margin = margin(6, 6, 2, 2)) +
    # guide_axis(cap) has no effect on a discrete axis, so the theme x-axis line is turned off and a segment at y = 0 spanning the panel is drawn instead, meeting the y axis at the origin
    theme(axis.line.x = element_blank()) + annotate("segment", x = .4, xend = 4.6, y = 0, yend = 0, colour = INK, linewidth = .35, lineend = "square")
  if (rdlab) p0 <- p0 + annotate("text", x = 4.5, y = ym * .99, hjust = 1, vjust = 1, label = rd, size = (FS - .5) / .pt, family = "Arial", colour = INK)
  p0
}
# one outcome per row, bars left and curves right; panel letters A–F in reading order
YTR <- c(allcause = "Cumulative mortality (%)", cvd = "Cumulative incidence (%)", cancer = "Cumulative incidence (%)")  # panels too short for long titles; the outcome is named by the y-axis title of the bar chart in the same row
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
NOM <- theme(plot.margin = margin(0, 0, 0, 0))  # remove the default 5.5-pt margin of the outer wrapper
warm <- function(x, w, h) ggsave(tempfile(fileext = ".png"), x, width = w, height = h, dpi = 300, bg = "white")  # the first render measures text sizes wrongly, so warm up first
p <- wrap_elements(full = top) + NOM + plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
warm(p, 5.2, 5.2)
save_fig(p, "Figure2", 5.2, 5.2)

# Supplementary figure: subgroup forest plot (all-cause | cardiovascular blocks share the left subgroup-name and n columns)
bottom <- pdD + plot_spacer() + pdE + plot_layout(widths = c(BW[["d"]], GAP, BW[["e"]])) + ZM
psg <- wrap_elements(full = bottom) + NOM + plot_annotation(theme = theme(plot.margin = margin(4, 12, 4, 4)))  # right margin so the header is not clipped
warm(psg, 7.2, 3.2)
save_fig(psg, "FigureS3", 7.2, 3.2)
for (o in names(EV)) { x <- unique(SGA[SGA$outcome == o, c("fam", "p_het")]); cat(sprintf("subgroups %-8s full-cohort HR %.2f; P heterogeneity %s\n", o, REFH$hr[REFH$outcome == o], paste(sprintf("%s %.2f", x$fam, x$p_het), collapse = "; "))) }
print(SGA[, c("outcome", "fam", "lev", "n", "dth", "hr", "lo", "hi")], digits = 3, row.names = FALSE)
