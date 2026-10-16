# Figure 3: incremental prognostic value of KSAA beyond Phenotypic Age acceleration.
#   A, per-SD HRs with both in the model; B, Harrell's C and R²_D before and after adding KSAA, with relative gains;
#   C-E, HRs in cells cross-classified by the quartiles of both (reference Q1 x Q1).
# n = 4,522 (PhenoAge available). C and R²_D are apparent values (no optimism correction); CIs from replicate-weight bootstrap.
# Input: results/analysis/weighted/{figure3a, figure3b_cindex, figure3_cells, figure3_within}.csv
#   or, with KSAA_UNWEIGHTED=1, results/analysis/{phenoage_combined_cindex, phenoage_combined_improvement, phenoage_wald,
#   joint_cells, joint_within}.csv
# Output: results/figures/Figure3.{png,pdf}
source(here::here("code", "analysis", "fig_common.R"))
GREEN <- "#2e7d55"   # color of the "gain from adding KSAA" annotation, not a data encoding
OUT3 <- c("All-cause", "Cardiovascular", "Cancer")
norm_out <- function(x) factor(sub(" mortality$", "", sub("^CVD", "Cardiovascular", x)), OUT3)

CI7 <- read.csv(file.path(V3, "phenoage_combined_cindex.csv")); CI7$o <- norm_out(CI7$outcome)
IMP <- read.csv(file.path(V3, "phenoage_combined_improvement.csv")); IMP$o <- norm_out(IMP$outcome)
WD <- read.csv(file.path(V3, "phenoage_wald.csv")); WD$o <- norm_out(WD$outcome)
JC <- read.csv(file.path(V3, "joint_cells.csv")); JC <- JC[JC$panel %in% c("main", "main_cvd", "main_cancer"), ]   # three heat maps: all-cause / cardiovascular / cancer
JW <- read.csv(file.path(V3, "joint_within.csv")); JW <- JW[JW$panel %in% c("main", "main_cvd"), ]   # per-SD HR of KSAA within each quartile (for the legend and console output)
# survey-weighted results by default; KSAA_UNWEIGHTED=1 draws the unweighted version (same sample and quartile cut points)
if (Sys.getenv("KSAA_UNWEIGHTED") != "1") { WD_ <- file.path(V3, "weighted")
  OM <- c(allcause = "All-cause mortality", cvd = "CVD mortality", cancer = "Cancer mortality"); PN <- c(allcause = "main", cvd = "main_cvd", cancer = "main_cancer")
  x <- read.csv(file.path(WD_, "figure3a.csv")); x$outcome <- OM[x$outcome]; x$o <- norm_out(x$outcome); WD <- x
  b <- read.csv(file.path(WD_, "figure3b_cindex.csv"))
  CI7 <- data.frame(outcome = OM[b$outcome], C_pa = b$C0, C_pa_lo = b$C0_lo, C_pa_hi = b$C0_hi, C_cb = b$C1, C_cb_lo = b$C1_lo, C_cb_hi = b$C1_hi,
                    dC_vs_pa = b$dC, dC_vs_pa_lo = b$dC_lo, dC_vs_pa_hi = b$dC_hi, relC = b$relC, relC_lo = b$relC_lo, relC_hi = b$relC_hi,
                    R2_pa = b$R0, R2_pa_lo = b$R0_lo, R2_pa_hi = b$R0_hi, R2_cb = b$R1, R2_cb_lo = b$R1_lo, R2_cb_hi = b$R1_hi, relR2 = b$relR2, relR2_lo = b$relR2_lo, relR2_hi = b$relR2_hi)
  CI7$o <- norm_out(CI7$outcome)
  c3 <- read.csv(file.path(WD_, "figure3_cells.csv")); c3$hr_lo[is.na(c3$hr_lo)] <- 1; c3$hr_hi[is.na(c3$hr_hi)] <- 1
  JC <- data.frame(panel = PN[c3$outcome], y = c3$y, x = c3$x, risk = NA, hr = c3$hr, hr_lo = c3$hr_lo, hr_hi = c3$hr_hi, ref = c3$y == "Q1" & c3$x == "Q1", n = c3$n, deaths = c3$deaths)
  w3 <- read.csv(file.path(WD_, "figure3_within.csv")); w3 <- w3[w3$outcome %in% c("allcause", "cvd"), ]
  JW <- data.frame(panel = PN[w3$outcome], stratum = w3$stratum, n = NA, deaths = NA, hr = w3$hr, lo = w3$lo, hi = w3$hi) }

# Both measures are residuals on chronological age and are reported per SD.
FS <- 5.2
AL <- as.numeric(Sys.getenv("FIG3_AL", "24"))   # left indent of panel letter A (pt): aligned with letter C in the bottom row
PA <- "Phenotypic Age acceleration"; PK <- "Phenotypic Age acceleration + KSAA"; KS <- "KSAA"
NAVY <- "#0d366b"
PAL <- setNames(c(ORANGE, NAVY, NAVY), c(PA, PK, KS))
SHP <- setNames(c(16, 18, 15), c(PA, PK, KS))
SZ <- setNames(c(1.3, 1.9, 1.35), c(PA, PK, KS))
thk <- theme_classic(base_family = "Arial", base_size = FS) +
  theme(axis.line = element_blank(), panel.border = element_rect(colour = MUTED, fill = NA, linewidth = .35),   # mid-gray border and tick marks
        axis.ticks = element_line(colour = MUTED, linewidth = .35), axis.ticks.length = unit(2, "pt"),
        axis.text = element_text(colour = INK, size = FS), axis.title = element_text(colour = INK, size = FS),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)), plot.title.position = "plot",
        legend.position = "top", legend.location = "plot", legend.justification = "center", legend.direction = "horizontal", legend.title = element_blank(),   # one row, centered on the whole panel (including the y axis)
        legend.text = element_text(colour = INK, size = FS), legend.key.size = unit(6, "pt"), legend.key.spacing.x = unit(5, "pt"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"), plot.background = element_blank(), panel.background = element_blank(),
        plot.margin = margin(2, 6, 2, 2))
cap <- NULL   # do not cap the axes of the estimate plots (error bars would extend beyond capped axes)
num <- function(x) { r <- round(x); ifelse(r == 0, "0", mi(as.character(r))) }                       # no sign on 0
sgn <- function(x) { r <- round(x); ifelse(r == 0, "0", ifelse(r > 0, paste0("+", r), mi(as.character(r)))) }
pct <- function(est, d) sprintf("%s%%", ifelse(round(est, d) > 0, paste0("+", f2(est, d)), f2(est, d)))

# y-axis upper limit: the annotation area (2% top margin, about 7% per line, plus 3%) must lie above the highest CI limit
ytop <- function(hi, bot, lines, step) { z <- 0.05 + 0.07 * lines; ceiling(((max(hi) - z * bot) / (1 - z)) / step) * step }
mkpanel <- function(L, ann, title, ylab, ylim, breaks, ref = NA, ann_y = ylim[2] - 0.02 * diff(ylim)) {   # 2% between the top of the green text and the top border
  lev <- unique(L$series); L$series <- factor(L$series, lev)
  L$x <- as.integer(L$o) + ifelse(L$series == lev[1], -0.15, 0.15)
  p <- ggplot(L, aes(x, est, colour = series))
  if (!is.na(ref)) p <- p + geom_hline(yintercept = ref, linetype = "22", colour = MUTED, linewidth = 0.3)
  p + geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.45, show.legend = FALSE) +
    geom_point(aes(shape = series, size = series)) +
    geom_text(data = ann, aes(x = as.integer(o), y = ann_y, label = lab), inherit.aes = FALSE,
              size = FS / .pt, family = "Arial", colour = GREEN, lineheight = 1.0, vjust = 1) +
    scale_colour_manual(values = PAL[lev], name = NULL) + scale_shape_manual(values = SHP[lev], name = NULL) +
    scale_size_manual(values = SZ[lev], guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = OUT3, limits = c(0.5, 3.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = breaks, labels = function(v) mi(formatC(v, format = "f", digits = if (max(breaks) < 1) 2 else 1)), expand = c(0, 0)) +
    coord_cartesian(ylim = ylim, clip = "off") +
    labs(title = title, x = NULL, y = ylab) + thk +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.6)))
}

# ---- Harrell's C (not plotted) ----
La <- rbind(data.frame(o = CI7$o, est = CI7$C_pa, lo = CI7$C_pa_lo, hi = CI7$C_pa_hi, series = PA),
            data.frame(o = CI7$o, est = CI7$C_cb, lo = CI7$C_cb_lo, hi = CI7$C_cb_hi, series = PK))
Aa <- data.frame(o = CI7$o, lab = pct(CI7$relC, 2))
pa <- mkpanel(La, Aa, "A", "Harrell's C (95% CI)", c(0.63, ytop(La$hi, 0.63, 1, 0.005)), seq(0.65, 0.77, 0.03))

# ---- explained variation R²_D (Royston–Sauerbrei; not plotted) ----
Lr <- rbind(data.frame(o = CI7$o, est = CI7$R2_pa, lo = CI7$R2_pa_lo, hi = CI7$R2_pa_hi, series = PA),
            data.frame(o = CI7$o, est = CI7$R2_cb, lo = CI7$R2_cb_lo, hi = CI7$R2_cb_hi, series = PK))
Ar <- data.frame(o = CI7$o, lab = pct(CI7$relR2, 1))
pr <- mkpanel(Lr, Ar, "B", "Explained variation, R² (95% CI)", c(0.08, ytop(Lr$hi, 0.08, 1, 0.01)), seq(0.1, 0.4, 0.1))   # no plotmath subscript (cairo_pdf fails with a CID font error); R²_D is defined in the legend

# ---- A: per-SD HRs with both in the model, for three outcomes (horizontal forest plot, values on the right) ----
XLO <- 0.88; XHI <- 1.68; XT <- 1.73; XMAX <- 2.35   # forest area + value column on the right (log axis)
XBRK <- c(0.9, 1, 1.2, 1.4, 1.6)   # axis line is capped at the outer ticks, so the lowest tick must be below the lowest CI limit (KSAA, cancer: 0.95)
hrtxt <- function(e, l, h) sprintf("%s (%s\u2013%s)", f2(e, 2), f2(l, 2), f2(h, 2))
fbase <- function(df, ylim) ggplot(df, aes(hr, y)) +
  annotate("rect", xmin = XT / 1.04, xmax = XMAX, ymin = ylim[1], ymax = ylim[2], fill = "white", colour = NA) +
  annotate("segment", x = 1, xend = 1, y = ylim[1], yend = ylim[2], linetype = "22", colour = MUTED, linewidth = 0.3) +   # a segment rather than geom_vline, so it meets the x-axis line
  scale_x_continuous(trans = "log", breaks = XBRK, labels = function(v) mi(formatC(v, format = "f", digits = 1)), expand = c(0, 0)) +
  coord_cartesian(xlim = c(XLO, XMAX), ylim = ylim, expand = FALSE, clip = "off") +
  thk + guides(x = guide_axis(cap = "both")) +   # axis line capped between the first and last ticks (not crossing the value column) and joined to the ticks
  theme(axis.line.x = element_line(colour = INK, linewidth = 0.35),
        axis.line.y = element_line(colour = INK, linewidth = 0.35),
        axis.ticks.y = element_blank(), panel.border = element_blank(), panel.grid = element_blank())

LA <- rbind(data.frame(o = WD$o, hr = WD$hr_pa, lo = WD$lo_pa, hi = WD$hi_pa, series = PA),
            data.frame(o = WD$o, hr = WD$hr_ksaa, lo = WD$lo_ksaa, hi = WD$hi_ksaa, series = KS))
LA$series <- factor(LA$series, c(PA, KS)); LA$row <- 4 - as.integer(factor(LA$o, OUT3))
LA$y <- LA$row + ifelse(LA$series == PA, 0.18, -0.18); LA$txt <- hrtxt(LA$hr, LA$lo, LA$hi)
pc <- fbase(LA, c(0.45, 3.75)) +
  geom_linerange(aes(xmin = lo, xmax = hi, colour = series), linewidth = 0.45) +
  geom_point(aes(colour = series, shape = series, size = series)) +
  geom_text(aes(x = XT, label = txt, colour = series), hjust = 0, size = FS / .pt, family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = PAL[c(PA, KS)], name = NULL) + scale_shape_manual(values = SHP[c(PA, KS)], name = NULL) +
  scale_size_manual(values = SZ[c(PA, KS)], guide = "none") +
  scale_y_continuous(breaks = 3:1, labels = OUT3, expand = c(0, 0)) +
  labs(title = "A", x = "HR per SD (95% CI)", y = NULL) +
  theme(plot.title = element_text(margin = margin(l = AL, b = 1)),   # shift the letter right to align with letter C in the bottom row
        axis.title.x = element_text(hjust = log(XHI / XLO) / (2 * log(XMAX / XLO)), margin = margin(t = 3.2)))   # centered on the forest area; top margin aligns it with the axis title of B

QLV <- paste0("Q", 1:4)
qhet_p <- function(w) { b <- log(w$hr); se <- (log(w$hi) - log(w$lo)) / (2 * qnorm(.975)); wt <- 1 / se^2   # Cochran Q, df = 3
  pchisq(sum(wt * (b - sum(wt * b) / sum(wt))^2), length(b) - 1, lower.tail = FALSE) }
# ---- B: Harrell's C and R²_D of the base model (age + sex + PhenoAge acceleration) and with KSAA; relative gains on the right ----
#   CIs from 200 bootstrap replicates, refitting both models in each replicate.
MC <- "Harrell\u2019s C"; MR <- "Explained variation, R\u00b2"
mkb <- function(est, lo, hi, m, ser) data.frame(o = CI7$o, est = est, lo = lo, hi = hi, m = m, s = ser)
LB2 <- rbind(mkb(CI7$C_pa, CI7$C_pa_lo, CI7$C_pa_hi, MC, PA), mkb(CI7$C_cb, CI7$C_cb_lo, CI7$C_cb_hi, MC, PK),
             mkb(CI7$R2_pa, CI7$R2_pa_lo, CI7$R2_pa_hi, MR, PA), mkb(CI7$R2_cb, CI7$R2_cb_lo, CI7$R2_cb_hi, MR, PK))
LB2$m <- factor(LB2$m, c(MC, MR)); LB2$s <- factor(LB2$s, c(PA, PK))
LB2$y <- (4 - as.integer(factor(LB2$o, OUT3))) + ifelse(LB2$s == PA, 0.16, -0.16)
ANN <- rbind(data.frame(o = CI7$o, m = MC, lab = pct(CI7$relC, 2)), data.frame(o = CI7$o, m = MR, lab = pct(CI7$relR2, 1)))
ANN$m <- factor(ANN$m, c(MC, MR)); ANN$y <- 4 - as.integer(factor(ANN$o, OUT3))
pd <- ggplot(LB2, aes(est, y, colour = s)) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.45) +
  geom_point(aes(shape = s, size = s)) +
  geom_text(data = ANN, aes(x = Inf, y = y, label = lab), inherit.aes = FALSE, hjust = 1.02,
            size = (FS - 0.5) / .pt, family = "Arial", colour = PAL[[PK]]) +
  facet_wrap(~ m, nrow = 1, scales = "free_x", strip.position = "bottom") +
  scale_colour_manual(values = PAL[c(PA, PK)], name = NULL) + scale_shape_manual(values = SHP[c(PA, PK)], name = NULL) +
  scale_size_manual(values = SZ[c(PA, PK)], guide = "none") +
  scale_x_continuous(breaks = scales::breaks_pretty(3), labels = function(v) mi(formatC(v, format = "f", digits = 2)),
                     expand = expansion(mult = c(0.10, 0.34))) +   # room for the percentages on the right
  scale_y_continuous(breaks = 3:1, labels = OUT3, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.45, 3.75), clip = "off") +
  labs(title = "B", x = NULL, y = NULL) + thk +
  theme(axis.line.x = element_line(colour = INK, linewidth = 0.35),
        axis.line.y = element_line(colour = INK, linewidth = 0.35), axis.ticks.y = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank(), panel.spacing.x = unit(10, "pt"),
        strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = FS, colour = INK, margin = margin(t = 2)))

LB <- rbind(transform(JW[JW$panel == "main", ], blk = 1), transform(JW[JW$panel == "main_cvd", ], blk = 2))
stopifnot(nrow(LB) == 8)
LB$y <- ifelse(LB$blk == 1, 9.4, 4.0) - (as.integer(factor(LB$stratum, QLV)) - 1)
LB$txt <- hrtxt(LB$hr, LB$lo, LB$hi)
PH <- sapply(c("main", "main_cvd"), function(tg) qhet_p(JW[JW$panel == tg, ]))
if (Sys.getenv("KSAA_UNWEIGHTED") == "1") stopifnot(round(PH, 2) == c(0.57, 0.50))   # check against the numbers in the text (unweighted version)
HD <- data.frame(y = c(10.5, 5.1), lab = c("All-cause mortality", "Cardiovascular mortality"),
                 p = sprintf("P for heterogeneity = %s", f2(PH, 2)))
pw <- fbase(LB, c(0.35, 11.2)) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.45, colour = NAVY) +
  geom_point(shape = SHP[[KS]], size = SZ[[KS]], colour = NAVY) +
  geom_text(aes(x = XT, label = txt), hjust = 0, size = FS / .pt, family = "Arial", colour = INK) +
  geom_text(data = HD, aes(x = XLO, y = y, label = lab), inherit.aes = FALSE, hjust = 0, vjust = 0.5,
            size = FS / .pt, family = "Arial", fontface = "bold", colour = INK) +
  geom_text(data = HD, aes(x = XT, y = y, label = p), inherit.aes = FALSE, hjust = 0, vjust = 0.5,
            size = (FS - 0.5) / .pt, family = "Arial", colour = INK) +
  scale_y_continuous(breaks = LB$y, labels = LB$stratum, expand = c(0, 0)) +
  labs(title = "B", subtitle = "KSAA within quartiles of Phenotypic Age acceleration",
       x = "HR per SD of KSAA (95% CI)", y = NULL) +
  theme(plot.subtitle = element_text(size = FS, colour = INK, margin = margin(b = 2)))

# ---- C-E: HRs of the cross-classified cells (reference Q1 x Q1; colored by log HR) ----
QL <- paste0("Q", 1:4)
JC$y <- factor(JC$y, QL); JC$x <- factor(JC$x, QL)
LIM <- log(5); BRK <- c(0.5, 1, 2, 4)   # covers the highest cell (cardiovascular Q4 x Q4, HR 4.49)
FW <- as.numeric(Sys.getenv("FIG3_W", "6.27")); FH <- as.numeric(Sys.getenv("FIG3_H", "5.5")); R3 <- as.numeric(Sys.getenv("FIG3_R3", "0.62")); HM_ASP <- as.numeric(Sys.getenv("FIG3_ASP", "0"))
stopifnot(max(abs(log(JC$hr))) < LIM)   # the heat maps share one color scale; no cell may exceed it
JC$dark <- abs(log(JC$hr)) > 0.8
HM_FS <- as.numeric(Sys.getenv("FIG3_HMFS", "4.2"))   # axis labels / axis titles / color-bar text of the heat maps (the rest of the figure stays at FS = 5.2)
HM_TY <- as.numeric(Sys.getenv("FIG3_HMTY", "0.97"))   # vertical position of panel letters (npc of the whole panel, 1 = top): smaller is closer to the heat map
HM_RM <- as.numeric(Sys.getenv("FIG3_HMRM", "12"))    # right margin of the heat maps (pt): the color bar goes here, and the value also sets the right edge of E
HM_LX <- as.numeric(Sys.getenv("FIG3_HMLX", "1.005"))  # color-bar position relative to the right edge of the tiles (1 = adjacent); too large overlaps the y-axis title of the next panel
HM_RME <- as.numeric(Sys.getenv("FIG3_HMRME", "0"))   # right margin of E (pt); 0 = tiles reach the cell edge and the color bar sits in the canvas margin
heat <- function(tag, letter, outcome, showleg = FALSE, rm = HM_RM) {   # the heat maps share one color scale, drawn once on the right map
  g <- JC[JC$panel == tag, ]
  ggplot(g, aes(x, y, fill = log(hr))) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = ifelse(ref, "1.00", f2(hr)), colour = dark), size = (FS - 2) / .pt, family = "Arial", fontface = "bold", vjust = -0.25) +
    geom_text(aes(label = ifelse(ref, "Reference", sprintf("(%s–%s)", f2(hr_lo), f2(hr_hi))), colour = dark), size = (FS - 2) / .pt, family = "Arial", vjust = 1.35) +
    scale_fill_gradient2(low = BLUE, mid = "white", high = ORANGE, midpoint = 0, limits = c(-LIM, LIM), breaks = log(BRK), labels = c("0.5", "1", "2", "4"), name = "HR",
                         guide = if (showleg) "colourbar" else "none") +
    scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "white"), guide = "none") +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(limits = QL, expand = c(0, 0)) +   # Q1 at the bottom, Q4 at the top
    annotate("text", x = 2.5, y = 4.5, label = outcome, vjust = -0.9, size = FS / .pt, family = "Arial", colour = INK) +   # outcome name centered on the tiles (subtitle is only a placeholder)
    coord_cartesian(clip = "off") +
    labs(tag = letter, subtitle = " ", x = "Phenotypic Age acceleration quartile", y = "KSAA quartile") + thk +
    theme(axis.line = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
          # panel letters via plot.tag (gtable clamps a negative plot.title margin); y position HM_TY
          plot.tag = element_text(size = 8, face = "bold", colour = INK, hjust = 0, vjust = 1),
          plot.tag.position = c(0.002, HM_TY), plot.tag.location = "plot",
          axis.text = element_text(colour = INK, size = HM_FS), axis.title = element_text(colour = INK, size = HM_FS),
          axis.title.x = element_text(margin = margin(t = 1)), axis.title.y = element_text(margin = margin(r = 1)),
          plot.subtitle = element_text(size = FS, margin = margin(b = 3)),
          # color bar drawn as an inside legend in the space right of the tiles, taking no layout space
          legend.position = "inside", legend.position.inside = c(HM_LX, 0.5), legend.justification.inside = c(0, 0.5), legend.direction = "vertical",
          plot.margin = margin(2, rm, 1, 2),   # room on the right for the color bar (otherwise clipped); E uses a smaller value so the color bar ends at the right end of the x axis of B
          legend.title = element_text(size = HM_FS, colour = INK, margin = margin(b = 2)), legend.text = element_text(size = HM_FS, colour = INK),
          legend.key.width = unit(3.5, "pt"), legend.key.height = unit(10, "pt"), legend.ticks.length = unit(1, "pt")) +
    # no fixed aspect ratio (HM_ASP = 0), so the panel edges can be aligned with B; tile shape follows the row height
    (if (HM_ASP > 0) theme(aspect.ratio = HM_ASP) else NULL)
}
pe <- heat("main", "C", "All-cause mortality", showleg = TRUE)
pf <- heat("main_cvd", "D", "Cardiovascular mortality", showleg = TRUE)
pg <- heat("main_cancer", "E", "Cancer mortality", showleg = TRUE, rm = HM_RME)

# shared legend taken from a dummy plot and placed in the top row (guides = "collect" would also collect the heat-map color bar)
leg_grob <- function(lv, sizes, w) {
  LG <- data.frame(series = factor(lv, lv), x = seq_along(lv), y = 1)
  pl <- ggplot(LG, aes(x, y, colour = series, shape = series)) + geom_point() + scale_colour_manual(values = PAL[lv], name = NULL) + scale_shape_manual(values = SHP[lv], name = NULL) + thk +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = sizes)))
  ragg::agg_png(tempfile(fileext = ".png"), width = w, height = 2, units = "in", res = 72); gl <- ggplotGrob(pl); invisible(dev.off())
  gl$grobs[[grep("guide-box-top", gl$layout$name)]]
}
NOLEG <- theme(legend.position = "none")
FH <- as.numeric(Sys.getenv("FIG3_H", "3.6"))   # two rows; without a fixed aspect ratio the row height is the height of the tiles
# pw (HR of KSAA within each PhenoAge quartile) is computed but not plotted; its heterogeneity P is cited in the text.
# A and B are wrapped with wrap_elements; otherwise patchwork aligns their axis areas and the axis of A detaches from its panel
toprow <- wrap_elements(full = pc + NOLEG) + wrap_elements(full = pd + NOLEG) + plot_layout(widths = c(1.15, 1.25))
BL <- as.numeric(Sys.getenv("FIG3_BL", "0.11")); BR <- as.numeric(Sys.getenv("FIG3_BR", "0"))   # side margins: align C with the y axis of A and E with the end of the x axis of B
botrow <- plot_spacer() + pe + pf + pg + plot_spacer() + plot_layout(widths = c(BL, 1, 1, 1, BR))   # three heat maps side by side (cancer cells have few events and wide CIs)
p <- wrap_plots(L = wrap_elements(full = leg_grob(c(PA, KS, PK), c(1.6, 1.5, 1.9), FW)), A = wrap_elements(full = toprow), B = wrap_elements(full = botrow),
                design = "L\nA\nB", heights = c(.06, .68, R3)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
ggsave(tempfile(fileext = ".png"), p, width = FW, height = FH, dpi = 300, bg = "white")   # warm-up render: the first render can mis-measure text sizes
save_fig(p, "Figure3", FW, FH)

cat("\nPanel values (n = 4,522; combined-score values are apparent, without optimism correction)\n")
for (i in seq_len(nrow(CI7))) with(CI7[i, ], cat(sprintf("  a %-15s C PhenoAge %.4f -> combined %.4f  ΔC %+.4f (%+.4f, %+.4f)  relative %+.2f%% (%+.2f%%, %+.2f%%)\n",
  as.character(o), C_pa, C_cb, dC_vs_pa, dC_vs_pa_lo, dC_vs_pa_hi, relC, relC_lo, relC_hi)))
for (i in seq_len(nrow(CI7))) with(CI7[i, ], cat(sprintf("  b %-15s R²_D PhenoAge %.3f -> combined %.3f  relative %+.1f%% (%+.1f%%, %+.1f%%)\n",
  as.character(o), R2_pa, R2_cb, relR2, relR2_lo, relR2_hi)))
for (i in which(IMP$metric == "per_sd")) with(IMP[i, ], cat(sprintf("  c %-15s gain %+.1f%% (%+.1f%%, %+.1f%%)\n", as.character(o), pct, pct_lo, pct_hi)))
for (i in seq_len(nrow(WD))) with(WD[i, ], cat(sprintf("  d %-15s KSAA %s | PhenoAge %s | Wald P %.2g\n", as.character(o),
  ci(hr_ksaa, lo_ksaa, hi_ksaa), ci(hr_pa, lo_pa, hi_pa), p_diff)))
qhet <- function(w) { b <- log(w$hr); se <- (log(w$hi) - log(w$lo)) / (2 * qnorm(.975)); wt <- 1 / se^2   # heterogeneity of the per-SD HR of KSAA across the four quartiles (Cochran Q, df = 3)
  Q <- sum(wt * (b - sum(wt * b) / sum(wt))^2); c(Q = Q, p = pchisq(Q, length(b) - 1, lower.tail = FALSE)) }
for (tg in c("main", "main_cvd")) { g <- JC[JC$panel == tg, ]; w <- JW[JW$panel == tg, ]
  h <- qhet(w); cat(sprintf("  %s heterogeneity across quartiles Cochran Q = %.2f, P = %.2f\n", if (tg == "main") "e" else "f", h[["Q"]], h[["p"]]))
  q4 <- g[g$y == "Q4" & g$x == "Q4", ]
  cat(sprintf("  %s %-9s Q1×Q1 %s, Q4×Q4 %s, reference Q1×Q1; %d–%d participants and %d–%d deaths per cell; KSAA per SD within quartiles: %s\n",
              if (tg == "main") "e" else "f", if (tg == "main") "all-cause" else "cardiovascular", f2(g$hr[g$y == "Q1" & g$x == "Q1"]),
              ci(q4$hr, q4$hr_lo, q4$hr_hi), min(g$n), max(g$n), min(g$deaths), max(g$deaths), paste(ci(w$hr, w$lo, w$hi), collapse = " / "))) }
