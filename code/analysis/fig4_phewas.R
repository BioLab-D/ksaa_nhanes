# Figure 4B (forest plot) and the data for 4A; fig4_phewas_clarite.py draws 4A with CLARITE and combines the panels.
#   B: top NTOP independent signals (FDR < 5% hits clumped at |r| > 0.7, ordered by P), per SD of KSAA; binary/ordinal
#   phenotypes as OR (top block), continuous as SD difference (bottom block); numbers = rank, colors = health class.
# Input: results/analysis/phewas/phewas_results.csv
# Output: results/figures/Figure4A.csv, Figure4B.{png,pdf}
source(here::here("code", "analysis", "fig_common.R"))
source(here::here("code", "analysis", "phewas_labels.R"))
# full figure is 183 mm wide, font size 5.2 pt (as in Figures 2 and 3); B is 3.15 in wide
TS <- 5.2 / .pt
forest1 <- function(D, xl, fr, log = FALSE, breaks, blabels = breaks, ref, xlab, head) {
  tf <- if (log) log10 else identity; itf <- if (log) function(z) 10^z else identity
  L0 <- tf(xl[1]); L1 <- tf(xl[2]); W <- (L1 - L0) / fr[2]
  XMIN <- L0 - fr[1] * W; XT <- L1 + 0.02 * W; XMAX <- L1 + fr[3] * W
  srt <- sort(D$row[!D$hdr], decreasing = TRUE); blk <- split(srt, cumsum(c(1, diff(srt) != -1)))
  BL <- do.call(rbind, lapply(blk, function(r) data.frame(y0 = min(r) - 0.45, y1 = max(r) + 0.45)))
  oth <- setdiff(breaks, ref)
  VL <- merge(data.frame(x = c(ref, oth), lt = c("22", rep("solid", length(oth))), col = c(MUTED, rep(GRID, length(oth)))), BL)
  G <- D[!D$hdr & !is.na(D$est), ]
  ggplot() +
    geom_segment(data = VL, aes(x = x, xend = x, y = y0, yend = y1, linetype = lt, colour = I(col)), linewidth = 0.3) + scale_linetype_identity() +
    geom_segment(aes(x = itf(L0), xend = itf(L1), y = 0.4, yend = 0.4), colour = AXIS, linewidth = 0.3) +
    geom_linerange(data = G, aes(xmin = pmax(lo, xl[1]), xmax = pmin(hi, xl[2]), y = row, colour = I(col)), linewidth = 0.5) +
    geom_point(data = G, aes(est, row, colour = I(col)), shape = 16, size = 1.6) +
    geom_text(data = D, aes(x = itf(XMIN + ifelse(hdr, 0, 0.015 * W)), y = row, label = lab, fontface = ifelse(hdr, "bold", "plain")), hjust = 0, size = TS, family = "Arial", colour = INK) +
    geom_text(data = D, aes(x = itf(XT), y = row, label = txt), hjust = 0, size = TS - 0.1, family = "Arial", colour = INK) +
    annotate("text", x = itf(XT), y = max(D$row) + 0.95, label = head, hjust = 0, size = TS, fontface = "bold", family = "Arial", colour = INK) +
    (if (log) scale_x_log10(limits = itf(c(XMIN, XMAX)), breaks = breaks, labels = blabels, expand = c(0, 0))
     else scale_x_continuous(limits = c(XMIN, XMAX), breaks = breaks, labels = blabels, expand = c(0, 0))) +
    coord_cartesian(clip = "off") + scale_y_continuous(limits = c(0.4, max(D$row) + 1.4), expand = c(0, 0), breaks = NULL) +
    labs(x = xlab, y = NULL) + th +
    theme(panel.grid.major = element_blank(), axis.line = element_blank(), axis.ticks.x = element_line(colour = AXIS, linewidth = 0.3),
          axis.title.x = element_text(hjust = ((L0 + L1) / 2 - XMIN) / (XMAX - XMIN), size = 5.2), axis.text.x = element_text(size = 5.2), plot.margin = margin(2, 30, 2, 4))
}
R <- read.csv(file.path(V3, "phewas", "phewas_results.csv")); R <- R[R$set == "main" & !is.na(R$p), ]
stopifnot(all(R$health_class %in% CLS_ORDER))
L <- R[!is.na(R$clump_lead) & R$clump_lead == R$variable, ]; L <- L[order(L$p), ]
NTOP <- as.numeric(Sys.getenv("PHEWAS_NTOP", "25"))
TOP <- head(L$variable, NTOP); rk <- setNames(seq_along(TOP), TOP)   # labels in A show the same ranks
SHORT <- c(BDPTRBMD = "Trochanter BMD", PEP1 = "Abnormal locomotion (examination)", HAY11D = "Walks slowly (observed)", PFPNSTPA = "Steps to walk 8 ft", HAH9 = "Chair-rise difficulty",
           HAH3 = "Stooping difficulty", HAG34 = "Chronic limp", PEP13E5A = "Heavy-housework difficulty (physician)", PFPRFLEX = "Limited hip and knee flexion", PEP13A = "Poorer physician-rated health", HAH2 = "10-step walk difficulty", HAG22 = "Knee pain",
           APPSI = "Alkaline phosphatase", HAH1 = "\u00bc-mile walk difficulty", HAH15 = "Uses walking aid (self-reported)", HAC1A = "Arthritis", CRP = "C-reactive protein",
           # the three walking-aid items name their source so they do not look like duplicates
           HAY11C = "Uses walking aid (observed)", FPP1090 = "Hard exudates (retina)", PFPDEVIC = "Type of walking aid", HAH5 = "Household-chore difficulty",
           HAH6 = "Meal-preparation difficulty", HAH8 = "Room-to-room walk difficulty", HAG5A = "Self-reported hip fracture", PEP10B4 = "Knee pain on passive motion", HAH4 = "Lifting 10 lb difficulty",
           PFP_steps_8ft_mean = "Steps to walk 8 ft", PEP13E1A = "\u00bc-mile walk difficulty (physician)", PEP13E2A = "100-yard run difficulty (physician)")
short <- function(v, lab) ifelse(v %in% names(SHORT), SHORT[v], nice(v, lab))
# ---- data for A (for CLARITE) ----
A <- data.frame(variable = R$variable, label = ifelse(R$variable %in% TOP, sprintf("%d  %s", rk[R$variable], short(R$variable, R$phenotype)), short(R$variable, R$phenotype)),
                category = CLS_NAME[R$health_class], colour = CLS_COL[R$health_class], p = R$p, fdr = R$fdr, replicated = R$replicated, is_top = R$variable %in% TOP)
write.csv(A, file.path(FIGD, "Figure4A.csv"), row.names = FALSE)
# ---- B ----
mrow <- function(v) { r <- R[R$variable == v, ]; data.frame(lab = sprintf("%d  %s", rk[[v]], short(v, r$phenotype)), hdr = FALSE, est = r$est, lo = r$lo, hi = r$hi, txt = ci(r$est, r$lo, r$hi), col = CLS_COL[[r$health_class]]) }
hrow <- function(lab) data.frame(lab = lab, hdr = TRUE, est = NA, lo = NA, hi = NA, txt = "", col = INK)
VOR <- TOP[R$model[match(TOP, R$variable)] %in% c("binomial", "ordinal")]; VSD <- TOP[R$model[match(TOP, R$variable)] == "gaussian"]; stopifnot(length(VOR) + length(VSD) == NTOP)
# within each block of B, order by health domain (= point color), then by P
ordcls <- function(vs) { i <- match(vs, R$variable)   # alphabetical, as in the legend of A (CLARITE's order), not CLS_ORDER
  vs[order(CLS_NAME[R$health_class[i]], R$p[i], method = "radix")] }
VOR <- ordcls(VOR); VSD <- ordcls(VSD)
B1 <- rbind(hrow("Binary and ordinal phenotypes"), do.call(rbind, lapply(VOR, mrow))); B1$row <- rev(seq_len(nrow(B1)))
B2 <- rbind(hrow("Continuous phenotypes"), do.call(rbind, lapply(VSD, mrow))); B2$row <- rev(seq_len(nrow(B2)))
XL1 <- c(min(0.8, floor(min(B1$lo, na.rm = TRUE) * 10) / 10), max(1.8, ceiling(max(B1$hi, na.rm = TRUE) * 10) / 10)); BR1 <- c(0.5, 0.8, 1, 1.5, 2, 3); BR1 <- BR1[BR1 >= XL1[1] & BR1 <= XL1[2]]
XL2 <- max(0.25, ceiling(max(abs(c(B2$lo, B2$hi)), na.rm = TRUE) * 10) / 10); XL2 <- c(-XL2, XL2); BR2 <- pretty(XL2, 4); BR2 <- BR2[BR2 > XL2[1] & BR2 < XL2[2]]
pb1 <- forest1(B1, xl = XL1, fr = c(0.57, 0.21, 0.22), log = TRUE, breaks = BR1, blabels = as.character(BR1), ref = 1, xlab = "Odds ratio per SD of KSAA", head = "OR (95% CI)")
pb2 <- forest1(B2, xl = XL2, fr = c(0.57, 0.21, 0.22), breaks = BR2, blabels = mi(as.character(BR2)), ref = 0, xlab = "SD difference per SD of KSAA", head = "Difference, SD (95% CI)")
pb <- pb1 / pb2 + plot_layout(heights = c(nrow(B1) + 2.6, nrow(B2) + 2.6))
PF <- as.numeric(Sys.getenv("PHEWAS_PITCH", "0.78"))   # row-spacing compression factor
BH <- PF * 4.2 * (nrow(B1) + nrow(B2) + 5.2) / (17 + 5.2)
save_fig(pb, "Figure4B", 3.15, BH)
cat(sprintf("A data: %d phenotypes; B: top %d (OR %d / SD %d), height %.2f in\n", nrow(A), NTOP, length(VOR), length(VSD), BH))
for (i in seq_len(NTOP)) cat(sprintf("  %2d %-9s %-40s %s\n", i, L$variable[i], short(L$variable[i], L$phenotype[i]), ci(L$est[i], L$lo[i], L$hi[i])))
