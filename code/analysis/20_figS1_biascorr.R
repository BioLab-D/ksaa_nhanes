suppressMessages({library(ggplot2); library(patchwork)})
res <- here::here("results"); fig <- here::here("results", "figures")
d <- read.csv(file.path(res,"ksaa_oof.csv"))
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN)
d$sex <- sl$sex[i]; d$phase <- sl$phase[i]; d$gap <- d$gap_corr
BLUE <- "#2a78d6"; RED <- "#e34948"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
th <- theme_minimal(base_family="Arial", base_size=7) +
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(colour=GRID, linewidth=0.25),
  axis.text=element_text(colour=INK2, size=6.5), axis.title=element_text(colour=INK, size=7), axis.line=element_line(colour=AXIS, linewidth=0.3),
  plot.title=element_text(colour=INK, size=7.5, hjust=0.5), plot.tag=element_text(face="bold", size=9), plot.margin=margin(2,4,2,2))
YM <- c(1.75, 1.4)  # left/right margins of y-axis tick labels (pt), set explicitly so the two rows align; ggplot2 4.0 needs axis.text.y.left
th <- th + theme(axis.text.y.left=element_text(margin=margin(l=YM[1], r=YM[2])))
HH <- c(0.8, 5)   # height ratio of marginal histogram : scatter plot
lim <- c(58, 92); YL <- c(-22, 27)  # extra headroom on the y axis for the top-right annotation (r, n) so it does not cover points
AGELAB <- c(60, 70, 80, "90+"); cm <- function(x) formatC(x, format="d", big.mark=",")
gapsc <- function(dd, y, ttl, corrected=FALSE){ r <- cor(dd[[y]], dd$age); f <- lm(dd[[y]]~dd$age)
  p <- ggplot(dd, aes(age, .data[[y]])) + geom_hline(yintercept=0, colour=MUTED, linewidth=0.3) +
    geom_point(colour=BLUE, alpha=0.18, size=0.5, stroke=0, shape=16) +
    geom_abline(intercept=coef(f)[1], slope=coef(f)[2], colour=BLUE, linewidth=0.9) +
    annotate("text", x=lim[2]-1, y=YL[2]-0.8, hjust=1, vjust=1, size=2.3, family="Arial", colour=INK, lineheight=1.05,
    label=sprintf("r = %s%.2f\nn = %s", ifelse(r<0,"−",""), abs(r), cm(nrow(dd)))) +
    scale_x_continuous(limits=lim, breaks=seq(60,90,10), labels=AGELAB, expand=c(0,0)) + scale_y_continuous(limits=YL, breaks=seq(-20,20,10)) +
    labs(title=ttl, x="Chronological age, years", y=if(corrected) "KSAA, years" else "Age gap, years") + th
  top <- ggplot(dd, aes(age)) + geom_histogram(binwidth=1, fill=BLUE, alpha=0.35, colour=NA, boundary=59.5) +
      scale_x_continuous(limits=lim, expand=c(0,0)) + theme_void()
  right <- ggplot(dd, aes(.data[[y]])) + geom_histogram(binwidth=1, fill=BLUE, alpha=0.35, colour=NA) +
      scale_x_continuous(limits=YL) + coord_flip() + theme_void()
  top + plot_spacer() + p + right + plot_layout(ncol=2, widths=c(5,1), heights=HH) }
errsc <- function(y, ttl, corrected=FALSE){ dd <- d; dd$ae <- abs(dd[[y]])
  p <- ggplot(dd, aes(age, ae)) +
    geom_point(colour=BLUE, alpha=0.18, size=0.5, stroke=0, shape=16) +
    geom_smooth(method="loess", formula=y~x, span=0.4, se=FALSE, colour=RED, linewidth=1) +
    scale_x_continuous(limits=lim, breaks=seq(60,90,10), labels=AGELAB, expand=c(0,0)) +
    scale_y_continuous(breaks=seq(0,20,5)) + coord_cartesian(ylim=c(0,20)) +
    labs(title=ttl, x="Chronological age, years", y="Absolute error, years") + th +
    theme(axis.text.y.left=element_text(margin=margin(l=YM[1]+2.17, r=YM[2])))  # pad by one hyphen width (2.17 pt) so the y axes of both rows align
  void <- ggplot() + theme_void()  # empty placeholder with the layout of gapsc() (plot_spacer gives unequal panel widths)
  void + plot_spacer() + p + void + plot_layout(ncol=2, widths=c(5,1), heights=HH) }
# A
A1 <- wrap_elements(full=gapsc(d, "gap_raw", "Before bias correction"))
A2 <- wrap_elements(full=gapsc(d, "gap", "After bias correction", TRUE))
A3 <- wrap_elements(full=errsc("gap_raw", "Before bias correction"))
A4 <- wrap_elements(full=errsc("gap", "After bias correction", TRUE))
# B: phase
B <- lapply(1:2, function(p){ dd <- d[d$phase==p,]
  list(wrap_elements(full=gapsc(dd,"gap_raw",sprintf("Phase %d, before",p))),
      wrap_elements(full=gapsc(dd,"gap",sprintf("Phase %d, after",p),TRUE))) })
ttl <- function(x) ggplot() + annotate("text", x=0, y=0.5, label=x, hjust=0, size=3.1, fontface="bold",
    family="Arial", colour=INK) + coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off") + theme_void()
TH <- grid::unit(c(5, 1, 1), c("mm", "null", "null"))  # title row fixed at 5 mm so text is not clipped
left <- ttl("A. Bias correction, study population") / (A1 | A2) / (A3 | A4) + plot_layout(heights=TH)
right <- ttl("B. By survey phase") / (B[[1]][[1]] | B[[1]][[2]]) / (B[[2]][[1]] | B[[2]][[2]]) + plot_layout(heights=TH)
div <- ggplot() + geom_vline(xintercept=0, linetype="longdash", linewidth=0.9, colour=INK) + theme_void()
p <- wrap_elements(full=left) + div + wrap_elements(full=right) + plot_layout(widths=c(0.494, 0.012, 0.494))  # A and B equal width
# total figure height can be overridden by FIGS1_H (in)
FH <- as.numeric(Sys.getenv("FIGS1_H", "4.0"))
ggsave(file.path(fig,"FigureS1.png"), p, width=7.2, height=FH, dpi=300, bg="white")
pdf_f <- file.path(fig,"FigureS1.pdf"); if (file.exists(pdf_f)) file.remove(pdf_f)
tryCatch(suppressWarnings(ggsave(pdf_f, p, width=7.2, height=FH, device=cairo_pdf, bg="white")), error=function(e) NULL)
if (!file.exists(pdf_f) && capabilities("aqua")) tryCatch(ggsave(pdf_f, p, width=7.2, height=FH,
    device=grDevices::quartz, type="pdf", bg="white"), error=function(e) NULL)
cat(if (file.exists(pdf_f)) "  PDF written\n" else "  PDF not written\n")
# console: bias-correction formulas (overall, by sex) and correlations
f0 <- lm(pred_raw~age, d); fm <- lm(pred_raw~age, d[d$sex==1,]); ff <- lm(pred_raw~age, d[d$sex==2,])
cat(sprintf("gap~age r before correction %.3f, after %+.4f; E[raw|age] = %.2f + %.3fxage; men %.2f + %.3fxage; women %.2f + %.3fxage; corrected gap SD %.2f\n",
  cor(d$gap_raw,d$age), cor(d$gap,d$age), coef(f0)[1],coef(f0)[2], coef(fm)[1],coef(fm)[2], coef(ff)[1],coef(ff)[2], sd(d$gap)))
for(p_ in 1:2){ dd <- d[d$phase==p_,]
  cat(sprintf("phase %d: r before correction %.3f, after %+.3f; corrected MAE %.2f\n", p_, cor(dd$gap_raw,dd$age),
      cor(dd$gap,dd$age), mean(abs(dd$gap)))) }
cat("-> FigureS1.png/.pdf\n")
