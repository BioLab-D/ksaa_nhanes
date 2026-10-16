# 图的公共部分:路径、配色、主题、森林图、结果文件解析、保存
suppressMessages({library(ggplot2); library(patchwork); library(survival)})
RES <- here::here("results"); V3 <- file.path(RES, "analysis"); FIGD <- here::here("results", "figures")
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; SEQ_LO <- "#cde2fb"; SEQ_HI <- "#0d366b"
th <- theme_minimal(base_family = "Arial", base_size = 8) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
  axis.text = element_text(colour = INK2, size = 7), axis.title = element_text(colour = INK, size = 7.5), axis.line = element_line(colour = AXIS, linewidth = 0.3),
  plot.title = element_text(colour = INK, size = 8.5, face = "bold"), plot.margin = margin(4, 6, 2, 4), strip.text = element_text(colour = INK, size = 7.5, face = "bold"),
  legend.text = element_text(colour = INK, size = 7), legend.title = element_text(colour = INK, size = 7))
mi <- function(x) gsub("-", "\u2212", x, fixed = TRUE)
f2 <- function(x, d = 2) mi(sub("^-(0\\.0+)$", "\\1", formatC(round(x + sign(x) * 1e-9, d), format = "f", digits = d)))  # 常规四舍五入:1.035 → 1.04(formatC 按二进制会给 1.03)
ci <- function(h, l, u, d = 2) ifelse(is.na(h), "", sprintf("%s (%s%s%s)", f2(h, d), f2(l, d), ifelse(l < 0 | u < 0, " to ", "\u2013"), f2(u, d)))
cm <- function(x) formatC(round(x), format = "d", big.mark = ",")

# 森林图:F 需列 lab、hdr(逻辑)、est、lo、hi、txt;可选 series(着色)、filled。三栏(标签 | 森林 | 数值)各自成图再拼,避免文字压到相邻面板
forest <- function(F, title, xlab, xlim, breaks, ref = 1, pal = c(main = BLUE), widths = c(1.35, 1, 0.62), lab_size = 2.4) {
  F$y <- rev(seq_len(nrow(F))); if (is.null(F$series)) F$series <- names(pal)[1]; if (is.null(F$filled)) F$filled <- TRUE
  yl <- c(0.5, nrow(F) + 0.5); G <- F[!F$hdr, ]; many <- length(pal) > 1
  pl <- ggplot(F, aes(0, y, label = lab)) + geom_text(aes(fontface = ifelse(hdr, "bold", "plain"), x = ifelse(hdr, 0, 0.04)), hjust = 0, size = lab_size, family = "Arial", colour = INK) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = yl, expand = c(0, 0)) + labs(title = title) +
    theme_void(base_family = "Arial") + theme(plot.title = element_text(colour = INK, size = 8.5, face = "bold", margin = margin(b = if (many) 16 else 4)), plot.margin = margin(4, 0, 22, 4))
  pf <- ggplot(G) + geom_vline(xintercept = ref, linetype = 2, colour = MUTED, linewidth = 0.3) +
    geom_linerange(aes(xmin = pmax(lo, xlim[1]), xmax = pmin(hi, xlim[2]), y = y, colour = series), linewidth = 0.55) +
    geom_point(aes(est, y, colour = series), shape = ifelse(G$filled, 16, 21), fill = "white", size = 1.9, stroke = 0.6) +
    scale_colour_manual(values = pal, name = NULL) + scale_x_continuous(breaks = breaks, limits = xlim, expand = c(0, 0), labels = function(x) mi(format(x))) +
    scale_y_continuous(limits = yl, expand = c(0, 0), breaks = NULL) + labs(title = " ", x = xlab, y = NULL) + th +
    theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(), plot.margin = margin(4, 4, 2, 2),
          legend.position = if (many) "top" else "none", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"))
  pv <- ggplot(F, aes(0, y, label = txt)) + geom_text(hjust = 0, size = lab_size - 0.05, family = "Arial", colour = INK) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = yl, expand = c(0, 0)) + labs(title = " ") +
    theme_void(base_family = "Arial") + theme(plot.title = element_text(size = 8.5, margin = margin(b = if (many) 16 else 4)), plot.margin = margin(4, 4, 22, 2))
  pl + pf + pv + plot_layout(widths = widths)
}
hdrrow <- function(lab) data.frame(lab = lab, hdr = TRUE, est = NA, lo = NA, hi = NA, txt = "")
estrow <- function(lab, est, lo, hi, ...) data.frame(lab = lab, hdr = FALSE, est = est, lo = lo, hi = hi, txt = ci(est, lo, hi), ...)

# 解析结果文件里的 "HR x (lo–hi)",解析失败即报错
grab <- function(file, key, after = "HR ") {
  l <- grep(key, readLines(file.path(RES, file), encoding = "UTF-8"), value = TRUE, fixed = TRUE); if (length(l) < 1) stop("找不到: ", key)
  s <- sub(paste0("^.*?", after), "", l[1], perl = TRUE); x <- regmatches(s, regexec("([0-9.]+) \\(([0-9.]+)[-\u2013]([0-9.]+)\\)", s))[[1]]
  if (length(x) < 4) stop("解析失败: ", key); as.numeric(x[2:4])
}
save_fig <- function(p, name, w, h) {  # PNG 必出;PDF 尽力而为
  ggsave(file.path(FIGD, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  # cairo 加载失败时 ggsave 可能不报错也不写文件:先删旧 PDF,cairo 没写出就改用 quartz
  pdf_f <- file.path(FIGD, paste0(name, ".pdf")); if (file.exists(pdf_f)) file.remove(pdf_f)
  try_dev <- function(dev) tryCatch({ suppressWarnings(dev()); file.exists(pdf_f) && file.size(pdf_f) > 0 },
                                    error = function(e) { cat("  PDF(", conditionMessage(e), ")\n"); FALSE })
  ok <- try_dev(function() ggsave(pdf_f, p, width = w, height = h, device = cairo_pdf, bg = "white"))
  if (!ok && capabilities("aqua")) ok <- try_dev(function() ggsave(pdf_f, p, width = w, height = h, device = grDevices::quartz, type = "pdf", bg = "white"))
  cat("->", name, if (ok) ".png/.pdf\n" else ".png(PDF 未写出)\n")
}
