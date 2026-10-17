# Shared figure code: paths, colors, theme, number formatting, saving
suppressMessages({library(ggplot2); library(patchwork); library(survival)})
RES <- here::here("results"); ANA <- file.path(RES, "analysis"); FIGD <- here::here("results", "figures")
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
th <- theme_minimal(base_family = "Arial", base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
  axis.text = element_text(colour = INK2, size = 7), axis.title = element_text(colour = INK, size = 7.5),
  axis.line = element_line(colour = AXIS, linewidth = 0.3),
  plot.title = element_text(colour = INK, size = 8.5, face = "bold"), plot.margin = margin(4, 6, 2, 4),
  strip.text = element_text(colour = INK, size = 7.5, face = "bold"),
  legend.text = element_text(colour = INK, size = 7), legend.title = element_text(colour = INK, size = 7))
mi <- function(x) gsub("-", "\u2212", x, fixed = TRUE)
# round half up: 1.035 -> 1.04 (formatC rounds the binary value and gives 1.03)
f2 <- function(x, d = 2) mi(sub("^-(0\\.0+)$", "\\1", formatC(round(x + sign(x) * 1e-9, d), format = "f", digits = d)))
ci <- function(h, l, u, d = 2) ifelse(is.na(h), "", sprintf("%s (%s%s%s)", f2(h, d), f2(l, d), ifelse(l < 0 | u < 0, " to ", "\u2013"), f2(u, d)))
cm <- function(x) formatC(round(x), format = "d", big.mark = ",")

save_fig <- function(p, name, w, h) {  # PNG always; PDF best effort
  ggsave(file.path(FIGD, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  # if cairo fails to load, ggsave may silently write nothing: remove the old PDF first and fall back to quartz
  pdf_f <- file.path(FIGD, paste0(name, ".pdf")); if (file.exists(pdf_f)) file.remove(pdf_f)
  try_dev <- function(dev) tryCatch({ suppressWarnings(dev()); file.exists(pdf_f) && file.size(pdf_f) > 0 },
                                    error = function(e) { cat("  PDF(", conditionMessage(e), ")\n"); FALSE })
  ok <- try_dev(function() ggsave(pdf_f, p, width = w, height = h, device = cairo_pdf, bg = "white"))
  if (!ok && capabilities("aqua")) ok <- try_dev(function() ggsave(pdf_f, p, width = w, height = h,
      device = grDevices::quartz, type = "pdf", bg = "white"))
  cat("->", name, if (ok) ".png/.pdf\n" else ".png (PDF not written)\n")
}
