suppressMessages(library(duckdb))
res <- here::here("results"); ew <- file.path(res, "analysis", "phewas"); out <- here::here("data", "derived")
write_mode <- Sys.getenv("WRITE_DERIVED") == "1"
bad <- character()

# ---------- build each derived table from results/ ----------
d <- readRDS(file.path(res, "analysis", "analysis_data.rds"))
fac2chr <- function(x) { for (v in names(x)) if (is.factor(x[[v]])) x[[v]] <- as.character(x[[v]]); x }
e <- read.csv(file.path(res, "ksaa_oof.csv")); names(e)[names(e) == "gap_corr"] <- "ksaa"
A5 <- c(resnet18 = "ResNet-18", resnet50 = "ResNet-50", efficientnetb0 = "EfficientNet-B0", convnexttiny = "ConvNeXt-Tiny",
        densenet121 = "DenseNet-121")
W <- do.call(rbind, lapply(names(A5), function(a) {
  f <- read.csv(file.path(res, "predictions", sprintf("cv_%s_folds.csv", a)))
  data.frame(fold = f$fold, backbone = A5[[a]], inner_val_mae = f$inner_val_mae)
}))
W$weight <- ave(1 / W$inner_val_mae, W$fold, FUN = function(x) x / sum(x))
W <- W[order(W$fold, match(W$backbone, A5)), ]; rownames(W) <- NULL
B <- do.call(rbind, lapply(0:4, function(k) {
  f <- lm(I(pred_raw - age) ~ age, e[e$fold != k, ])
  data.frame(fold = k, intercept = coef(f)[[1]], slope = coef(f)[[2]])
}))
# the fold weights and bias-correction coefficients must reproduce pred_raw and ksaa
M <- sapply(names(A5), function(a) {
  x <- read.csv(file.path(res, "predictions", sprintf("cv_%s_oof.csv", a)))
  x$pred_age[match(e$SEQN, x$SEQN)]
})
raw <- ks <- rep(NA_real_, nrow(e))
for (k in 0:4) {
  i <- e$fold == k
  raw[i] <- M[i, ] %*% W$weight[W$fold == k]
  ks[i] <- raw[i] - e$age[i] - (B$intercept[k + 1] + B$slope[k + 1] * e$age[i])
}
if (max(abs(raw - e$pred_raw)) > 1e-6 || max(abs(ks - e$ksaa)) > 1e-6)
  bad <- c(bad, "fold weights / bias correction do not reproduce pred_raw / ksaa")
# calibration for the full-data checkpoints in models/: weights averaged over folds, bias correction fitted on all
# out-of-fold predictions (not on the full-data models' in-sample predictions)
fb <- coef(lm(pred_raw ~ age, e))
CAL <- rbind(data.frame(item = "ensemble_weight", backbone = unname(A5), checkpoint = sprintf("full_%s.pt", names(A5)),
                        value = sapply(unname(A5), function(b) mean(W$weight[W$backbone == b]))),
             data.frame(item = c("expected_pred_intercept", "expected_pred_slope", "ksaa_sd"), backbone = NA, checkpoint = NA,
                        value = c(fb[[1]], fb[[2]], sd(e$ksaa))))
if (abs(sum(CAL$value[1:5]) - 1) > 1e-12) bad <- c(bad, "new-data ensemble weights do not sum to 1")
Y <- readRDS(file.path(ew, "phewas_y.rds"))
con <- dbConnect(duckdb(), here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"), read_only = TRUE)
dbWriteTable(con, "coh", data.frame(SEQN = d$SEQN), temporary = TRUE)
S <- dbGetQuery(con, "select SEQN, WTPFEX6, SDPPSU6, SDPSTRA6, WTPFEX1, SDPPSU1, SDPSTRA1, WTPFEX2, SDPPSU2, SDPSTRA2,
                      WTPFSD6, WTPFSD1, WTPFSD2 from nhanes3_adult join coh using (SEQN)")
dbDisconnect(con, shutdown = TRUE)
S <- S[match(d$SEQN, S$SEQN), ]; for (v in setdiff(names(S), "SEQN")) S[[v]][is.na(S[[v]])] <- 0
S <- data.frame(SEQN = S$SEQN, phase = d$phase, S[, -1])

TAB <- list("analysis_dataset.csv" = fac2chr(d), "ensemble_predictions.csv" = e, "ensemble_weights_by_fold.csv" = W,
            "bias_correction_by_fold.csv" = B, "new_data_calibration.csv" = CAL,
            "phewas/variable_audit_table.csv" = read.csv(file.path(ew, "phewas_audit_table.csv"), check.names = FALSE),
            "phewas/phewas_results.csv" = read.csv(file.path(ew, "phewas_results.csv"), check.names = FALSE),
            "phewas/phenotype_values.csv" = data.frame(SEQN = d$SEQN, Y, check.names = FALSE), "phewas/survey_design.csv" = S)

# ---------- write or compare ----------
same <- function(x, y, tag) {
  if (!identical(names(x), names(y)) || nrow(x) != nrow(y))
    return(sprintf("%s: columns or rows differ (%d x %d vs %d x %d)", tag, nrow(x), ncol(x), nrow(y), ncol(y)))
  for (v in names(x)) {
    a <- x[[v]]; b <- y[[v]]
    if (!identical(is.na(a), is.na(b))) return(sprintf("%s: missing pattern differs in %s", tag, v))
    if (is.numeric(a) && is.numeric(b)) {
      dm <- suppressWarnings(max(abs(a - b), na.rm = TRUE))
      # 1e-6 relative: ordinal fits (svyolr) differ at ~1e-7 between runs from optimizer noise
      if (is.finite(dm) && dm > 1e-6 * max(1, max(abs(b), na.rm = TRUE))) return(sprintf("%s: %s max diff %.3g", tag, v, dm))
    } else if (!identical(as.character(a), as.character(b))) return(sprintf("%s: %s differs", tag, v))
  }
  NULL
}
for (f in names(TAB)) {
  p <- file.path(out, f)
  if (write_mode) {
    dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
    write.csv(TAB[[f]], p, row.names = FALSE, na = ""); next
  }
  if (!file.exists(p)) { bad <- c(bad, paste("missing", f)); next }
  r <- read.csv(p, check.names = FALSE, na.strings = c("", "NA"))   # PheWAS tables are copied as is, with NA for missing
  x <- TAB[[f]]   # empty strings are read back as missing
  for (v in names(x)) if (is.character(x[[v]])) x[[v]][!is.na(x[[v]]) & x[[v]] == ""] <- NA
  m <- same(r, x, f); if (!is.null(m)) bad <- c(bad, m)
}
sf <- file.path(out, "phewas", "phewas_summary.txt")
if (write_mode) file.copy(file.path(ew, "phewas_summary.txt"), sf, overwrite = TRUE) else
  if (!identical(readLines(sf), readLines(file.path(ew, "phewas_summary.txt")))) bad <- c(bad, "phewas/phewas_summary.txt differs")

if (write_mode) { cat("data/derived written from results/\n"); quit(save = "no") }
cat(sprintf("compared %d derived tables and the PheWAS summary with the regenerated results\n", length(TAB)))
if (length(bad)) { cat("MISMATCHES:\n", paste(" ", bad, collapse = "\n"), "\n"); quit(status = 1) }
cat("ALL IDENTICAL: the regenerated results match data/derived/\n")
