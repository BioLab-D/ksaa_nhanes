# Export per-participant data for the interpretability scripts (gap = out-of-fold KSAA; pred_raw = uncorrected ensemble prediction).
# Output: results/analysis/saliency/ksaa_for_saliency.csv
# Run after run_all.sh: Rscript code/interpretability/saliency_input.R
d <- readRDS(here::here("results", "analysis", "analysis_data.rds"))
out <- here::here("results", "analysis", "saliency"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("SEQN", "age", "sex", "gap", "pred_raw", "fold", "race"); stopifnot(all(k %in% names(d)), nrow(d) == 4875)
write.csv(d[, k], file.path(out, "ksaa_for_saliency.csv"), row.names = FALSE); cat(nrow(d), "rows\n")
