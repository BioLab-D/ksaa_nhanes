d <- readRDS(here::here("results", "analysis", "analysis_data.rds"))
out <- here::here("results", "analysis", "saliency"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("SEQN", "age", "sex", "gap", "pred_raw", "fold", "race"); stopifnot(all(k %in% names(d)), nrow(d) == 4875)
write.csv(d[, k], file.path(out, "ksaa_for_saliency.csv"), row.names = FALSE); cat(nrow(d), "rows\n")
