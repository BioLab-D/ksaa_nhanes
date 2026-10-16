# 解释性分析 · 第 0 步:导出 saliency*.py 需要的逐人数据(gap = 交叉拟合 OOF 的 KSAA;pred_raw = 未校正 OOF 集成预测)。
# 跑法(仓库根目录下):Rscript code/interpretability/saliency_input.R
d <- readRDS(here::here("results", "analysis", "analysis_data.rds"))
out <- here::here("results", "analysis", "saliency"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
k <- c("SEQN", "age", "sex", "gap", "pred_raw", "fold", "race"); stopifnot(all(k %in% names(d)), nrow(d) == 4875)
write.csv(d[, k], file.path(out, "ksaa_for_saliency.csv"), row.names = FALSE); cat(nrow(d), "rows\n")
