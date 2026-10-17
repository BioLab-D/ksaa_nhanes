suppressMessages(library(survival))
res <- here::here("results")
d <- read.csv(file.path(res,"ksaa_oof.csv"))
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); i <- match(d$SEQN, sl$SEQN)
d$t <- sl$PERMTH_EXM[i]; d$e <- sl$MORTSTAT[i]

cat("=========== 1. Follow-up ===========\n")
rk <- survfit(Surv(t, 1-e) ~ 1, d)
cat(sprintf("n %d; deaths %d (%.1f%%); person-years %.0f; median follow-up %.1f years (max %.1f); reverse-KM median %.1f years\n",
            nrow(d), sum(d$e), 100*mean(d$e), sum(d$t)/12, median(d$t)/12, max(d$t)/12, summary(rk)$table["median"]/12))

cat("\n=========== 2. Trivial baseline ===========\n")
bm <- sapply(0:4, function(k){ tr <- d$fold!=k; te <- d$fold==k; mean(abs(median(d$age[tr]) - d$age[te])) })
rm <- sapply(0:4, function(k){ te <- d$fold==k; mean(abs(d$pred_raw[te]-d$age[te])) })
cat(sprintf("  Training-fold median prediction MAE %.3f +/- %.3f; raw ensemble %.3f; relative reduction %.1f%%\n",
            mean(bm), sd(bm), mean(rm), 100*(1-mean(rm)/mean(bm))))
cat(sprintf("  Pooled OOF: raw r %.3f; corrected r %.3f\n", cor(d$pred_raw, d$age), cor(d$pred_corr, d$age)))
