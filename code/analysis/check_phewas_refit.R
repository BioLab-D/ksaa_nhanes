suppressMessages({ library(survey); library(MASS) }); options(survey.lonely.psu = "adjust")
d <- read.csv(here::here("data", "derived", "analysis_dataset.csv"), na.strings = "")
Y <- read.csv(here::here("data", "derived", "phewas", "phenotype_values.csv"), check.names = FALSE, na.strings = "")
W <- read.csv(here::here("data", "derived", "phewas", "survey_design.csv"))
R <- read.csv(here::here("data", "derived", "phewas", "phewas_results.csv"))
stopifnot(identical(d$SEQN, Y$SEQN), identical(d$SEQN, W$SEQN))
d$z <- d$gap / sd(d$gap)                                           # KSAA per SD (full-cohort SD)
for (v in c("sex", "race", "edu3")) d[[v]] <- factor(d[[v]])
COV <- c("age", "sex", "race", "edu3", "pir")
d <- cbind(d, W[, setdiff(names(W), c("SEQN", "phase", "WTPFEX6", "WTPFEX2", "SDPPSU6", "SDPSTRA6", "SDPPSU2", "SDPSTRA2"))])
des_all <- svydesign(id = ~SDPPSU6, strata = ~SDPSTRA6, weights = ~WTPFEX6, nest = TRUE, data = d)
des_fast <- svydesign(id = ~SDPPSU6, strata = ~SDPSTRA6, weights = ~WTPFSD6, nest = TRUE, data = d[d$WTPFSD6 > 0, ])

refit <- function(v) {
  r <- R[R$variable == v, ]; fast <- grepl("^WTPFSD6", r$weight); des <- if (fast) des_fast else des_all
  y <- Y[[v]][if (fast) d$WTPFSD6 > 0 else TRUE]; des$variables$y <- y
  ds <- subset(des, !is.na(y) & complete.cases(des$variables[, COV])); df <- degf(ds)
  f <- y ~ z + age + sex + race + edu3 + pir
  if (r$model_used == "gaussian") {
    s <- summary(svyglm(f, ds))$coefficients; b <- s["z", 1]; se <- s["z", 2]
  } else if (r$model_used == "binomial") {
    s <- summary(svyglm(f, ds, family = quasibinomial()))$coefficients; b <- s["z", 1]; se <- s["z", 2]
  } else if (r$model_used == "ordinal") {
    ds$variables$y <- factor(ds$variables$y, ordered = TRUE); m <- svyolr(f, ds)
    b <- coef(m)[["z"]]; se <- sqrt(vcov(m)["z", "z"])
  } else return(NULL)
  est <- if (r$model_used == "gaussian") b else exp(b); p <- 2 * pt(-abs(b / se), df)
  data.frame(variable = v, model = r$model_used, n = nrow(ds), n_reported = r$n, est = est, est_reported = r$est,
             p = p, p_reported = r$p, df = df)
}
# one gaussian, one binomial, one ordinal and one fasting-subsample phenotype from the main scan, plus the clump leads of the top 5 signals
pick <- unique(c("BDPTRBMD", "HAH9", "PEP13A", "G1P", head(unique(na.omit(R$clump_lead[R$set == "main"])), 5)))
out <- do.call(rbind, lapply(pick, refit)); print(out, digits = 4, row.names = FALSE)
stopifnot(all(out$n == out$n_reported),
          max(abs(out$est - out$est_reported) / abs(out$est_reported)) < 1e-6,
          max(abs(log(out$p) - log(out$p_reported))) < 1e-6)
cat("All refitted estimates and P values match phewas_results.csv.\n")
