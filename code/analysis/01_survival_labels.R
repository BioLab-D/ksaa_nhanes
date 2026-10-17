# Demographics, covariates and linked mortality of the imaging cohort from the NHANES III DuckDB
# Input: data/nhanes3/nhanes3_knee_cohort.duckdb; output: results/work/surv_labels.csv
suppressMessages(library(duckdb))
sp <- here::here("results", "work")
con <- dbConnect(duckdb::duckdb(),here::here("data", "nhanes3", "nhanes3_knee_cohort.duckdb"),read_only=TRUE)
dbWriteTable(con,"img",data.frame(SEQN=as.numeric(readLines(
  file.path(sp,"cdc_knee_seqn.txt"),warn=FALSE))),temporary=TRUE)
d <- dbGetQuery(con,"
 SELECT i.SEQN, a.SDPPHASE AS phase, a.HSAGEIR AS age, a.HSSEX AS sex,
        a.DMARETHN AS race, a.DMPPIR AS pir, a.HFA8R AS edu,
        e.BMPBMI AS bmi, e.BDPFNBMD AS fnbmd,
        a.WTPFEX2, a.WTPFEX6, a.SDPPSU2, a.SDPSTRA2, a.SDPPSU6, a.SDPSTRA6,
        m.MORTSTAT, m.PERMTH_INT, m.PERMTH_EXM, m.UCOD_LEADING   -- PERMTH_EXM: follow-up from the examination
 FROM img i JOIN nhanes3_adult a USING(SEQN)
 LEFT JOIN nhanes3_exam e USING(SEQN)
 LEFT JOIN nhanes3_mortality m USING(SEQN)")
dbDisconnect(con,shutdown=TRUE)

# NHANES III missing codes (all 8s) -> NA
d$bmi[d$bmi   >= 100] <- NA
d$fnbmd[d$fnbmd >= 3] <- NA
d$pir[d$pir     >= 50] <- NA
d$edu[d$edu     >= 90] <- NA

cat("Raw:", nrow(d), "\n")
d <- subset(d, !is.na(MORTSTAT) & !is.na(PERMTH_INT) & PERMTH_INT > 0)
cat("With valid follow-up:", nrow(d), "  deaths", sum(d$MORTSTAT==1),
    sprintf("(%.1f%%)", 100*mean(d$MORTSTAT==1)), "\n")
cat("Median follow-up", round(median(d$PERMTH_INT)/12,1), "years\n")
print(table(phase=d$phase, dead=d$MORTSTAT))
write.csv(d, file.path(sp,"surv_labels.csv"), row.names=FALSE)
cat("-> surv_labels.csv\n")
