# Fracture, fall and physical function variables for Table 1, Table S10 (weighted.R / weighted_common.R) and phewas.R
#  Definitions: fracture HAG5A–C (8/9 missing), low-energy HAG10A–F, fall HAG14 (88/99 already removed in analysis_data), activity difficulty HAH ≥2, tandem stand < 10 s; counts checked by assertions
# Input results/analysis/analysis_data.rds, results/work/{fracture,physfun}.csv; output results/analysis/fxf_vars.csv
# Run by run_all.sh
res <- here::here("results"); out <- file.path(res, "analysis"); sp <- file.path(res, "work")
d <- readRDS(file.path(out, "analysis_data.rds"))
FR <- read.csv(file.path(sp, "fracture.csv")); PF <- read.csv(file.path(sp, "physfun.csv")); stopifnot(!anyDuplicated(FR$SEQN), !anyDuplicated(PF$SEQN))
g <- function(D, v) suppressWarnings(as.numeric(D[[v]][match(d$SEQN, D$SEQN)]))
yn <- function(x) ifelse(x %in% c(8, 9), NA, as.integer(x == 1))
same <- function(a, b) all((!is.na(a) & !is.na(b) & a == b) | (is.na(a) & is.na(b)))
hip <- yn(g(FR, "HAG5A")); d$fx_wri <- yn(g(FR, "HAG5B")); d$fx_spi <- yn(g(FR, "HAG5C"))
stopifnot(same(hip, d$fx_hip))
fx <- cbind(hip, d$fx_wri, d$fx_spi); any_ <- as.integer(rowSums(fx, na.rm = TRUE) > 0); any_[rowSums(!is.na(fx)) == 0] <- NA; stopifnot(same(any_, d$fx_any))
# Low-energy fracture: HAG10A–F give the injury mechanism of each fracture (1 = fall from standing height or less), asked only for fractures at age 50 or later (HAG9); 8/9 = don't know/refused.
#   Any site 1 -> positive; otherwise mechanism 8/9, or fracture after 50 with missing mechanism -> undeterminable (NA); else -> negative.
#   HAG10D–F are BOOLEAN in the database; as numeric they equal the original code 1.
num <- function(x) if (is.logical(x)) as.numeric(x) else suppressWarnings(as.numeric(x))
idx <- match(d$SEQN, FR$SEQN)
MEC <- sapply(paste0("HAG10", LETTERS[1:6]), function(v) num(FR[[v]])[idx])
AGE <- sapply(paste0("HAG9", LETTERS[1:6], "R"), function(v) { a <- num(FR[[v]])[idx]; a[a >= 500] <- NA; a })   # special codes such as 555 / 888 are not ages
n_low <- rowSums(MEC == 1 & !is.na(MEC))
n_unk <- rowSums((MEC %in% c(8, 9) & !is.na(MEC)) | (is.na(MEC) & !is.na(AGE) & AGE >= 50))
d$fx_frag <- ifelse(is.na(d$fx_any), NA_integer_, ifelse(d$fx_any == 0, 0L, ifelse(n_low > 0, 1L, ifelse(n_unk > 0, NA_integer_, 0L))))
cat(sprintf("  Low-energy fracture: positive %d; set to missing as mechanism undeterminable %d (explicit 8/9 %d; fracture at 50+ with missing mechanism %d)\n",
            sum(d$fx_frag == 1, na.rm = TRUE), sum(is.na(d$fx_frag)) - sum(is.na(d$fx_any)),
            sum(d$fx_any == 1 & rowSums(MEC %in% c(8, 9) & !is.na(MEC)) > 0, na.rm = TRUE),
            sum(d$fx_any == 1 & rowSums(is.na(MEC) & !is.na(AGE) & AGE >= 50) > 0, na.rm = TRUE)))
for (v in c("HAH1", "HAH2", "HAH8", "HAH9")) { x <- g(PF, v); x[x %in% c(8, 9)] <- NA; d[[paste0("dif_", v)]] <- ifelse(x >= 2, 1L, ifelse(x == 1, 0L, NA)) }
d$tand_fail <- as.integer(d$tandem < 10)
cnt <- function(x) c(sum(x == 1, na.rm = TRUE), sum(!is.na(x)))
EXP <- list(fx_hip = c(143, 4875), fx_wri = c(422, 4875), fx_spi = c(103, 4874), fx_frag = c(245, 4838), fx_any = c(614, 4875), fall = c(1113, 4852),   # fx_frag denominator 4838: undeterminable mechanism set to missing
            dif_HAH1 = c(1401, 4720), dif_HAH2 = c(1341, 4716), dif_HAH8 = c(344, 4855), dif_HAH9 = c(1026, 4853))
# Reference counts for the four activity-difficulty items are among those with non-missing BMI (16 fewer), so compare in that subset
for (v in names(EXP)) { w <- if (grepl("^dif_", v)) !is.na(d$bmi) else rep(TRUE, nrow(d)); k <- cnt(d[[v]][w]); kr <- cnt(d[[v]])
  cat(sprintf("  %-9s raw %5d / %5d; comparison basis %5d / %5d (old Table 4: %d / %d)\n", v, kr[1], kr[2], k[1], k[2], EXP[[v]][1], EXP[[v]][2])); stopifnot(all(k == EXP[[v]])) }
cat(sprintf("  tand_fail %5d / %5d; gait n %d\n", cnt(d$tand_fail)[1], cnt(d$tand_fail)[2], sum(!is.na(d$gait))))
write.csv(d[, c("SEQN", "fx_wri", "fx_spi", "fx_frag", "tand_fail", "dif_HAH1", "dif_HAH2", "dif_HAH8", "dif_HAH9")], file.path(out, "fxf_vars.csv"), row.names = FALSE)
