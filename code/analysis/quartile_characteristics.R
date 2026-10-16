# Baseline characteristics by knee age gap quartile (used by Table S3 and Table 1): unweighted means / percentages; design-based P values: Rao–Scott F (svychisq) for categorical,
#   design-based Kruskal–Wallis (svyranktest) for continuous variables; design on all adult examinees with the cohort as a domain (weighted_common.R), phase 2 design for the K–L reading.
# Input results/analysis/analysis_data.rds, fxf_vars.csv (read by weighted_common.R). Run by run_all.sh; output quartile_characteristics.csv
res <- here::here("results"); out <- file.path(res, "analysis")
d <- readRDS(file.path(out, "analysis_data.rds"))
qc <- quantile(d$gap, 0:4 / 4)
d$q <- cut(d$gap, qc, include.lowest = TRUE, labels = paste0("Q", 1:4)); d$qi <- as.integer(d$q)
mi_ <- function(x) gsub("-", "−", x, fixed = TRUE)
f <- function(x, dg = 1) mi_(sub("^-(0\\.0+)$", "\\1", formatC(round(x + sign(x) * 1e-9, dg), format = "f", digits = dg)))
cm <- function(x) formatC(round(x), format = "d", big.mark = ",")
fp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "&lt;0.001", formatC(p, format = "f", digits = 3)))

WP <- local({   # design-based P values (by variable name); participants missing a variable are excluded from its domain
  source(here::here("code", "analysis", "weighted_common.R"), local = TRUE)
  D6 <- update(D6, q = factor(q, QL), sexf = factor(HSSEX), racef = factor(DMARETHN)); D2 <- update(D2, q = factor(q, QL))
  stopifnot(all(F$HSSEX[F$inc] == F$sex[F$inc]), all(F$DMARETHN[F$inc] == F$race[F$inc]))   # original NHANES variables agree with cohort variables
  kw <- function(v) svyranktest(as.formula(paste(v, "~ q")), dom6(F$inc & !is.na(F[[v]])), test = "KruskalWallis")$p.value[[1]]
  rs <- function(v, ok = F[[v]], des = dom6) svychisq(as.formula(paste("~", v, "+ q")), des(F$inc & !is.na(ok)), statistic = "F")$p.value[[1]]
  c(sapply(c(age = "age", bmi = "bmi", tc = "tc", hdl = "hdl", egfr = "egfr", pa_accel = "pa_accel", fn_bmd = "fn_bmd", tot_bmd = "tot_bmd", gait = "gait"), kw),
    sex = rs("sexf", F$sex), race = rs("racef", F$race),
    sapply(c(edu3 = "edu3", poor = "poor", smk3 = "smk3", alc3 = "alc3", pa3 = "pa3", htn = "htn", dm = "dm", cvd_hist = "cvd_hist",
             cancer_hist = "cancer_hist", pain = "pain", anylimit = "anylimit", fx_hip = "fx_hip", fx_any = "fx_any", fall = "fall"), rs),
    oa = rs("oa", des = dom2))
}); cat("Design-based P values:", paste(sprintf("%s %.3g", names(WP), WP), collapse = "; "), "\n")

cat("=========== Baseline characteristics by knee age gap quartile (design-based P values) ===========\n")
QL <- levels(d$q); T <- list()
addr <- function(lab, type, ind, cells, p) T[[length(T) + 1]] <<- data.frame(lab = lab, type = type, ind = ind, overall = cells[1], Q1 = cells[2], Q2 = cells[3], Q3 = cells[4], Q4 = cells[5], p = p)
by_q <- function(fun, x) c(fun(x), sapply(QL, function(l) fun(x[d$q == l])))
fm <- function(dg) function(x) sprintf("%s (%s)", f(mean(x, na.rm = TRUE), dg), f(sd(x, na.rm = TRUE), dg))
fb <- function(x) sprintf("%s (%s)", cm(sum(x == 1, na.rm = TRUE)), f(100 * mean(x == 1, na.rm = TRUE), 1))
hdr <- function(lab) addr(lab, "hdr", 0, rep("", 5), "")
pw <- function(key) if (is.na(key)) "" else fp(WP[[key]])
rowc <- function(lab, v, dg = 1, ind = 0, key = v) addr(lab, "row", ind, by_q(fm(dg), d[[v]]), pw(key))
rowb <- function(lab, x, ind = 0, key = NA) addr(lab, "row", ind, by_q(fb, x), pw(key))
rowcat <- function(lab, v, levs) { addr(lab, "parent", 0, rep("", 5), pw(v))
  for (i in seq_along(levs)) rowb(levs[i], as.integer(d[[v]] == names(levs)[i]), ind = 1) }
addr("Participants, n", "row", 0, by_q(function(x) cm(length(x)), d$gap), "")
addr("Knee age gap, range, y", "row", 0, by_q(function(x) sprintf("%s to %s", f(min(x)), f(max(x))), d$gap), "")
rowc("Knee age gap, mean (SD), y", "gap", 2, key = NA)
hdr("Demographics")
rowc("Age, years", "age"); rowb("Sex, female", as.integer(d$sex == 2), key = "sex")
d$race_f <- factor(d$race); WP[["race_f"]] <- WP[["race"]]; rowcat("Race and ethnicity", "race_f", c("1" = "Non-Hispanic White", "2" = "Non-Hispanic Black", "3" = "Mexican American", "4" = "Other"))
rowcat("Education", "edu3", c("<12" = "&lt;12 years", "12" = "12 years", ">12" = "&gt;12 years"))
rowb("Poverty–income ratio &lt;1.0<sup>a</sup>", d$poor, key = "poor")
hdr("Lifestyle")
rowcat("Smoking", "smk3", c(never = "Never", former = "Former", current = "Current"))
rowcat("Alcohol use<sup>b</sup>", "alc3", c(never = "Never", former_light = "Former or light", current = "Current"))
rowcat("Leisure-time physical activity<sup>b</sup>", "pa3", c(active = "Active", insufficient = "Insufficiently active", inactive = "Inactive"))
rowc("Body mass index, kg/m²", "bmi")
hdr("Cardiometabolic factors and disease history")
# cholesterol reported in SI units (mmol/L = mg/dL ÷ 38.67); model covariates tc / hdl stay in mg/dL
d$tc_mmol <- d$tc / 38.67; d$hdl_mmol <- d$hdl / 38.67
rowb("Hypertension<sup>b</sup>", d$htn, key = "htn"); rowb("Diabetes<sup>b</sup>", d$dm, key = "dm"); rowc("Total cholesterol, mmol/L", "tc_mmol", 2, key = "tc"); rowc("HDL cholesterol, mmol/L", "hdl_mmol", 2, key = "hdl")
rowc("eGFR, mL/min/1.73 m²<sup>b</sup>", "egfr", 1); rowb("Self-reported history of cardiovascular disease<sup>b</sup>", d$cvd_hist, key = "cvd_hist"); rowb("Self-reported history of cancer<sup>b</sup>", d$cancer_hist, key = "cancer_hist")
rowc("Phenotypic Age acceleration, years<sup>c</sup>", "pa_accel", key = "pa_accel")   # winsorized at the 1st / 99th percentiles (as in the regressions)
hdr("Musculoskeletal status")
rowc("Femoral neck BMD, g/cm²<sup>d</sup>", "fn_bmd", 3); rowc("Total hip BMD, g/cm²<sup>d</sup>", "tot_bmd", 3)
rowb("Radiographic knee OA<sup>e</sup>", d$oa, key = "oa"); rowb("Knee pain", d$pain, key = "pain"); rowc("Gait speed, m/s", "gait", 2)
rowb("Any functional limitation", d$anylimit, key = "anylimit"); rowb("Self-reported history of hip fracture", d$fx_hip, key = "fx_hip"); rowb("Self-reported history of hip, wrist or vertebral fracture", d$fx_any, key = "fx_any"); rowb("≥1 fall in the past 12 months", d$fall, key = "fall")
hdr("Mortality through 31 December 2019")
rowb("All-cause deaths", d$e); rowb("Cardiovascular deaths", d$e_cvd, ind = 1); rowb("Cancer deaths", d$e_cancer, ind = 1)
T1 <- do.call(rbind, T); write.csv(T1, file.path(out, "quartile_characteristics.csv"), row.names = FALSE)
print(T1[T1$type != "hdr", c("lab", "overall", "Q1", "Q4", "p")], row.names = FALSE)
cat("\n-> quartile_characteristics.csv\n")
