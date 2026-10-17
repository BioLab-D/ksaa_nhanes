res <- here::here("results"); ana <- file.path(res, "analysis"); td <- file.path(ana, "tables")
dir.create(td, showWarnings = FALSE)
mi <- function(x) gsub("-", "−", x, fixed = TRUE)
# conventional rounding (as in fig_common.R)
f <- function(x, d = 2) ifelse(is.na(x), "", mi(sub("^-(0\\.0+)$", "\\1", formatC(round(x + sign(x) * 1e-9, d), format = "f", digits = d))))
cm <- function(x) formatC(round(x), format = "d", big.mark = ",")
th <- function(h) paste0("<thead><tr>", paste0(sprintf("<th%s>%s</th>", ifelse(seq_along(h) == 1, "", ' class="num"'), h), collapse = ""), "</tr></thead>")
row <- function(cells, ind = FALSE, cls = "") {
  paste0(sprintf('<tr class="%s">', cls),
         paste0(sprintf('<td class="%s">%s</td>', c(if (ind) "ind" else "", rep("num", length(cells) - 1)), cells), collapse = ""), "</tr>")
}
hdr <- function(txt, k) sprintf('<tr class="hdr"><td colspan="%d">%s</td></tr>', k, txt)
put <- function(name, head, body, notes, thead = NULL) {   # thead: optional custom (two-level) header, default th(head)
  html <- paste0('<div class="twrap"><table>', if (is.null(thead)) th(head) else thead, "<tbody>", paste0(body, collapse = ""), "</tbody></table></div>",
                 '<div class="tnote">', paste0("<p>", notes, "</p>", collapse = ""), "</div>")
  html <- gsub("<(?=[0-9=\\s])", "&lt;", html, perl = TRUE)
  html <- gsub("[Kk]nee age gap", "KSAA", html)   # label the exposure KSAA in the HTML only; the line above escapes "<" outside tags (BMI <25, etc.)
  writeLines(html, file.path(td, paste0(name, ".html")), useBytes = TRUE); cat("->", name, "\n")
}
d <- readRDS(file.path(ana, "analysis_data.rds")); qc <- quantile(d$gap, 1:3 / 4)
PAW <- quantile(d$pa_accel_raw, c(.01, .99), na.rm = TRUE)   # winsorization limits of Phenotypic Age acceleration

# ---------- prediction performance of the backbones and the ensemble ----------
A5 <- c("ResNet-18" = "resnet18", "ResNet-50" = "resnet50", "EfficientNet-B0" = "efficientnetb0", "ConvNeXt-Tiny" = "convnexttiny",
        "DenseNet-121" = "densenet121")
met <- function(p, a) {
  e <- p - a
  c(MAE = mean(abs(e)), RMSE = sqrt(mean(e^2)), CS5 = 100 * mean(abs(e) <= 5), r = cor(p, a), R2 = 1 - sum(e^2) / sum((a - mean(a))^2))
}
fm <- function(p) {   # mean +/- SD over the five held-out folds
  m <- sapply(0:4, function(k) met(p[d$fold == k], d$age[d$fold == k]))
  c(sprintf("%s ± %s", f(mean(m["MAE", ]), 2), f(sd(m["MAE", ]), 2)), sprintf("%s ± %s", f(mean(m["RMSE", ]), 2), f(sd(m["RMSE", ]), 2)),
    sprintf("%s ± %s", f(mean(m["CS5", ]), 1), f(sd(m["CS5", ]), 1)), sprintf("%s ± %s", f(mean(m["r", ]), 3), f(sd(m["r", ]), 3)))
}
body <- c()
for (nm in names(A5)) {
  x <- read.csv(file.path(res, "predictions", sprintf("cv_%s_oof.csv", A5[[nm]])))
  p <- x$pred_age[match(d$SEQN, x$SEQN)]; stopifnot(!anyNA(p))
  body <- c(body, row(c(nm, fm(p))))
}
body <- c(body, row(c("Weighted ensemble<sup>a</sup>", fm(d$pred_raw)), cls = "bold"))
# Metrics use the uncorrected prediction, because the bias correction itself uses chronological age
me <- sapply(0:4, function(k) met(d$pred_raw[d$fold == k], d$age[d$fold == k])["MAE"])
base <- rep(NA_real_, nrow(d)); for (k in 0:4) base[d$fold == k] <- median(d$age[d$fold != k])
bm <- sapply(0:4, function(k) { i <- d$fold == k; mean(abs(base[i] - d$age[i])) })
put("tableS1", c("Model", "MAE, years", "RMSE, years", "CS (α = 5 years), %", "Pearson r"), body, c(
  sprintf(paste("Mean ± SD across the five held-out folds of participant-level cross-fitting (n = %s). All metrics are for the uncorrected",
                "prediction, because the age-bias correction uses chronological age (Methods). The ensemble MAE ranged from %s to %s years across folds;",
                "assigning every participant the median age of the training folds gave an MAE of %s years."),
          cm(nrow(d)), f(min(me), 2), f(max(me), 2), f(mean(bm), 2)),
  "<sup>a</sup>Backbone predictions weighted by the inverse internal-validation MAE within the training folds; the primary prediction."))

# ---------- baseline characteristics by KSAA quartile ----------
T1 <- read.csv(file.path(ana, "quartile_characteristics.csv"), stringsAsFactors = FALSE); T1[is.na(T1)] <- ""
nq <- T1[T1$lab == "Participants, n", ]
# drop section headers and the deaths block here; gq() below still reads the death rows from the CSV
drop_rows <- c("All-cause deaths", "Cardiovascular deaths", "Cancer deaths")
T1 <- T1[T1$type != "hdr" & !(T1$lab %in% drop_rows), ]
lab3 <- function(x) {   # harmonize units and row labels (quartile_characteristics.csv uses different wording)
  x <- sub(", y$", ", years", x); x <- sub("Poverty-income", "Poverty–income", x, fixed = TRUE)
  sub("hip, wrist, or spine fracture", "hip, wrist or vertebral fracture", x, fixed = TRUE)
}
body <- c()
for (i in seq_len(nrow(T1))) {
  r <- T1[i, ]
  body <- c(body, if (r$type == "hdr") hdr(r$lab, 6) else row(c(lab3(gsub("<sup>[^<]*</sup>", "", r$lab)), r$Q1, r$Q2, r$Q3, r$Q4, r$p), r$ind == 1))
}
body <- body[!grepl("Participants, n", body, fixed = TRUE)]
put("tableS3", c("Characteristic", sprintf("Q1 (n = %s)", nq$Q1), sprintf("Q2 (n = %s)", nq$Q2), sprintf("Q3 (n = %s)", nq$Q3),
                 sprintf("Q4 (n = %s)", nq$Q4), "P"), body, c(
  sprintf(paste("Data are presented as mean (SD) or n (%%) and are unweighted. Percentages were calculated among participants with available",
                "data. Phenotypic Age acceleration was winsorized at the 1st and 99th percentiles (%s and %s years), as in the regression analyses.",
                "KSAA quartile cut points were %s, %s, and %s years. P values were obtained using survey-weighted Rao–Scott tests for categorical",
                "variables and design-based Kruskal–Wallis tests for continuous variables."),
          f(PAW[1]), f(PAW[2]), f(qc[1]), f(qc[2]), f(qc[3])),
  paste("Abbreviations: BMD, bone mineral density; eGFR, estimated glomerular filtration rate; HDL, high-density lipoprotein; KSAA, knee skeletal",
        "age acceleration; OA, osteoarthritis; Q1–Q4, first through fourth quartiles; SD, standard deviation.")))

# ---------- Table 1 (fxf_vars.csv is written by 04_fracture_function_vars.R) ----------
V <- read.csv(file.path(ana, "fxf_vars.csv")); d1 <- d
for (v in setdiff(names(V), "SEQN")) d1[[v]] <- V[[v]][match(d1$SEQN, V$SEQN)]
msd <- function(x, dig) { x <- x[!is.na(x)]; sprintf("%s (%s)", f(mean(x), dig), f(sd(x), dig)) }
pb <- function(x) { x <- x[!is.na(x)]; sprintf("%s (%s)", cm(sum(x == 1)), f(100 * mean(x == 1), 1)) }
# the overall column of quartile_characteristics.csv must agree with Table 1
T1Q <- read.csv(file.path(ana, "quartile_characteristics.csv"))
gq <- function(lab) { v <- T1Q$overall[T1Q$lab == lab]; stopifnot(length(v) == 1); v }
chk <- list(c("Age, years", msd(d1$age, 1)), c("Sex, female", pb(as.integer(d1$sex == 2))), c("Body mass index, kg/m²", msd(d1$bmi, 1)),
  c("Diabetes<sup>b</sup>", pb(d1$dm)), c("Total cholesterol, mmol/L", msd(d1$tc / 38.67, 2)), c("HDL cholesterol, mmol/L", msd(d1$hdl / 38.67, 2)),
  c("eGFR, mL/min/1.73 m²<sup>b</sup>", msd(d1$egfr, 1)), c("Phenotypic Age acceleration, years<sup>c</sup>", msd(d1$pa_accel, 1)),
  c("Total hip BMD, g/cm²<sup>d</sup>", msd(d1$tot_bmd, 3)), c("Radiographic knee OA<sup>e</sup>", pb(as.integer(d1$klmax >= 2))), c("Knee pain", pb(d1$pain)),
  c("Gait speed, m/s", msd(d1$gait, 2)), c("≥1 fall in the past 12 months", pb(d1$fall)), c("Knee age gap, mean (SD), y", msd(d1$gap, 2)),
  c("Cardiovascular deaths", pb(d1$e_cvd)))
for (z in chk) if (gq(z[1]) != z[2]) stop("Table 1 and the Table S3 overall column disagree: ", z[1], " ", gq(z[1]), " vs ", z[2])
# Table 1: missing data in a footnote; education in years (cleaned HFA8R from phewas_clean.rds); race as non-Hispanic White
CLN <- readRDS(file.path(ana, "phewas", "phewas_clean.rds")); stopifnot(identical(CLN$SEQN, d1$SEQN))
d1$edu_y <- CLN$x[["HFA8R"]]
NMISS <- list()   # rows with missing data, listed in a footnote
mrow <- function(lab, x, dig, ind = FALSE) {
  if (sum(!is.na(x)) < nrow(d1)) NMISS[[lab]] <<- sum(!is.na(x))
  row(c(lab, msd(x, dig)), ind)
}
brow <- function(lab, x, ind = FALSE) {
  if (sum(!is.na(x)) < nrow(d1)) NMISS[[lab]] <<- sum(!is.na(x))
  row(c(lab, pb(x)), ind)
}
body <- c(hdr("Demographics and socioeconomic status", 2),
  mrow("Age, years", d1$age, 1), brow("Female sex", as.integer(d1$sex == 2)),
  brow("Non-Hispanic White", as.integer(d1$race == 1)),
  mrow("Education, years", d1$edu_y, 1),
  brow("Family income below the poverty line", d1$poor),
  hdr("Lifestyle and cardiometabolic factors (Model 2 and 3 covariates)", 2),
  brow("Former or current smoking", as.integer(d1$smk3 %in% c("former", "current"))),
  brow("Current alcohol use", as.integer(d1$alc3 == "current")),
  brow("Physically inactive in leisure time", as.integer(d1$pa3 == "inactive")),
  mrow("Body mass index, kg/m²", d1$bmi, 1), brow("Hypertension", d1$htn), brow("Diabetes", d1$dm),
  mrow("Total cholesterol, mmol/L", d1$tc / 38.67, 2), mrow("HDL cholesterol, mmol/L", d1$hdl / 38.67, 2),   # SI units (mg/dL / 38.67)
  mrow("eGFR, mL/min/1.73 m²", d1$egfr, 1),
  brow("History of cardiovascular disease", d1$cvd_hist), brow("History of cancer", d1$cancer_hist),
  hdr("Knee skeletal age and comparison measures", 2),
  row(c("KSAA, years, range (quartile cut points)", sprintf("%s to %s (%s, %s, %s)", f(min(d1$gap), 1), f(max(d1$gap), 1),
        f(quantile(d1$gap, .25), 2), f(quantile(d1$gap, .5), 2), f(quantile(d1$gap, .75), 2)))),
  mrow("Phenotypic Age acceleration, years", d1$pa_accel_raw, 1),   # raw values for description (analyses use values winsorized at 1/99)
  brow("Radiographic knee OA, K–L grade ≥2 (phase 2)", as.integer(d1$klmax >= 2)), brow("Knee pain", d1$pain),
  mrow("Total hip BMD, g/cm²", d1$tot_bmd, 3),
  hdr("Fracture, falls and physical function", 2),
  brow("Prior hip fracture (self-reported)", d1$fx_hip), brow("Any prior hip, wrist or vertebral fracture (self-reported)", d1$fx_any),
  brow("Fall in the past 12 months", d1$fall), mrow("Gait speed, m/s", d1$gait, 2),
  hdr("Mortality through 31 December 2019", 2),
  row(c("Follow-up time, years, median (IQR)",   # observed follow-up (PERMTH_EXM/12)
        sprintf("%s (%s–%s)", f(median(d1$t / 12), 1), f(quantile(d1$t / 12, .25), 1), f(quantile(d1$t / 12, .75), 1)))),
  brow("All-cause deaths", d1$e), brow("Cardiovascular deaths", d1$e_cvd, TRUE), brow("Cancer deaths", d1$e_cancer, TRUE))
miss_note <- sprintf("Numbers with data where fewer than %s: %s.", cm(nrow(d1)),
                     paste(sprintf("%s, %s", sub(",.*$", "", sub(" \\(.*$", "", names(NMISS))), cm(unlist(NMISS))), collapse = "; "))
put("table1", c("Characteristic", sprintf("Overall (n = %s)", cm(nrow(d1)))), body, c(
  paste("Data are presented as mean (SD), n (%), or median (IQR). All descriptive statistics are unweighted. Percentages were calculated",
        "among participants with available data. Observed follow-up was measured from the examination to death or censoring; the median",
        "potential follow-up (reverse Kaplan–Meier method) was 28.1 years. Radiographic knee OA was defined as a K–L grade ≥2 in the more",
        "severely affected knee; K–L readings were available only in NHANES III phase 2."),
  paste("Abbreviations: BMD, bone mineral density; eGFR, estimated glomerular filtration rate; HDL, high-density lipoprotein; IQR, interquartile",
        "range; K–L, Kellgren–Lawrence; OA, osteoarthritis; SD, standard deviation.")))
