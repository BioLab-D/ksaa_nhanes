ew <- here::here("results", "analysis", "phewas")
A <- read.csv(file.path(ew, "phewas_audit_counts.csv"), stringsAsFactors = FALSE)
cand <- read.csv(file.path(ew, "phewas_candidates.csv"), stringsAsFactors = FALSE)
A$in_pool <- A$variable %in% cand$variable[cand$analysed == "TRUE" | cand$analysed == TRUE]
lab <- tolower(A$label)
sec <- A$section
sub <- A$subsection

# ---------- 1. Health class: by questionnaire section, then by label ----------
OBJ_SEC <- c("HEMATOLOGY", "GENERAL BIOCHEMISTRY TESTS", "BIOCHEMISTRY PROFILE", "DIABETES TESTING PROFILE", "URINE TESTS", "BODY MEASUREMENTS",
             "BONE DENSITOMETRY", "SPIROMETRY", "GALLBLADDER ULTRASONOGRAPHY", "FUNDUS PHOTOGRAPHY", "DENTAL EXAMINATION", "PHYSICIAN'S EXAMINATION",
             "CENTRAL NERVOUS SYSTEM FUNCTION EVALUATION", "ALLERGY SKIN TEST", "ANTIBODY TESTS")
cls <- rep("unclassified", nrow(A))
cls[sec %in% OBJ_SEC] <- "objective"
cls[sec == "ANTIBODY TESTS"] <- "disease"
cls[sec == "PHYSICAL FUNCTION EVALUATION"] <- "function"
cls[grepl("DIAGNOSTIC INTERVIEW SCHEDULE", sec)] <- "disease"
cls[grepl("SECTION H\\. PHYSICAL FU", sub)] <- "function"
cls[grepl("SECTION [CDEJK]\\. (SELECTED CO|DIABETES|HIGH BLOOD|GALLBLADDER|KIDNEY)", sub)] <- "disease"
cls[grepl("SECTION G\\. MUSCULOSKEL|SECTION L\\. RESPIRATORY|SECTION F\\. CARDIOVASCU|SECTION P\\. VISION", sub)] <- "symptom"
cls[grepl("SECTION [NM]\\. (FOOD|DIET)|SECTION R\\. TOBACCO|SECTION T\\. EXERCISE|MEC YOUTH", sub)] <- "exposure"
cls[grepl("SECTION E\\. HOUSING", sub)] <- "exposure"                                     # housing environment -> extended module
cls[grepl("SECTION S\\. OCCUPATION|SECTION F\\. FAMILY CHA|SECTION B\\. HEALTH INS|GENERAL INFORMATION|SECTION Y|SECTION A\\.", sub) |
      grepl("HOUSEHOLD FAMILY", sec)] <- "covariate"
cls[grepl("REPRODUCTI", sub) | grepl("REPRODUCTIVE HEALTH", sub)] <- "reproductive"
cls[grepl("SECTION Q\\. DENTAL CARE|SECTION B\\. (HEALTH CARE|HEALTH SERV|GENERAL)", sub)] <- "healthcare"
# cognitive tests (orientation, serial subtraction, story recall) = objective performance
cls[grepl("COGNITIVE FUNCTION", sub) | grepl("^HAP18|^HAA[456]", A$variable)] <- "objective"
# home-visit blood pressure/pulse readings
cls[grepl("SECTION Z\\. BLOOD PRESSURE", sub)] <- "objective"
# social support: not SES; extended module
cls[grepl("SECTION V\\. SOCIAL SUPPORT", sub)] <- "psychosocial"
cls[grepl("PHLEBOTOMY", sub)] <- "non_outcome"
cls[grepl("^HAY11", A$variable)] <- "objective"                                                    # interviewer-observed signs
# physician-estimated activity difficulty: physical function
cls[grepl("^PEP13E[1-5]A$", A$variable)] <- "function"
# cognitive test items (registration, serial subtraction)
cls[grepl("^HAP17|^HAP18", A$variable)] <- "objective"
# hearing: deafness, hearing aid, hearing difficulty
cls[grepl("^HAP1[0-3]$", A$variable)] <- "symptom"
# label-based overrides within questionnaires
q <- !(sec %in% OBJ_SEC) & !(cls %in% c("objective", "psychosocial", "non_outcome")) & !grepl("SECTION E\\. HOUSING", sub)
# behaviour/diet keywords first; later disease/function/symptom rules may override (e.g. HAH6, difficulty preparing meals, is a function item)
cls[q & grepl("smok|cigar|tobacco|snuff|chew|alcohol|drink|beer|wine|liquor|vitamin|supplement|calcium tablet|exercise|jog|swim|garden|bicycl|dance|calisthen|lift weights|physical activ|active compared|walk .*mile|milk|\\bdiet|\\beat\\b|\\bfood|\\bmeal|coffee|\\btea\\b", lab)] <- "exposure"
cls[q & grepl("doctor (ever )?told|told by a doctor|ever been told|ever told you|diagnos|ever had|have you had|ever have", lab)] <- "disease"
cls[q & grepl("difficult|limited in|unable to|need help|use device|cane|walker|wheelchair|change work|housework, health", lab)] <- "function"
# "\bache" has a word boundary so "peaches" does not match
cls[q & grepl("\\bpain\\b|\\bache|stiff|swell|numbness|\\bnumb\\b|tingl|cough|phlegm|wheez|short of breath|breathless|dizz|limp\\b", lab)] <- "symptom"
cls[grepl("health in general|health status|excellent,", lab)] <- "self_rated"
cls[q & grepl("nicotine gum", lab)] <- "exposure"
cls[grepl("cotinine|\\blead\\b|cadmium|mercury|carotene|lycopene|lutein|cryptoxanthin|vitamin [acde]\\b|vitamin b12|folate|retinyl|selenium", lab)] <- "exposure_biomarker"
cls[q & grepl("income|educat|school|insur|medicare|medicaid|\\bjob\\b|occupation|employ|marital|born|language|household|family size|military|lived|address", lab)] <- "covariate"
cls[grepl("cuff|number of (bp|blood pressure) readings|position|examiner|proxy|who answered|record type|time of (exam|day)|exam status|reason not|arm circumference for", lab)] <- "non_outcome"
# physician prognosis: not a current-state phenotype
cls[grepl("prognosis|prog over 1 yr", lab)] <- "review_prognosis"
cls[grepl("consider self over|like to weigh", lab)] <- "review_self_perception"
cls[A$variable %in% c("PEP13A", "HAB1")] <- "self_rated"
cls[A$variable == "PHPINSU"] <- "disease"
cls[grepl("^HAE8", A$variable)] <- "exposure"
# HAE advice items: lifestyle advice -> exposure; medication (HAE4A/5A antihypertensives, HAE8D/9D lipid-lowering) -> disease
cls[A$variable %in% c("HAE4B", "HAE4C", "HAE4D", "HAE5B", "HAE5C", "HAE9B", "HAE9C")] <- "exposure"
cls[A$variable %in% c("HAE4A", "HAE5A", "HAE8D", "HAE9D")] <- "disease"
# home age/number of rooms/years at residence: not exposures
cls[A$variable %in% c("HFE1", "HFE2", "HFE3")] <- "covariate"
# "wearing street clothes" observation: unclear meaning
cls[A$variable == "HAY11I"] <- "review_observation"
A$health_class <- cls

# ---------- 2. Duplicates / compositional overlap ----------
fam <- gsub("\\s*(trial [ab]|1st|2nd|3rd|first|second|third|reading|replicate|average|mean of)\\s*", " ", lab)
fam <- gsub("[[:punct:]]+", " ", fam)
fam <- gsub("\\s+", " ", trimws(fam))   # strip punctuation before grouping ("walk, trial A" and "walk--trial B" are one family)
A$dup_family <- ifelse(ave(seq_along(fam), fam, FUN = length) > 1 & grepl("trial [ab]|1st|2nd|3rd|reading|replicate", lab), fam, "")
A$threshold_version <- grepl(">=|≥|\\bcut\\b|categor|group\\b|abnormal|\\bhigh\\b|\\blow\\b", lab) & A$type == "categorical" & sec %in% OBJ_SEC
# no readings of the same radiographs in the pool (K-L etc. are not in it); KSAA is built only from knee films and chronological age
A$overlap_ksaa <- "no"
# ---------- 3. Eligible population / phase / weight ----------
AGE_RX <- "\\((age|ages?|yrs?)?\\s*[0-9]+\\s*(\\+|and over|-\\s*[0-9]+|to [0-9]+|years and over|yrs and over| years)?[^)]*\\)|\\(all ages\\)"
m <- regexpr(AGE_RX, A$label, ignore.case = TRUE)
A$target_pop <- ""
A$target_pop[m > 0] <- regmatches(A$label, m)
A$phase <- ifelse(A$n_valid_phase1 < 50 & A$n_valid_phase2 >= 50, "phase 2 only",
    ifelse(A$n_valid_phase2 < 50 & A$n_valid_phase1 >= 50, "phase 1 only", "both"))
A$weight_basis <- ifelse(sec == "DIABETES TESTING PROFILE", "fasting subsample weight (WTPFSD6)",
    "WTPFEX6 (examined cohort), PSU SDPPSU6, strata SDPSTRA6")
# ---------- 4. Data-quality thresholds (not based on P values) ----------
size_ok <- with(A, ifelse(type == "continuous", n_valid >= 500 & n_distinct >= 5,
                   ifelse(n_levels == 2, n_valid >= 500 & min_level_n >= 100, n_valid >= 500)))
A$sparse_level <- A$type == "categorical" & A$n_levels > 2 & A$min_level_n < 30
A$validity_ok <- !is.na(A$validity_applicable) & A$validity_applicable >= 0.70
size_ok[is.na(size_ok)] <- FALSE
A$size_ok <- size_ok
# ---------- 5. Inclusion status ----------
core <- A$health_class %in% c("objective", "disease", "symptom", "function", "self_rated", "healthcare")
# reproductive history -> extended module (women only); physician 1-year prognosis, self-rated weight, "street clothes" observation -> excluded
A$health_class[A$health_class == "reproductive"] <- "reproductive_exposure"
A$target_pop[A$health_class == "reproductive_exposure"] <- "women only"
A$health_class[A$health_class %in% c("review_prognosis", "review_self_perception", "review_observation")] <- "non_outcome"
# repeated measures: single home BP readings excluded (means HAZMNK1R/K5R used); story-recall segments and address
#   orientation items are summed into scores in 10_phewas.R
# single readings and separate home/MEC means: use HAZMNK1R/K5R (combined home + MEC means)
A$health_class[grepl("^HAZA8[ABC]K[15]|^HAZA11[AB]R|^PEPMNK[15]R$|^PEP6[GHI][13]$", A$variable)] <- "non_outcome"
A$summarise_to <- ifelse(grepl("1st story recall", lab), "MAPB_story_recall_immediate_total",
                  ifelse(grepl("2nd story recall", lab), "MAPD_story_recall_delayed_total",
                  ifelse(grepl("^HAA6", A$variable), "HAA_orientation_score", "")))
# other summaries built in 10_phewas.R: gait speed (faster trial), mean steps and pulse, number correct in serial subtraction
#   and three-word registration/recall
SUMM <- c(PFPTWLKA = "PFP_gait_speed", PFPTWLKB = "PFP_gait_speed",
          PFPNSTPA = "PFP_steps_8ft_mean", PFPNSTPB = "PFP_steps_8ft_mean",
          HAZA5R = "pulse_mean", PEP6DR = "pulse_mean",
          HAP18A = "HAP18_serial_subtraction", HAP18B = "HAP18_serial_subtraction", HAP18C = "HAP18_serial_subtraction",
          HAP18D = "HAP18_serial_subtraction", HAP18E = "HAP18_serial_subtraction",
          HAP17A1 = "HAP17_word_registration", HAP17A2 = "HAP17_word_registration", HAP17A3 = "HAP17_word_registration",
          HAP19A = "HAP19_word_recall", HAP19B = "HAP19_word_recall", HAP19C = "HAP19_word_recall")
A$summarise_to[A$variable %in% names(SUMM)] <- SUMM[A$variable[A$variable %in% names(SUMM)]]
# in-place binary recodes matching the Table S10 definitions (each variable once; applied in 10_phewas.R)
S8V <- c(setNames(rep("binary: any difficulty vs none", 17), c(paste0("HAH", 1:12), paste0("PEP13E", 1:5, "A"))),
         HAG14 = "binary: at least one fall", PFPTTAND = "binary: unable to hold 10 s")
A$s8_version <- ifelse(A$variable %in% names(S8V), S8V[A$variable], "")
A$status <- ifelse(A$health_class %in% c("exposure", "exposure_biomarker", "psychosocial", "reproductive_exposure"), "exposure_module",
            ifelse(A$health_class == "covariate", "covariate",
            ifelse(A$health_class == "non_outcome", "excluded",
            ifelse(grepl("^review", A$health_class) | A$health_class %in% c("reproductive", "unclassified"), "review",
            ifelse(core & size_ok & A$validity_ok, "main", "supplementary")))))
A$reason <- ifelse(A$status == "supplementary",
  trimws(paste(ifelse(!size_ok, "sample/case size below threshold", ""),
               ifelse(!A$validity_ok, sprintf("validity %.0f%% < 70%% of applicable", 100 * A$validity_applicable), ""))),
  ifelse(A$status == "review", A$health_class,
  ifelse(A$status == "excluded",
    ifelse(grepl("^HAZA|^PEP6|^PEPMNK", A$variable), "single BP reading / partial mean: combined mean HAZMNK1R/K5R used instead",
    ifelse(grepl("^PEP13E", A$variable), "physician prognosis: judgement, not a state",
    ifelse(grepl("^HAM1[12]|HAY11I", A$variable), "self-perception / unclear observation item",
           "not an outcome (QC/process field)"))),
    "")))
A$reason[A$status == "main" & A$sparse_level] <- "sparse category (<30): consider merging if meaningful"
A$reason[A$status == "main" & A$dup_family != ""] <- paste0("repeated measure family: summarise per protocol [",
    A$dup_family[A$status == "main" & A$dup_family != ""], "]")
cg <- grepl("can't grade|could not be assessed", A$levels, ignore.case = TRUE)
ud <- grepl("unable to do|can't do", A$levels, ignore.case = TRUE)
A$note_cant_do <- ifelse(cg, sprintf("RECODE: 'can't grade' kept as a level in cleaning (n=%d) -> must be missing", A$n_cant_do),
                  ifelse(ud, sprintf("'unable to do' retained as top level (n=%d): real outcome, keep", A$n_cant_do),
                  ifelse(A$n_cant_do > 0,
                         sprintf("continuous: 'cannot do' set to missing (n=%d); consider unable-as-outcome sensitivity", A$n_cant_do),
                         "")))

OUT <- A[, c("variable", "label", "health_class", "type", "n_levels", "levels", "target_pop", "n_applicable", "phase",
             "overlap_ksaa", "validity_applicable", "validity_cohort", "n_valid", "min_level_n", "n_cant_do", "size_ok",
             "validity_ok", "weight_basis", "in_pool", "status", "reason", "sparse_level", "dup_family", "threshold_version",
             "summarise_to", "s8_version", "note_cant_do", "section", "subsection")]
OUT$audit_version <- "v1.2"
write.csv(OUT, file.path(ew, "phewas_audit_table.csv"), row.names = FALSE)

P <- OUT[OUT$in_pool, ]
cat(sprintf("=== Candidate pool (n = %d) by status x health class ===\n", nrow(P)))
print(table(P$health_class, P$status))
cat("\n=== Status totals in the pool ===\n")
print(table(P$status))
cat("\n=== Outside the pool (500 <= n_valid < 1000, within rules) by status ===\n")
O <- OUT[!OUT$in_pool & OUT$n_valid >= 500, ]
print(table(O$status))
cat("Of which would enter the main analysis:", sum(O$status == "main"), "\n")
cat("\n=== Reasons for supplementary in the pool ===\n")
print(sort(table(P$reason[P$status == "supplementary"]), decreasing = TRUE))
cat("\n=== Pending review in the pool ===\n")
print(table(P$health_class[P$status == "review"]))
cat("\n=== Top 15 main-class variables in the pool with valid rate < 70% (within the eligible population) ===\n")
x <- P[P$health_class %in% c("objective","disease","symptom","function","self_rated") & !P$validity_ok, ]
x <- x[order(x$validity_applicable), ]
print(head(x[, c("variable", "label", "n_applicable", "n_valid", "validity_applicable")], 15), row.names = FALSE)
cat("\n=== Continuous measures with 'cannot do' cleaned to missing ===\n")
print(P[P$note_cant_do != "", c("variable", "label", "note_cant_do")], row.names = FALSE)
cat("\n=== Repeated-measure families ===\n")
print(unique(P$dup_family[P$dup_family != ""]))
cat("\n=== Phase-specific ===\n")
print(table(P$phase))
cat("\n=== Unclassified in the pool ===\n")
print(P[P$health_class == "unclassified", c("variable","label","subsection")], row.names = FALSE)
cat("\n=== Can't-grade variables to set to missing ===\n")
print(P[grepl("^RECODE", P$note_cant_do), c("variable","label","n_cant_do")], row.names = FALSE)
cat("\n=== Main analysis: class x type ===\n")
print(table(P$health_class[P$status=="main"], P$type[P$status=="main"]))
cat("\n=== Review items in the pool ===\n")
print(P[P$status == "review" & P$health_class != "unclassified", c("variable","label","health_class")], row.names = FALSE)
cat("\n=== Covariate items in the pool ===\n")
print(P[P$status == "covariate", c("variable","label")], row.names = FALSE)
cat("\n=== exposure_biomarker in the pool ===\n")
print(P[P$health_class == "exposure_biomarker", c("variable","label")], row.names = FALSE)
