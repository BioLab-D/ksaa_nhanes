# PheWAS 变量审核表(表 S8):按 健康含义 → 去重与构成重叠 → 适用人群 → 数据质量 逐变量判定纳入状态;不看效应量与 P 值。
# 输入:results/analysis/phewas/ewas_audit_counts.csv(phewas_clean.R EWAS_STAGE=audit 的逐变量计数,无模型)、ewas_candidates.csv(825 池标记)
# 输出:results/analysis/phewas/ewas_audit_table.csv,每行一个候选(健康领域、类型、适用人群、相位、有效率、样本量、纳入状态、理由等)
# 纳入状态:main(主分析)/ supplementary(补充)/ exposure_module(扩展关联模块)/ covariate / excluded / review(待人工复核)
ew <- here::here("results", "analysis", "phewas")
A <- read.csv(file.path(ew, "ewas_audit_counts.csv"), stringsAsFactors = FALSE)
cand <- read.csv(file.path(ew, "ewas_candidates.csv"), stringsAsFactors = FALSE)
A$in_pool <- A$variable %in% cand$variable[cand$analysed == "TRUE" | cand$analysed == TRUE]
lab <- tolower(A$label); sec <- A$section; sub <- A$subsection

# ---------- 1. 健康含义分类(先按模块,再按标签覆盖)----------
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
cls[grepl("SECTION E\\. HOUSING", sub)] <- "exposure"                                     # 家居环境 → 扩展模块
cls[grepl("SECTION S\\. OCCUPATION|SECTION F\\. FAMILY CHA|SECTION B\\. HEALTH INS|GENERAL INFORMATION|SECTION Y|SECTION A\\.", sub) | grepl("HOUSEHOLD FAMILY", sec)] <- "covariate"
cls[grepl("REPRODUCTI", sub) | grepl("REPRODUCTIVE HEALTH", sub)] <- "reproductive"
cls[grepl("SECTION Q\\. DENTAL CARE|SECTION B\\. (HEALTH CARE|HEALTH SERV|GENERAL)", sub)] <- "healthcare"
cls[grepl("COGNITIVE FUNCTION", sub) | grepl("^HAP18|^HAA[456]", A$variable)] <- "objective"      # 认知测验(定向、连减、故事回忆)= 客观表现
cls[grepl("SECTION Z\\. BLOOD PRESSURE", sub)] <- "objective"                                     # 家庭访视血压/脉搏读数
cls[grepl("SECTION V\\. SOCIAL SUPPORT", sub)] <- "psychosocial"                                  # 社会支持:非 SES,作扩展模块
cls[grepl("PHLEBOTOMY", sub)] <- "non_outcome"
cls[grepl("^HAY11", A$variable)] <- "objective"                                                    # 访员观察体征
cls[grepl("^PEP13E[1-5]A$", A$variable)] <- "function"                                             # 医师评估的活动困难度:身体功能
cls[grepl("^HAP17|^HAP18", A$variable)] <- "objective"                                              # 认知测验项(复述、连减)
cls[grepl("^HAP1[0-3]$", A$variable)] <- "symptom"                                                 # 听力:耳聋、助听器、听力困难
# 标签覆盖(问卷内更细的含义)
q <- !(sec %in% OBJ_SEC) & !(cls %in% c("objective", "psychosocial", "non_outcome")) & !grepl("SECTION E\\. HOUSING", sub)
cls[q & grepl("smok|cigar|tobacco|snuff|chew|alcohol|drink|beer|wine|liquor|vitamin|supplement|calcium tablet|exercise|jog|swim|garden|bicycl|dance|calisthen|lift weights|physical activ|active compared|walk .*mile|milk|\\bdiet|\\beat\\b|\\bfood|\\bmeal|coffee|\\btea\\b", lab)] <- "exposure"   # 行为/饮食关键词先行,随后的疾病/功能/症状规则可覆盖(如 HAH6 备餐困难是功能项)
cls[q & grepl("doctor (ever )?told|told by a doctor|ever been told|ever told you|diagnos|ever had|have you had|ever have", lab)] <- "disease"
cls[q & grepl("difficult|limited in|unable to|need help|use device|cane|walker|wheelchair|change work|housework, health", lab)] <- "function"
cls[q & grepl("\\bpain\\b|\\bache|stiff|swell|numbness|\\bnumb\\b|tingl|cough|phlegm|wheez|short of breath|breathless|dizz|limp\\b", lab)] <- "symptom"   # "\bache" 带词界,避免匹配 "peaches"
cls[grepl("health in general|health status|excellent,", lab)] <- "self_rated"
cls[q & grepl("nicotine gum", lab)] <- "exposure"
cls[grepl("cotinine|\\blead\\b|cadmium|mercury|carotene|lycopene|lutein|cryptoxanthin|vitamin [acde]\\b|vitamin b12|folate|retinyl|selenium", lab)] <- "exposure_biomarker"
cls[q & grepl("income|educat|school|insur|medicare|medicaid|\\bjob\\b|occupation|employ|marital|born|language|household|family size|military|lived|address", lab)] <- "covariate"
cls[grepl("cuff|number of (bp|blood pressure) readings|position|examiner|proxy|who answered|record type|time of (exam|day)|exam status|reason not|arm circumference for", lab)] <- "non_outcome"
cls[grepl("prognosis|prog over 1 yr", lab)] <- "review_prognosis"                      # 医师预后判断:非当前状态表型
cls[grepl("consider self over|like to weigh", lab)] <- "review_self_perception"
cls[A$variable %in% c("PEP13A", "HAB1")] <- "self_rated"
cls[A$variable == "PHPINSU"] <- "disease"
cls[grepl("^HAE8", A$variable)] <- "exposure"
# HAE 医嘱 / 行为题:生活方式类(减重、限盐、运动等)→ 暴露;用药类 → 疾病(HAE4A/5A 降压药,HAE8D/9D 降脂药)
cls[A$variable %in% c("HAE4B", "HAE4C", "HAE4D", "HAE5B", "HAE5C", "HAE9B", "HAE9C")] <- "exposure"
cls[A$variable %in% c("HAE4A", "HAE5A", "HAE8D", "HAE9D")] <- "disease"
cls[A$variable %in% c("HFE1", "HFE2", "HFE3")] <- "covariate"                                     # 房龄/房间数/居住年限:非暴露
cls[A$variable == "HAY11I"] <- "review_observation"                                               # "穿街装"观察项:含义不明
A$health_class <- cls

# ---------- 2. 去重 / 构成重叠 ----------
fam <- gsub("\\s*(trial [ab]|1st|2nd|3rd|first|second|third|reading|replicate|average|mean of)\\s*", " ", lab); fam <- gsub("[[:punct:]]+", " ", fam); fam <- gsub("\\s+", " ", trimws(fam))   # 先去标点再分组("walk, trial A" 与 "walk--trial B" 同族)
A$dup_family <- ifelse(ave(seq_along(fam), fam, FUN = length) > 1 & grepl("trial [ab]|1st|2nd|3rd|reading|replicate", lab), fam, "")
A$threshold_version <- grepl(">=|≥|\\bcut\\b|categor|group\\b|abnormal|\\bhigh\\b|\\blow\\b", lab) & A$type == "categorical" & sec %in% OBJ_SEC
A$overlap_ksaa <- "no"   # 池中无同片影像判读(K–L 等不在池内);KSAA 仅由膝片与实足年龄构成
# ---------- 3. 适用人群 / 相位 / 权重 ----------
age_q <- regmatches(A$label, regexpr("\\((age|ages?|yrs?)?\\s*[0-9]+\\s*(\\+|and over|-\\s*[0-9]+|to [0-9]+|years and over|yrs and over| years)?[^)]*\\)|\\(all ages\\)", A$label, ignore.case = TRUE))
A$target_pop <- ""; A$target_pop[regexpr("\\((age|ages?|yrs?)?\\s*[0-9]+\\s*(\\+|and over|-\\s*[0-9]+|to [0-9]+|years and over|yrs and over| years)?[^)]*\\)|\\(all ages\\)", A$label, ignore.case = TRUE) > 0] <- age_q
A$phase <- ifelse(A$n_valid_phase1 < 50 & A$n_valid_phase2 >= 50, "phase 2 only", ifelse(A$n_valid_phase2 < 50 & A$n_valid_phase1 >= 50, "phase 1 only", "both"))
A$weight_basis <- ifelse(sec == "DIABETES TESTING PROFILE", "fasting subsample weight (WTPFSD6; verify)", "WTPFEX6 (examined cohort), PSU SDPPSU6, strata SDPSTRA6")
# ---------- 4. 数据质量门槛(不依据 P 值)----------
size_ok <- with(A, ifelse(type == "continuous", n_valid >= 500 & n_distinct >= 5, ifelse(n_levels == 2, n_valid >= 500 & min_level_n >= 100, n_valid >= 500)))
A$sparse_level <- A$type == "categorical" & A$n_levels > 2 & A$min_level_n < 30
A$validity_ok <- !is.na(A$validity_applicable) & A$validity_applicable >= 0.70
size_ok[is.na(size_ok)] <- FALSE; A$size_ok <- size_ok
# ---------- 5. 纳入状态 ----------
core <- A$health_class %in% c("objective", "disease", "symptom", "function", "self_rated", "healthcare")
# 生殖史 → 扩展模块(仅女性);医师一年预后、自评体重、"穿街装"观察 → 排除
A$health_class[A$health_class == "reproductive"] <- "reproductive_exposure"
A$target_pop[A$health_class == "reproductive_exposure"] <- "women only"
A$health_class[A$health_class %in% c("review_prognosis", "review_self_perception", "review_observation")] <- "non_outcome"
# 重复测量族:家庭血压单次读数排除(用合并均值 HAZMNK1R/K5R);故事回忆 6 段与地址定向 4 题汇总为得分(phewas.R 构造)
A$health_class[grepl("^HAZA8[ABC]K[15]|^HAZA11[AB]R|^PEPMNK[15]R$|^PEP6[GHI][13]$", A$variable)] <- "non_outcome"   # 单次读数与家访/MEC 各自均值:统一用 HAZMNK1R/K5R(家访+MEC 合并均值)
A$summarise_to <- ifelse(grepl("1st story recall", lab), "MAPB_story_recall_immediate_total", ifelse(grepl("2nd story recall", lab), "MAPD_story_recall_delayed_total", ifelse(grepl("^HAA6", A$variable), "HAA_orientation_score", "")))
# 其余重复测量与多题测验的汇总(构造见 phewas.R):步行 → 步速(较快一次);步数、脉搏 → 均值;连减、三词复述 / 回忆 → 答对题数
SUMM <- c(PFPTWLKA = "PFP_gait_speed", PFPTWLKB = "PFP_gait_speed", PFPNSTPA = "PFP_steps_8ft_mean", PFPNSTPB = "PFP_steps_8ft_mean",
          HAZA5R = "pulse_mean", PEP6DR = "pulse_mean", HAP18A = "HAP18_serial_subtraction", HAP18B = "HAP18_serial_subtraction", HAP18C = "HAP18_serial_subtraction",
          HAP18D = "HAP18_serial_subtraction", HAP18E = "HAP18_serial_subtraction", HAP17A1 = "HAP17_word_registration", HAP17A2 = "HAP17_word_registration",
          HAP17A3 = "HAP17_word_registration", HAP19A = "HAP19_word_recall", HAP19B = "HAP19_word_recall", HAP19C = "HAP19_word_recall")
A$summarise_to[A$variable %in% names(SUMM)] <- SUMM[A$variable[A$variable %in% names(SUMM)]]
# 与表 S10 定义统一的就地二分类重编码(每个变量只出现一次;在 phewas.R 中执行)
S8V <- c(setNames(rep("binary: any difficulty vs none", 17), c(paste0("HAH", 1:12), paste0("PEP13E", 1:5, "A"))),
         HAG14 = "binary: at least one fall", PFPTTAND = "binary: unable to hold 10 s")
A$s8_version <- ifelse(A$variable %in% names(S8V), S8V[A$variable], "")
A$status <- ifelse(A$health_class %in% c("exposure", "exposure_biomarker", "psychosocial", "reproductive_exposure"), "exposure_module",
            ifelse(A$health_class == "covariate", "covariate",
            ifelse(A$health_class == "non_outcome", "excluded",
            ifelse(grepl("^review", A$health_class) | A$health_class %in% c("reproductive", "unclassified"), "review",
            ifelse(core & size_ok & A$validity_ok, "main", "supplementary")))))
A$reason <- ifelse(A$status == "supplementary", trimws(paste(ifelse(!size_ok, "sample/case size below threshold", ""), ifelse(!A$validity_ok, sprintf("validity %.0f%% < 70%% of applicable", 100 * A$validity_applicable), ""))),
            ifelse(A$status == "review", A$health_class, ifelse(A$status == "excluded", ifelse(grepl("^HAZA|^PEP6|^PEPMNK", A$variable), "single BP reading / partial mean: combined mean HAZMNK1R/K5R used instead", ifelse(grepl("^PEP13E", A$variable), "physician prognosis: judgement, not a state (user 2026-09-16)", ifelse(grepl("^HAM1[12]|HAY11I", A$variable), "self-perception / unclear observation item (user 2026-09-16)", "not an outcome (QC/process field)"))), "")))
A$reason[A$status == "main" & A$sparse_level] <- "sparse category (<30): consider merging if meaningful"
A$reason[A$status == "main" & A$dup_family != ""] <- paste0("repeated measure family: summarise per protocol [", A$dup_family[A$status == "main" & A$dup_family != ""], "]")
cg <- grepl("can't grade|could not be assessed", A$levels, ignore.case = TRUE)
ud <- grepl("unable to do|can't do", A$levels, ignore.case = TRUE)
A$note_cant_do <- ifelse(cg, sprintf("RECODE: 'can't grade' kept as a level in cleaning (n=%d) -> must be missing", A$n_cant_do),
                  ifelse(ud, sprintf("'unable to do' retained as top level (n=%d): real outcome, keep", A$n_cant_do),
                  ifelse(A$n_cant_do > 0, sprintf("continuous: 'cannot do' set to missing (n=%d); consider unable-as-outcome sensitivity", A$n_cant_do), "")))

OUT <- A[, c("variable", "label", "health_class", "type", "n_levels", "levels", "target_pop", "n_applicable", "phase", "overlap_ksaa", "validity_applicable", "validity_cohort",
             "n_valid", "min_level_n", "n_cant_do", "size_ok", "validity_ok", "weight_basis", "in_pool", "status", "reason", "sparse_level", "dup_family", "threshold_version", "summarise_to", "s8_version", "note_cant_do", "section", "subsection")]
OUT$audit_version <- "v1.2 2026-10-06 (v1.1 + fixes to rule execution: ache word boundary, repeated-measure grouping and summaries, multi-item test scores, HAE advice items, Table S8-harmonized binary recoding (each variable once); no effect size/P consulted)"
write.csv(OUT, file.path(ew, "ewas_audit_table.csv"), row.names = FALSE)

P <- OUT[OUT$in_pool, ]
cat(sprintf("=== 825 池(实际 %d)按状态 × 健康类别 ===\n", nrow(P))); print(table(P$health_class, P$status))
cat("\n=== 池内状态合计 ===\n"); print(table(P$status))
cat("\n=== 池外(500 ≤ n_valid < 1000、规则内)按状态 ===\n"); O <- OUT[!OUT$in_pool & OUT$n_valid >= 500, ]; print(table(O$status)); cat("其中会进主分析的:", sum(O$status == "main"), "\n")
cat("\n=== 池内 supplementary 的原因 ===\n"); print(sort(table(P$reason[P$status == "supplementary"]), decreasing = TRUE))
cat("\n=== 池内待复核(review)===\n"); print(table(P$health_class[P$status == "review"]))
cat("\n=== 池内有效率 < 70%(在适用人群内)的主类变量前 15 ===\n"); x <- P[P$health_class %in% c("objective","disease","symptom","function","self_rated") & !P$validity_ok, ]; x <- x[order(x$validity_applicable), ]
print(head(x[, c("variable", "label", "n_applicable", "n_valid", "validity_applicable")], 15), row.names = FALSE)
cat("\n=== 'cannot do' 被清洗成缺失的连续指标 ===\n"); print(P[P$note_cant_do != "", c("variable", "label", "note_cant_do")], row.names = FALSE)
cat("\n=== 重复测量族 ===\n"); print(unique(P$dup_family[P$dup_family != ""]))
cat("\n=== 相位特异 ===\n"); print(table(P$phase))
cat("\n=== 池内 unclassified ===\n"); print(P[P$health_class == "unclassified", c("variable","label","subsection")], row.names = FALSE)
cat("\n=== 需改缺失的 can't grade 变量 ===\n"); print(P[grepl("^RECODE", P$note_cant_do), c("variable","label","n_cant_do")], row.names = FALSE)
cat("\n=== 主分析内主类 × 类型 ===\n"); print(table(P$health_class[P$status=="main"], P$type[P$status=="main"]))
cat("\n=== 池内 review 明细 ===\n"); print(P[P$status == "review" & P$health_class != "unclassified", c("variable","label","health_class")], row.names = FALSE)
cat("\n=== 池内 covariate 明细 ===\n"); print(P[P$status == "covariate", c("variable","label")], row.names = FALSE)
cat("\n=== 池内 exposure_biomarker ===\n"); print(P[P$health_class == "exposure_biomarker", c("variable","label")], row.names = FALSE)
