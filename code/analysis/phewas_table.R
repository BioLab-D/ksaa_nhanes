# PheWAS 补充表:tableS8(表 S8,进入主扫描的变量审核表)、tableS9(表 S9,主分析全部表型)
# 输入 results/analysis/phewas/{phewas_results.csv, ewas_audit_table.csv};输出 results/analysis/tables/{tableS8,tableS9}.html(样式同 manuscript_tables.R)
res <- here::here("results"); v3 <- file.path(res, "analysis"); td <- file.path(v3, "tables"); ew <- file.path(v3, "phewas")
source(here::here("code", "analysis", "phewas_labels.R"))
mi <- function(x) gsub("-", "−", x, fixed = TRUE)
f <- function(x, d = 2) ifelse(is.na(x), "", mi(sub("^-(0\\.0+)$", "\\1", formatC(round(x + sign(x) * 1e-9, d), format = "f", digits = d))))
# 分隔符按显示值选(同 manuscript_tables.R 的 ciD):任一端显示为负时用 " to ",否则用 "–"
ci <- function(h, l, u, d = 2) { lo <- f(l, d); hi <- f(u, d); ifelse(is.na(h), "—", sprintf("%s (%s%s%s)", f(h, d), lo, ifelse(grepl("^−", lo) | grepl("^−", hi), " to ", "–"), hi)) }
fpe <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, mi(sub("e", "×10<sup>", sprintf("%.1e", p)) |> paste0("</sup>") |> sub(pattern = "<sup>-0?", replacement = "<sup>−")), formatC(p, format = "f", digits = 3)))
cm <- function(x) formatC(round(x), format = "d", big.mark = ",")
th <- function(h, nl = 2) paste0("<thead><tr>", paste0(sprintf("<th%s>%s</th>", ifelse(seq_along(h) <= nl, "", ' class="num"'), h), collapse = ""), "</tr></thead>")
row <- function(cells, nl = 2, ind = FALSE) paste0("<tr>", paste0(sprintf('<td class="%s">%s</td>', c(if (ind) "ind" else "", rep("", nl - 1), rep("num", length(cells) - nl)), cells), collapse = ""), "</tr>")
hdr <- function(txt, k) sprintf('<tr class="hdr"><td colspan="%d">%s</td></tr>', k, txt)
esc <- function(x) gsub("<(?=[0-9=\\s])", "&lt;", x, perl = TRUE)
wrap <- function(head_, body, notes, nl = 2) esc(paste0('<div class="twrap"><table>', th(head_, nl), "<tbody>", paste0(body, collapse = ""), "</tbody></table></div>",
                                                       '<div class="tnote">', paste0("<p>", notes, "</p>", collapse = ""), "</div>"))
R <- read.csv(file.path(ew, "phewas_results.csv"), stringsAsFactors = FALSE)
A <- read.csv(file.path(ew, "ewas_audit_table.csv"), stringsAsFactors = FALSE)
MOD <- c(gaussian = "Linear", binomial = "Logistic", ordinal = "Ordinal logistic", nominal = "Joint test",
         `binomial (ordinal fallback)` = "Logistic (ordinal fit singular)")   # 写实际拟合的模型
est <- function(r) if (r$model == "gaussian") ci(r$est, r$lo, r$hi) else if (r$model == "nominal") "—" else ci(r$est, r$lo, r$hi)
unit <- function(r) if (r$model == "gaussian") "SD" else if (r$model == "nominal") "" else "OR"
pcell <- function(v) ifelse(is.na(v), "—", ifelse(v < 0.05, sprintf("<b>%s</b>", fpe(v)), fpe(v)))
rep_ <- function(r) if (is.na(r$p_phase1) || is.na(r$p_phase2)) "n/a" else if (r$replicated) "Yes" else "No"
# 合成表型显示其组成条目的 CDC 代码,而非内部变量名
CODE <- c(MAPB_story_recall_immediate_total = "MAPB2A–F", MAPD_story_recall_delayed_total = "MAPD2A–F",
          PFP_gait_speed = "PFPTWLKA–B", PFP_steps_8ft_mean = "PFPNSTPA–B", pulse_mean = "HAZA5R, PEP6DR", HAP18_serial_subtraction = "HAP18A–E",
          HAP17_word_registration = "HAP17A1–3", HAP19_word_recall = "HAP19A–C")
cd <- function(v) ifelse(v %in% names(CODE), CODE[v], v)
one <- function(r, rank = NA) row(c(sprintf("%s%s <span class=\"muted\">(%s)</span>", if (is.na(rank)) "" else sprintf("%d. ", rank), nice(r$variable, r$phenotype), cd(r$variable)),
                                     clsn(r$health_class), MOD[[if (!is.null(r$model_used) && !is.na(r$model_used)) r$model_used else r$model]], cm(r$n), est(r), unit(r), fpe(r$p), fpe(r$fdr),
                                     f(r$est_phase1), f(r$est_phase2), rep_(r),
                                      if (!is.na(r$clump_lead) && r$clump_lead != r$variable) sprintf("→ %s", cd(r$clump_lead)) else ""), nl = 3, ind = TRUE)
clsn <- function(h) if (h %in% names(CLS_NAME)) CLS_NAME[[h]] else h
HEAD <- c("Phenotype", "Health domain", "Model", "n", "Estimate (95% CI)", "Unit", "P", "FDR", "Phase 1 estimate", "Phase 2 estimate", "Replicated", "Clumped with")
K <- length(HEAD)

# ---- 主分析:按健康类别分组,组内按 P ----
M <- R[R$set == "main", ]; M <- M[order(factor(M$health_class, CLS_ORDER), M$p), ]
body <- c(); rk <- 0
for (cl in CLS_ORDER) { x <- M[M$health_class == cl, ]; if (!nrow(x)) next
  body <- c(body, hdr(sprintf("%s (%d phenotypes; %d at FDR &lt;0.05)", CLS_NAME[[cl]], nrow(x), sum(x$fdr < 0.05, na.rm = TRUE)), K))
  for (i in seq_len(nrow(x))) body <- c(body, one(x[i, ])) }
L <- M[!is.na(M$clump_lead) & M$clump_lead == M$variable, ]
NFB <- sum(M$model_used == "binomial (ordinal fallback)", na.rm = TRUE)
WORDS <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten")
notes_m <- c(
  sprintf("n: participants with the phenotype and complete covariates. Each of %s audited phenotypes (Table S8) was regressed on KSAA (per SD, 3.95 years) with NHANES III survey weights and design, adjusted for age, sex, race and ethnicity, education and poverty–income ratio: weighted linear regression for standardized continuous phenotypes (estimate, SD per SD of KSAA), logistic regression for binary and proportional-odds regression for ordinal phenotypes (odds ratio per SD). Further details are given in the Supplementary Methods. Items collected in only one survey phase were analyzed with that phase's design and weights, and cannot be replicated across phases (shown as n/a). %s ordinal phenotypes whose proportional-odds fit was singular were refitted as any level above the lowest versus the lowest and are labeled accordingly. Confidence intervals use the t distribution with the design degrees of freedom, as do the P values. Repeated measurements and multi-item tests were summarized to one phenotype, and difficulty with daily activities (any difficulty versus none), falls (at least one in the past 12 months) and the tandem stand (unable to hold for 10 s) were analyzed as binary outcomes defined as in Table S10 (Table S8).",
          cm(nrow(M)), WORDS[NFB]),
  sprintf("FDR, Benjamini–Hochberg across the %s tests; %s phenotypes had FDR &lt;0.05, %s passed Bonferroni correction for the %s effective independent tests (Li–Ji) and %s the conventional correction for all %s tests (the dashed line in Figure 4A). Replicated: same direction with P &lt;0.05 in both survey phases, each analyzed with its own weights; n/a, item asked in only one survey phase. Clumped: correlated at |r| &gt;0.7 with a stronger signal (%s independent signals, %s replicated).",
          cm(nrow(M)), cm(sum(M$fdr < 0.05, na.rm = TRUE)), cm(sum(M$p_liji < 0.05, na.rm = TRUE)), cm(round(M$meff[1])), cm(sum(M$bonf < 0.05, na.rm = TRUE)), cm(nrow(M)), cm(nrow(L)), cm(sum(L$replicated, na.rm = TRUE))),
  "Abbreviations: BMC, bone mineral content; BMD, bone mineral density; DIP, distal interphalangeal; FDR, false discovery rate; IP, interphalangeal; KSAA, knee skeletal age acceleration; MCP, metacarpophalangeal; MTP, metatarsophalangeal; NHANES III, Third National Health and Nutrition Examination Survey; OR, odds ratio; PIP, proximal interphalangeal.")
writeLines(wrap(HEAD, body, notes_m, 3), file.path(td, "tableS9.html"), useBytes = TRUE); cat("-> tableS9\n")

# ---- 审核表(只列进入主扫描的变量)----
stopifnot(all(A$overlap_ksaa[A$in_pool] == "no"))
P0 <- A[A$in_pool, ]; P <- P0[P0$status == "main", ]
P <- P[order(factor(P$health_class, CLS_ORDER), P$variable), ]
# 审核表:健康领域显示名、权重写变量名(WTPFEX6 / WTPFSD6),处理说明转为发表用语
hand_txt <- function(r) {
  out <- c(); rs <- r$reason; nt <- r$note_cant_do; st <- r$summarise_to; rec <- r$s8_version != ""
  S8T <- c(`binary: any difficulty vs none` = "Analyzed as any difficulty versus none, as in Table S10", `binary: at least one fall` = "Analyzed as at least one fall versus none, as in Table S10",
           `binary: unable to hold 10 s` = "Analyzed as inability to hold the tandem stand for 10 s, as in Table S10")
  if (rec) out <- c(out, S8T[[r$s8_version]])                       # 就地重编码为二分类:稀疏级合并与"unable to do"等级说明不再适用
  if (grepl("^sparse category", rs) && !rec) out <- c(out, "Categories with fewer than 30 participants merged with the adjacent category")
  if (grepl("story_recall", st)) out <- c(out, sprintf("Item summed into the %s story-recall total", if (grepl("immediate", st)) "immediate" else "delayed"))
  if (st == "PFP_gait_speed") out <- c(out, "Faster of the two trials converted to gait speed, as in Table S10")
  if (st == "PFP_steps_8ft_mean") out <- c(out, "Mean of the two trials analyzed as one phenotype")
  if (st == "pulse_mean") out <- c(out, "Mean of the household and examination readings analyzed as one phenotype")
  if (st == "HAP18_serial_subtraction") out <- c(out, "Item scored correct or incorrect (\u201ccan't do\u201d and \u201cdon't know\u201d incorrect) and summed into the serial-subtraction score")
  if (st == "HAP17_word_registration") out <- c(out, "Item summed into the three-word registration score")
  if (st == "HAP19_word_recall") out <- c(out, "Item summed into the three-word recall score")
  if (nzchar(rs) && !length(out)) stop("unmapped reason: ", rs)
  k <- length(out)
  if (grepl("cannot do' set to missing", nt) && !nzchar(st) && !rec) out <- c(out, "\u201cCannot do\u201d set to missing")
  if (grepl("unable to do' retained", nt) && !rec) out <- c(out, "\u201cUnable to do\u201d retained as the highest level")
  if (nzchar(nt) && length(out) == k && !nzchar(st) && !rec) stop("unmapped note: ", nt)
  paste(out, collapse = "; ") }
HA <- c("Variable", "Phenotype", "Health domain", "Type", "Applicable population", "Phase", "Validity, %", "n valid", "Smallest category", "Weight", "Handling")
body <- c()
for (i in seq_len(nrow(P))) { r <- P[i, ]
  typ <- if (r$type == "continuous") "continuous" else sprintf("%d levels", r$n_levels)
  j <- match(r$variable, R$variable)   # 实际所用权重(单期项目用该期权重,空腹指标用晨间亚样本权重);汇总为得分的单项按审核表
  wt <- if (!is.na(j)) sub(" .*$", "", R$weight[j]) else if (grepl("WTPFSD6", r$weight_basis)) "WTPFSD6" else if (grepl("^WTPFEX6", r$weight_basis)) "WTPFEX6" else stop("weight: ", r$weight_basis)
  body <- c(body, row(c(r$variable, tidy(r$label), clsn(r$health_class), typ, "all", r$phase,   # 研究样本均 ≥60 岁,codebook 的年龄适用域(如 "(2 months and over)")不适用
                       
                        f(100 * r$validity_applicable, 0), cm(r$n_valid), ifelse(is.na(r$min_level_n), "—", cm(r$min_level_n)), wt, hand_txt(r)), nl = 3)) }
nst <- table(factor(P0$status, c("main", "supplementary", "exposure_module", "covariate", "excluded"))); ex <- P0[P0$status == "excluded", ]
nx <- c(sum(grepl("single BP", ex$reason)), sum(grepl("prognosis", ex$reason)), sum(grepl("self-perception", ex$reason))); stopifnot(sum(nx) == nrow(ex))
notes_a <- c(sprintf("Audit by the pre-specified protocol, fixed before any association was examined (version 1.1); version 1.2 corrected the execution of the protocol after the initial scan (summary of repeated measurements and multi-item tests; classification of a food-frequency item and of physician-advice items) and recoded difficulty with daily activities, falls and the tandem stand as binary outcomes defined as in Table S10 (Supplementary Methods). Columns follow the protocol steps, and no candidate variable overlapped with KSAA. Of the %s candidate variables, the %s listed here entered the main scan (%s phenotypes after repeated measurements and multi-item tests were summarized); %s fell below the thresholds for sample size, category size or validity, %s were behaviors, diet, environmental or nutritional biomarkers, or psychosocial or reproductive factors, %s were covariates, and %s were excluded (%s single blood pressure readings replaced by the combined mean, %s physician prognosis items and %s self-perception or unclear observation items). Validity: proportion valid among participants to whom the item applied (blank-but-applicable, don't know and refused counted as invalid). Weight: WTPFEX6, examination weight; WTPFSD6, morning-session subsample weight, used for fasting assays.",
                     cm(nrow(P0)), cm(nst[["main"]]), cm(nrow(M)), cm(nst[["supplementary"]]), cm(nst[["exposure_module"]]), cm(nst[["covariate"]]), cm(nst[["excluded"]]), nx[1], nx[2], nx[3]),
             "Abbreviations: BMC, bone mineral content; BMD, bone mineral density; KSAA, knee skeletal age acceleration.")
writeLines(wrap(HA, body, notes_a, 3), file.path(td, "tableS8.html"), useBytes = TRUE); cat("-> tableS8\n")
