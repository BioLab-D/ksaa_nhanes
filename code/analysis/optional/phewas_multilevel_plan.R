# 多分类变量(>2 级)的处理计划(作结局时只能连续或二分类),phewas.R 读取其 class 列。
# 四类:ORD = 真有序量表 → 整数评分作连续结局(Gaussian);LAT = 无/右/左/双侧 → 任一侧 vs 无(二分类);
#       YNA = 是/否/不适用 → 第三档设缺失后 是 vs 否;NOM = 真无序 → 不能进 CLARITE,只能弃用或另立对照。
# 另标 drop_levels:需先设为缺失的档("Can't grade""Unable to evaluate""Doesn't work"等)。输出 results/analysis/phewas/ewas_multilevel_plan.csv
F <- read.csv(here::here("results", "analysis", "phewas", "ewas_results.csv"))
M <- F[F$type == "categorical" & F$levels > 2, c("variable", "label", "levels", "domain", "group")]
LAT <- c("HAG23", "PEP10B1", "PEP10B2", "PEP10B3", "PEP10B4", "PEP10C1", "PEP10C2", "PEP4A1", "PEP4A3", "PEP4CMCJ", "PEP4CMCP", "PEP4CMCT",
         "PEP4E1A", "PEP4E1B", "PEP4E1C", "PEP4E1D", "PEP4E1E", "PFPLSER", "PFPRSER")          # 最后两个:完全/部分/不能 → 任何受限 vs 完全
YNA <- c("HAL14D", "HAP15", "HAF2", "HAK4", "HAL17", "PFPUSTND", "PFPSCOOT")          # 第三档不适用/不能 → 缺失
NOM <- c("HAS22", "HAV12R", "HFA6XCR", "HAV10R", "HFE1", "HFE10", "HFE12", "HAN1D", "HAS17", "HAS17R", "HAS19R", "HAS28", "HAS29", "HAP2", "HFB3",
         "GUPFDX1R", "HFA12", "DEPLPTYP", "DEPUPTYP", "DEPDQL1", "DEPDQU1", "ANP", "MIP", "HZP")  # 血涂片分级只有 3 档且标签不清,保守归 NOM
M$class <- ifelse(M$variable %in% LAT, "LAT", ifelse(M$variable %in% YNA, "YNA", ifelse(M$variable %in% NOM, "NOM", "ORD")))
M$drop_levels <- ""
M$drop_levels[grepl("^FPP", M$variable)] <- "Can't grade"
M$drop_levels[M$variable %in% c("DEPLPTYP", "DEPUPTYP")] <- "Unable to evaluate"
M$drop_levels[M$variable == "HAQ1"] <- "Has no natural teeth"
M$drop_levels[M$variable %in% c("HAL14D")] <- "Doesn't work / not applicable"
M$handling <- c(ORD = "integer score as continuous outcome (Gaussian); reverse code so higher = worse where needed",
                LAT = "any side vs none (binary)", YNA = "yes vs no after setting the third level to missing (binary)",
                NOM = "cannot be an outcome in CLARITE: omit, or one pre-specified contrast")[M$class]
write.csv(M, here::here("results", "analysis", "phewas", "ewas_multilevel_plan.csv"), row.names = FALSE)
cat(sprintf("多分类 %d 个:", nrow(M))); print(table(M$class))
cat("\n可进 CLARITE 的:ORD + LAT + YNA =", sum(M$class != "NOM"), ";不能的 NOM =", sum(M$class == "NOM"), "\n")
cat("\nNOM 清单:", paste(M$variable[M$class == "NOM"], collapse = " "), "\n")
cat("\n全部 825 在 CLARITE 下的结局类型:连续 %d(原连续 291 + ORD %d);二分类 %d(原 412 + LAT %d + YNA %d);弃用 %d\n" |>
      sprintf(291 + sum(M$class == "ORD"), sum(M$class == "ORD"), 412 + sum(M$class %in% c("LAT", "YNA")), sum(M$class == "LAT"), sum(M$class == "YNA"), sum(M$class == "NOM")))
