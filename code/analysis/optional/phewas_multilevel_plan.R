# Handling plan for multi-level categorical variables (>2 levels; as outcomes they must be continuous or binary); phewas.R reads its class column.
# Four classes: ORD = truly ordinal scale -> integer score as continuous outcome (Gaussian); LAT = none/right/left/both -> either side vs none (binary);
#       YNA = yes/no/not applicable -> yes vs no after setting the third level to missing; NOM = truly nominal -> cannot enter CLARITE; omit or define a separate contrast.
# drop_levels marks levels to set to missing first ("Can't grade", "Unable to evaluate", "Doesn't work", etc.). Output: results/analysis/phewas/ewas_multilevel_plan.csv
F <- read.csv(here::here("results", "analysis", "phewas", "ewas_results.csv"))
M <- F[F$type == "categorical" & F$levels > 2, c("variable", "label", "levels", "domain", "group")]
LAT <- c("HAG23", "PEP10B1", "PEP10B2", "PEP10B3", "PEP10B4", "PEP10C1", "PEP10C2", "PEP4A1", "PEP4A3", "PEP4CMCJ", "PEP4CMCP", "PEP4CMCT",
         "PEP4E1A", "PEP4E1B", "PEP4E1C", "PEP4E1D", "PEP4E1E", "PFPLSER", "PFPRSER")          # last two: full/partial/unable -> any limitation vs full
YNA <- c("HAL14D", "HAP15", "HAF2", "HAK4", "HAL17", "PFPUSTND", "PFPSCOOT")          # third level (not applicable/unable) -> missing
NOM <- c("HAS22", "HAV12R", "HFA6XCR", "HAV10R", "HFE1", "HFE10", "HFE12", "HAN1D", "HAS17", "HAS17R", "HAS19R", "HAS28", "HAS29", "HAP2", "HFB3",
         "GUPFDX1R", "HFA12", "DEPLPTYP", "DEPUPTYP", "DEPDQL1", "DEPDQU1", "ANP", "MIP", "HZP")  # blood smear grades have only 3 levels with unclear labels; conservatively NOM
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
cat(sprintf("Multi-level variables: %d: ", nrow(M))); print(table(M$class))
cat("\nUsable in CLARITE: ORD + LAT + YNA =", sum(M$class != "NOM"), "; not usable (NOM) =", sum(M$class == "NOM"), "\n")
cat("\nNOM list:", paste(M$variable[M$class == "NOM"], collapse = " "), "\n")
cat("\nOutcome types of all 825 under CLARITE: continuous %d (291 originally continuous + ORD %d); binary %d (412 original + LAT %d + YNA %d); omitted %d\n" |>
      sprintf(291 + sum(M$class == "ORD"), sum(M$class == "ORD"), 412 + sum(M$class %in% c("LAT", "YNA")), sum(M$class == "LAT"), sum(M$class == "YNA"), sum(M$class == "NOM")))
