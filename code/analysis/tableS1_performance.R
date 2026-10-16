# Table S1: age-prediction performance (mean +/- SD over folds) of the five backbones and the ensemble; also by subgroup (not reported).
# Input: results/ksaa_oof.csv, results/predictions/m4_*_oof.csv, results/ksaa_ensemble_log.txt (cross-check)
# Output: results/tableS1_performance.xlsx (sheets keep the draft names "Table 2" / "Table S2"), results/work/tableS1*_en.csv
res <- here::here("results")
sp  <- here::here("results", "work")
suppressMessages(library(openxlsx))
B <- read.csv(file.path(res,"ksaa_oof.csv")); age <- B$age; fold <- B$fold; K <- 0:4
sl <- read.csv(file.path(res, "work", "surv_labels.csv")); ii <- match(B$SEQN, sl$SEQN); B$sex <- sl$sex[ii]; B$race <- sl$race[ii]; B$phase <- sl$phase[ii]
al <- function(f){ x <- read.csv(f); v <- x$pred_age[match(B$SEQN, x$SEQN)]; stopifnot(!anyNA(v)); v }
v2 <- function(a) al(file.path(res, "predictions", sprintf("m4_%s_oof.csv", a)))
M <- list(
  list(l="ResNet-18", p=v2("resnet18")), list(l="ResNet-50", p=v2("resnet50")), list(l="EfficientNet-B0", p=v2("efficientnetb0")),
  list(l="ConvNeXt-Tiny", p=v2("convnexttiny")), list(l="DenseNet-121", p=v2("densenet121")))
M <- c(M, list(
  list(l="Weighted ensemble{a}", p=B$pred_raw, bold=TRUE),
  list(l="Weighted ensemble, bias-corrected{b}", p=B$pred_corr, bold=TRUE)))
N <- format(nrow(B), big.mark=",")
f3 <- function(v, dg) if(all(is.na(v))) "—" else sprintf(paste0("%.",dg,"f ± %.",dg,"f"), mean(v), sd(v))
met <- function(p){ m <- sapply(K, function(k){ j <- fold==k; e <- p[j]-age[j]; ok <- sd(p[j])>0
    c(mean(abs(e)), 100*mean(abs(e)<=5), if(ok) cor(p[j],age[j]) else NA, if(ok) cor(p[j],age[j],method="spearman") else NA) })
  c(f3(m[1,],3), f3(m[2,],1), f3(m[3,],3), f3(m[4,],3)) }
T2 <- do.call(rbind, lapply(M, function(x) if(!is.null(x$g)) data.frame(lab=x$g, mae="", cs="", r="", rs="", type="hdr", bold=FALSE)
                                           else { v <- met(x$p); data.frame(lab=x$l, mae=v[1], cs=v[2], r=v[3], rs=v[4], type="row", bold=isTRUE(x$bold)) }))
chk <- T2[T2$type=="row", c("lab","mae")]
# cross-check the ensemble MAE against ksaa_ensemble_log.txt
kas <- readLines(file.path(res,"ksaa_ensemble_log.txt")); ref <- function(lab) sub("^.*MAE ([0-9.]+ ± [0-9.]+).*$", "\\1", grep(lab, kas, value=TRUE, fixed=TRUE)[1])
cat(sprintf("Check (should equal ksaa_ensemble_log.txt: raw %s / corrected %s):\n", ref("Ensemble (raw)"), ref("Ensemble (bias-corrected)"))); print(tail(chk, 7), row.names=FALSE)
stopifnot(tail(chk$mae,1)==ref("Ensemble (bias-corrected)"), chk$mae[nrow(chk)-1]==ref("Ensemble (raw)"))
T2_notes <- c(
  paste0("The five backbones were pre-specified following the Mammo-AGE protocol. Both knees of each participant were entered jointly and fused by the Mammo-AGE instance-bag transformer (448 × 448 pixels), trained with a probabilistic ordinal embedding head and mean–variance loss. Metrics were computed in each held-out fold of participant-level five-fold cross-fitting and are reported as mean ± SD across folds (n = ", N, ")."),
  "Abbreviations: CS (α = 5 y), cumulative score, the percentage of participants predicted within 5 years of chronological age; MAE, mean absolute error; r, Pearson correlation coefficient; rs, Spearman correlation coefficient.",
  "{a}Backbone predictions weighted in proportion to inverse MAE on an internal validation split within the training folds.",
  "{b}Linear age-bias correction fitted in the training folds and applied to the held-out fold; the primary knee age measure.")
# ---- subgroups (pooled out-of-fold predictions of the ensemble) ----
fd <- function(x,dg){ z <- formatC(x, format="f", digits=dg); z <- sub("^-(0(\\.0+)?)$", "\\1", z); sub("^-", "−", z) }
gap <- B$pred_corr - age
G <- list(list("Overall", rep(TRUE,nrow(B))),
  list(h="Survey phase"), list("Phase 1 (1988–1991)", B$phase==1, 1), list("Phase 2 (1991–1994)", B$phase==2, 1),
  list(h="Sex"), list("Male", B$sex==1, 1), list("Female", B$sex==2, 1),
  list(h="Race and ethnicity"), list("Non-Hispanic White", B$race==1, 1), list("Non-Hispanic Black", B$race==2, 1), list("Mexican American", B$race==3, 1), list("Other", B$race==4, 1),
  list(h="Age group, years"))
ab <- cut(age, c(59,64,69,74,79,84,120), labels=c("60–64","65–69","70–74","75–79","80–84","85–90{a}"))
for(lv in levels(ab)) G[[length(G)+1]] <- list(lv, ab==lv, 1)
S2 <- do.call(rbind, lapply(G, function(x){ if(!is.null(x$h)) return(data.frame(lab=x$h, n="", mae_raw="", err_raw="", mae_corr="", cs="", gap="", ind=0, type="hdr"))
  j <- x[[2]]; data.frame(lab=x[[1]], n=formatC(sum(j), format="d", big.mark=","), mae_raw=fd(mean(abs(B$pred_raw[j]-age[j])),2), err_raw=fd(mean(B$pred_raw[j]-age[j]),2), mae_corr=fd(mean(abs(gap[j])),2),
    cs=fd(100*mean(abs(gap[j])<=5),1), gap=sprintf("%s (%s)", fd(mean(gap[j]),2), fd(sd(gap[j]),2)), ind=if(length(x)>2) x[[3]] else 0, type="row") }))
S2_notes <- c(
  paste0("Pooled out-of-fold predictions of the ensemble (n = ", N, "). Mean error, uncorrected: predicted minus chronological age (positive, overestimation)."),
  "Abbreviations: KSAA, knee skeletal age acceleration; MAE, mean absolute error.",
  "{a}Age was top-coded at 90 years in NHANES III.")
names_s2 <- c("Subgroup","n","MAE, uncorrected, years","Mean error, uncorrected, years","KSAA, mean (SD), years")
write.csv(T2, file.path(sp,"tableS1_en.csv"), row.names=FALSE, fileEncoding="UTF-8"); write.csv(S2, file.path(sp,"tableS1_subgroups_en.csv"), row.names=FALSE, fileEncoding="UTF-8")
writeLines(T2_notes, file.path(sp,"tableS1_notes.txt"), useBytes=TRUE); writeLines(S2_notes, file.path(sp,"tableS1_subgroups_notes.txt"), useBytes=TRUE)
# ---- workbook ----
wb <- createWorkbook(); modifyBaseFont(wb, fontSize=10, fontName="Times New Roman")
sheet <- function(name, title, out, hdr_rows, ind_rows, bold_rows, widths, notes){
  addWorksheet(wb, name, gridLines=FALSE); nc <- ncol(out); last <- nrow(out)+3
  writeData(wb, name, title, startRow=1, startCol=1); addStyle(wb, name, createStyle(textDecoration="bold"), rows=1, cols=1)
  writeData(wb, name, out, startRow=3, startCol=1, headerStyle=createStyle(textDecoration="bold", halign="center", valign="center", wrapText=TRUE))
  addStyle(wb, name, createStyle(halign="left", textDecoration="bold"), rows=3, cols=1, stack=TRUE)
  addStyle(wb, name, createStyle(border="top", borderStyle="medium"), rows=3, cols=1:nc, gridExpand=TRUE, stack=TRUE)
  addStyle(wb, name, createStyle(border="bottom", borderStyle="thin"), rows=3, cols=1:nc, gridExpand=TRUE, stack=TRUE)
  addStyle(wb, name, createStyle(halign="center"), rows=4:last, cols=2:nc, gridExpand=TRUE, stack=TRUE)
  if(length(ind_rows))  addStyle(wb, name, createStyle(indent=2), rows=ind_rows+3, cols=1, stack=TRUE)
  if(length(hdr_rows))  addStyle(wb, name, createStyle(textDecoration="italic"), rows=hdr_rows+3, cols=1, stack=TRUE)
  if(length(bold_rows)) addStyle(wb, name, createStyle(textDecoration="bold"), rows=bold_rows+3, cols=1:nc, gridExpand=TRUE, stack=TRUE)
  addStyle(wb, name, createStyle(border="bottom", borderStyle="medium"), rows=last, cols=1:nc, gridExpand=TRUE, stack=TRUE)
  setColWidths(wb, name, cols=1:nc, widths=widths); setRowHeights(wb, name, rows=3, heights=30)
  writeData(wb, name, notes, startRow=last+1, startCol=1); addStyle(wb, name, createStyle(fontSize=9), rows=(last+1):(last+length(notes)), cols=1, stack=TRUE) }
o2 <- T2[,c("lab","mae","cs","r","rs")]; names(o2) <- c("Model","MAE, years","CS (α = 5 y), %","Pearson r","Spearman rs")
sheet("Table 2", "Table 2. Age prediction performance of the five pre-specified backbones and their weighted ensemble", o2,
      integer(0), integer(0), which(T2$bold), c(40,16,16,16,16), T2_notes)
os <- S2[,c("lab","n","mae_raw","err_raw","gap")]; names(os) <- names_s2
sheet("Table S2", "Supplementary Table S2. Performance of the multi-view ensemble by participant subgroup", os,
      which(S2$type=="hdr"), which(S2$ind==1), integer(0), c(30,9,16,16,20), S2_notes)
f <- file.path(res,"tableS1_performance.xlsx"); saveWorkbook(wb, f, overwrite=TRUE)
td <- tempfile(); dir.create(td); unzip(f, exdir=td)
for(rf in list.files(file.path(td,"xl/worksheets/_rels"), full.names=TRUE)){ x <- paste(readLines(rf, warn=FALSE), collapse="")
  writeLines(gsub('<Relationship [^>]*Target="\\.\\./drawings/[^"]*"/>', '', x), rf, useBytes=TRUE) }
for(sf in list.files(file.path(td,"xl/worksheets"), pattern="^sheet[0-9]+\\.xml$", full.names=TRUE)){ y <- paste(readLines(sf, warn=FALSE, encoding="UTF-8"), collapse="")
  y <- gsub('<drawing r:id="[^"]*"/>', '', y); writeLines(gsub('<legacyDrawing r:id="[^"]*"/>', '', y), sf, useBytes=TRUE) }
ssf <- file.path(td,"xl","sharedStrings.xml"); ss <- paste(readLines(ssf, warn=FALSE, encoding="UTF-8"), collapse="")
blocks <- regmatches(ss, gregexpr("<si>.*?</si>", ss, perl=TRUE))[[1]]; nconv <- 0
for(b in unique(blocks[grepl("\\{[a-j]\\}", blocks, perl=TRUE)])){
  txt <- sub("^<si><t[^>]*>(.*)</t></si>$", "\\1", b, perl=TRUE); sz <- if(grepl("^\\{[a-j]\\}", txt, perl=TRUE)) 9 else 10
  parts <- regmatches(txt, gregexpr("\\{[a-j]\\}|[^{]+|\\{", txt, perl=TRUE))[[1]]
  run <- function(s, sup) sprintf('<r><rPr>%s<sz val="%d"/><rFont val="Times New Roman"/></rPr><t xml:space="preserve">%s</t></r>', if(sup) '<vertAlign val="superscript"/>' else '', sz, s)
  nb <- paste0("<si>", paste(vapply(parts, function(p) if(grepl("^\\{[a-j]\\}$",p,perl=TRUE)) run(substr(p,2,2),TRUE) else run(p,FALSE), ""), collapse=""), "</si>")
  ss <- gsub(b, nb, ss, fixed=TRUE); nconv <- nconv+1 }
writeLines(ss, ssf, useBytes=TRUE)
owd <- setwd(td); unlink(f); zip(f, files=list.files(".", recursive=TRUE, all.files=TRUE), flags="-X -q"); setwd(owd)
cat(sprintf("-> %s | Table 2 %d rows, S2 %d rows, superscript strings %d\n", f, nrow(T2), nrow(S2), nconv))
