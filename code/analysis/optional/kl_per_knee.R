# 每膝 K–L 分级(NHANES III xr.dat,左右膝)→ results/work/kl_per_knee.csv;放射学 OA = K–L ≥ 2
sp <- here::here("results", "work")
L  <- readLines(here::here("data", "raw", "nhanes3_xray", "xr.dat"), warn=FALSE)
cc <- function(a,b) as.numeric(substr(L,a,b))
d <- data.frame(SEQN=cc(1,5), R=cc(6,6), L=cc(7,7))
d$R[d$R %in% c(8,9)] <- NA          # 8=空缺但适用, 9=不可读
d$L[d$L %in% c(8,9)] <- NA
long <- rbind(data.frame(seqn=d$SEQN, side="R", kl=d$R),
              data.frame(seqn=d$SEQN, side="L", kl=d$L))
long <- long[!is.na(long$kl), ]
long$oa <- as.integer(long$kl >= 2)   # 放射学 OA
write.csv(long, file.path(sp,"kl_per_knee.csv"), row.names=FALSE)
cat("按膝的 K-L 标签:", nrow(long), "个膝\n")
print(table(long$kl))
cat("放射学OA(K-L>=2):", sum(long$oa), sprintf("(%.1f%%)\n", 100*mean(long$oa)))
