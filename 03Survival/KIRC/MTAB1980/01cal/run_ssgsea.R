library(GSVA)
source("/proj/c.zihao/work3/function/ssGSEA.R")

EXPR <- "/proj/c.zihao/work3/03Survival/KIRC/MTAB1980/data/exprSet.csv"
CORUM <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT   <- "/proj/c.zihao/work3/03Survival/KIRC/MTAB1980/bench/ssgsea"

if (dir.exists(OUT)) unlink(OUT, recursive = TRUE, force = TRUE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[MTAB1980] ssgsea ...\n"))
score_ssgsea(EXPR, CORUM, OUT)
cat("Done.\n")
