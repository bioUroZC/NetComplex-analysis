
library(GSVA)
source("/proj/c.zihao/work3/function/ssGSEA.R")

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
STABLE_DIR <- "/proj/c.zihao/work3/01protein/02cor/01stable"
CORUM <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE <- "/proj/c.zihao/work3/01protein/02cor/02scored/ssgsea"

for (cancer in CANCERS) {
  cat(sprintf("[%s] protein-ssgsea ...\n", cancer))
  expr_path <- file.path(STABLE_DIR, paste0(cancer, "_protein_stable.csv"))
  out_dir <- file.path(OUT_BASE, cancer)
  score_ssgsea(expr_path, CORUM, out_dir)
}

cat("Done.\n")
