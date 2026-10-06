library(GSVA)
source("/proj/c.zihao/work3/function/ssGSEA.R")

CANCERS  <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
EXPR_DIR <- "/proj/c.zihao/work3/00data/cptacT/exprset"
CORUM    <- "/proj/c.zihao/work3/00data/complexDB/complexData.csv"
OUT_BASE <- "/proj/c.zihao/work3/04ComplexPortal/bench/ssgsea"

if (dir.exists(OUT_BASE)) unlink(OUT_BASE, recursive = TRUE, force = TRUE)
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

for (cancer in CANCERS) {
  cat(sprintf("[%s] ssgsea ...\n", cancer))
  expr_path <- file.path(EXPR_DIR, paste0(cancer, "_exprSet_filtered.csv"))
  out_dir   <- file.path(OUT_BASE, cancer)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  score_ssgsea(expr_path, CORUM, out_dir)
}

cat("Done.\n")
