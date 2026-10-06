library(GSVA)
source("/proj/c.zihao/work3/function/ssGSEA.R")

LABELS   <- c("frac0.20", "frac0.40", "frac0.60", "frac0.80", "frac1.00")
EXPR_DIR <- "/proj/c.zihao/work3/05independence/01data"
CORUM    <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE <- "/proj/c.zihao/work3/05independence/results/bench/ssgsea"

if (dir.exists(OUT_BASE)) unlink(OUT_BASE, recursive = TRUE, force = TRUE)
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

for (label in LABELS) {
  cat(sprintf("[%s] ssgsea ...\n", label))
  expr_path <- file.path(EXPR_DIR, paste0("KIRC_", label, "_exprSet.csv"))
  out_dir   <- file.path(OUT_BASE, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  score_ssgsea(expr_path, CORUM, out_dir)
}

cat("Done.\n")
