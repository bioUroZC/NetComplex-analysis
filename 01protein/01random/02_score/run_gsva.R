library(GSVA)
source("/proj/c.zihao/work3/function/GSVA.R")

CANCERS  <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
EXPR_DIR <- "/proj/c.zihao/work3/00data/cptacT/exprset"
CORUM    <- Sys.getenv("CORUM_PATH")
OUT_BASE <- file.path(Sys.getenv("RAND_OUT_ROOT"), "gsva")
if (CORUM == "" || Sys.getenv("RAND_OUT_ROOT") == "") {
  stop("CORUM_PATH and RAND_OUT_ROOT must be set in the environment")
}

if (dir.exists(OUT_BASE)) unlink(OUT_BASE, recursive = TRUE, force = TRUE)
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

for (cancer in CANCERS) {
  cat(sprintf("[%s] gsva ...\n", cancer))
  expr_path <- file.path(EXPR_DIR, paste0(cancer, "_exprSet_filtered.csv"))
  out_dir   <- file.path(OUT_BASE, cancer)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  score_gsva(expr_path, CORUM, out_dir)
}

cat("Done.\n")
