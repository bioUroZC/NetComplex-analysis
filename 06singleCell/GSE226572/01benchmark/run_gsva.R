library(GSVA)
library(data.table)
source("/proj/c.zihao/work3/function/GSVA.R")

SAMPLES  <- c("GSM7079694", "GSM7079695", "GSM7079696", "GSM7079702", "GSM7079703", "GSM7079705", "GSM7079710", "GSM7079711", "GSM7079713")
EXPR_DIR <- "/proj/c.zihao/work3/00data/singleCell/GSE226572/out"
CORUM    <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE <- "/proj/c.zihao/work3/06singleCell/GSE226572/bench/gsva"

if (dir.exists(OUT_BASE)) unlink(OUT_BASE, recursive = TRUE, force = TRUE)
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

corum <- read.table(CORUM, sep = ",", header = TRUE, stringsAsFactors = FALSE)

for (sample in SAMPLES) {
  cat(sprintf("[%s] gsva ...\n", sample))
  expr_path <- file.path(EXPR_DIR, sample, "logNorm.txt")
  out_dir   <- file.path(OUT_BASE, sample)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  data <- as.data.frame(fread(expr_path))
  names(data)[1] <- "gene"
  rownames(data) <- data$gene
  data$gene <- NULL
  expr <- as.matrix(data)
  mode(expr) <- "numeric"
  cat("expr dim:", dim(expr), "\n")

  scores <- round(run_gsva(expr, corum), 5)
  write.csv(scores, file.path(out_dir, "gsva.csv"), quote = FALSE)
  cat(sprintf("Saved gsva: %d x %d\n", nrow(scores), ncol(scores)))
}

cat("Done.\n")
