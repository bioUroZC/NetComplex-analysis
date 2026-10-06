library(UCell)
library(data.table)
source("/proj/c.zihao/work3/function/UCell.R")

SAMPLES  <- c("GSM2560248", "GSM2560249")
EXPR_DIR <- "/proj/c.zihao/work3/00data/singleCell/GSE96583/out"
CORUM    <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_BASE <- "/proj/c.zihao/work3/06singleCell/GSE96583/bench/ucell"

if (dir.exists(OUT_BASE)) unlink(OUT_BASE, recursive = TRUE, force = TRUE)
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

corum <- read.table(CORUM, sep = ",", header = TRUE, stringsAsFactors = FALSE)

for (sample in SAMPLES) {
  cat(sprintf("[%s] ucell ...\n", sample))
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

  scores <- round(run_ucell(expr, corum), 5)
  write.csv(scores, file.path(out_dir, "ucell.csv"), quote = FALSE)
  cat(sprintf("Saved ucell: %d x %d\n", nrow(scores), ncol(scores)))
}

cat("Done.\n")
