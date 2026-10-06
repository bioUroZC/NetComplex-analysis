set.seed(2024)

CELL_DIR <- "/proj/c.zihao/work3/08runtime/cell"
DATA_DIR <- file.path(CELL_DIR, "data")

for (d in c(DATA_DIR, file.path(CELL_DIR, "results"))) {
  if (dir.exists(d)) { unlink(d, recursive = TRUE); cat("Removed:", d, "\n") }
}
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

STRING_PATH      <- "/proj/c.zihao/work3/00data/string/links.csv"
CELL_COUNTS      <- c(500, 1000, 2000, 5000, 10000)
COMPLEX_COUNTS   <- c(100, 200, 500, 1000, 2000, 5000)
COMPLEX_SIZE_MIN <- 2
COMPLEX_SIZE_MAX <- 6
DROPOUT_RATE     <- 0.7



cat("Loading STRING...\n")
huri <- read.csv(STRING_PATH, row.names = 1)
keep <- intersect(c("protein1", "protein2", "score"), names(huri))
huri <- huri[, keep, drop = FALSE]
genes <- sort(unique(c(as.character(huri$protein1), as.character(huri$protein2))))
N     <- length(genes)
cat(sprintf("Gene pool: %d genes, %d edges\n", N, nrow(huri)))

write.csv(huri, file.path(DATA_DIR, "links.csv"), row.names = FALSE)
cat("Saved links.csv\n")



make_expr <- function(n_cells) {
  mat <- matrix(
    pmax(rnorm(N * n_cells, mean = 2.0, sd = 1.5), 0),
    nrow = N, ncol = n_cells,
    dimnames = list(genes, paste0("Cell_", seq_len(n_cells)))
  )
  dropout <- matrix(runif(N * n_cells) > DROPOUT_RATE, nrow = N)
  mat * dropout
}

make_complexes <- function(n) {
  data.frame(
    Complex = paste0("Complex_", seq_len(n)),
    Genes   = vapply(seq_len(n), function(i) {
      k <- sample(COMPLEX_SIZE_MIN:COMPLEX_SIZE_MAX, 1)
      paste(sample(genes, k), collapse = ";")
    }, character(1)),
    stringsAsFactors = FALSE
  )
}



cat("Cell-count sweep...\n")
for (n in CELL_COUNTS) {
  set.seed(10 + n)
  write.csv(round(make_expr(n), 5),
            file.path(DATA_DIR, sprintf("expr_n%d.csv", n)))
  cat(sprintf("  expr_n%d.csv\n", n))
}



cat("Complex-count sweep...\n")
for (n in COMPLEX_COUNTS) {
  set.seed(20 + n)
  write.csv(make_complexes(n),
            file.path(DATA_DIR, sprintf("complexes_n%d.csv", n)), row.names = FALSE)
  cat(sprintf("  complexes_n%d.csv\n", n))
}

cat("Done.\n")
