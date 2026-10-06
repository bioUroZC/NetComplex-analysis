library(UCell)
source("/proj/c.zihao/work3/function/UCell.R")

DATA_DIR       <- "/proj/c.zihao/work3/08runtime/cell/data"
RESULTS_DIR    <- "/proj/c.zihao/work3/08runtime/cell/results"
CELL_COUNTS    <- c(500, 1000, 2000, 5000, 10000)
COMPLEX_COUNTS <- c(100, 200, 500, 1000, 2000, 5000)

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
records <- list()

for (n in CELL_COUNTS) {
  cat(sprintf("[cell_count=%d] ucell ...\n", n))
  expr <- as.matrix(read.csv(file.path(DATA_DIR, sprintf("expr_n%d.csv", n)),
                             row.names = 1, check.names = FALSE))
  complexes <- read.csv(file.path(DATA_DIR, "complexes_n1000.csv"),
                        stringsAsFactors = FALSE)
  t0 <- proc.time()[["elapsed"]]
  run_ucell(expr, complexes)
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  %.3fs\n", elapsed))
  records[[length(records) + 1]] <- data.frame(
    sweep = "cell_count", value = n,
    n_cells = n, n_complexes = 1000,
    runtime_seconds = round(elapsed, 4), stringsAsFactors = FALSE
  )
}

for (n in COMPLEX_COUNTS) {
  cat(sprintf("[complex_count=%d] ucell ...\n", n))
  expr <- as.matrix(read.csv(file.path(DATA_DIR, "expr_n2000.csv"),
                             row.names = 1, check.names = FALSE))
  complexes <- read.csv(file.path(DATA_DIR, sprintf("complexes_n%d.csv", n)),
                        stringsAsFactors = FALSE)
  t0 <- proc.time()[["elapsed"]]
  run_ucell(expr, complexes)
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  %.3fs\n", elapsed))
  records[[length(records) + 1]] <- data.frame(
    sweep = "complex_count", value = n,
    n_cells = 2000, n_complexes = n,
    runtime_seconds = round(elapsed, 4), stringsAsFactors = FALSE
  )
}

write.csv(do.call(rbind, records),
          file.path(RESULTS_DIR, paste0("runtime_ucell_rep", Sys.getenv("RUNTIME_REPLICATE", "1"), ".csv")), row.names = FALSE)
cat("Done.\n")
