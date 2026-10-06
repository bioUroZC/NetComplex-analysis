
library(data.table)

cancer_dir <- "/proj/c.zihao/work3/03Survival/BRCA"
datasets   <- c("GSE11121", "GSE162228", "GSE17705", "GSE20685", "GSE21653",
                "GSE25055", "GSE25065", "GSE45255", "GSE61304", "TCGABRCA")
out_dir    <- file.path(cancer_dir, "02model", "data")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(out_dir, pattern = "\\.csv$", full.names = TRUE))



clin_list <- list()

for (ds in datasets) {
  pd <- fread(file.path(cancer_dir, ds, "data", "pd.csv"))

  if (all(c("OS", "OS_Time") %in% names(pd))) {
    event_col <- "OS"
    time_col  <- "OS_Time"
  } else if (all(c("EFS", "EFS_Time") %in% names(pd))) {
    event_col <- "EFS"
    time_col  <- "EFS_Time"
  } else {
    stop("no OS/OS_Time or EFS/EFS_Time columns in ", ds, "/data/pd.csv")
  }
  cat(ds, "endpoint:", event_col, "\n")

  clin_list[[ds]] <- data.frame(
    Sample   = as.character(pd$Sample),
    OS       = suppressWarnings(as.numeric(pd[[event_col]])),
    OS_Time  = suppressWarnings(as.numeric(pd[[time_col]])),
    dataset  = ds,
    endpoint = event_col
  )
}

clin <- do.call(rbind, clin_list)

keep <- !is.na(clin$Sample) & !is.na(clin$OS) & !is.na(clin$OS_Time) & clin$OS_Time > 0
clin <- clin[keep, ]

fwrite(clin, file.path(out_dir, "OS.csv"))
cat("OS.csv:", nrow(clin), "samples,", sum(clin$OS == 1), "events\n")



merge_and_save <- function(paths, id_name, out_file) {
  mats <- list()
  for (p in paths) {
    if (!file.exists(p)) stop("file not found: ", p)
    dt <- fread(p)
    m <- as.matrix(dt[, -1])
    rownames(m) <- dt[[1]]
    mats[[p]] <- m
  }

  shared <- rownames(mats[[1]])
  for (m in mats) {
    shared <- intersect(shared, rownames(m))
  }

  merged <- NULL
  for (m in mats) {
    merged <- cbind(merged, m[shared, , drop = FALSE])
  }

  bad <- !is.finite(merged)
  merged[bad] <- 0

  out <- data.frame(rownames(merged), merged, check.names = FALSE)
  names(out)[1] <- id_name
  fwrite(out, file.path(out_dir, out_file))
  cat(out_file, ":", nrow(merged), "rows x", ncol(merged), "samples,",
      sum(bad), "cells filled with 0\n")
}



paths <- file.path(cancer_dir, datasets, "data", "exprSet_filtered.csv")
merge_and_save(paths, "gene", "raw_expr.csv")

paths <- file.path(cancer_dir, datasets, "bench", "netcomplex",
                   paste0(datasets, "_netcomplex_complex_score.csv"))
merge_and_save(paths, "Complex", "netcomplex.csv")

paths <- file.path(cancer_dir, datasets, "bench", "ssgsea", "ssgsea.csv")
merge_and_save(paths, "Complex", "ssgsea.csv")

paths <- file.path(cancer_dir, datasets, "bench", "mean", "mean.csv")
merge_and_save(paths, "Complex", "mean.csv")

cat("done\n")
