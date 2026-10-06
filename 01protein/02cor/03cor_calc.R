
RNA_ROOT  <- "/proj/c.zihao/work3/01protein/01stringTest"
PROT_ROOT <- "/proj/c.zihao/work3/01protein/02cor/02scored"
OUT_DIR   <- "/proj/c.zihao/work3/01protein/02cor"

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHODS <- c("netcomplex", "mean", "min", "zscore", "netmean",
             "perpas_adapted", "gsva", "ssgsea", "plage")
MIN_SAMPLES <- 10

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load_rna <- function(method, cancer) {
  file <- if (method == "netcomplex") {
    file.path(RNA_ROOT, "netcomplex", "complex_score", cancer,
             sprintf("%s_netcomplex_complex_score.csv", cancer))
  } else {
    file.path(RNA_ROOT, "bench", method, cancer, paste0(method, ".csv"))
  }
  if (!file.exists(file)) return(NULL)
  as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
}

load_protein <- function(method, cancer) {
  file <- if (method == "netcomplex") {
    file.path(PROT_ROOT, "netcomplex", "complex_score", cancer,
             sprintf("%s_netcomplex_complex_score.csv", cancer))
  } else {
    file.path(PROT_ROOT, method, cancer, paste0(method, ".csv"))
  }
  if (!file.exists(file)) return(NULL)
  as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
}

spearman_per_row <- function(rna_mat, prot_mat) {
  n <- nrow(rna_mat)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    x <- rna_mat[i, ]
    y <- prot_mat[i, ]
    if (!all(is.finite(x)) || !all(is.finite(y))) next
    if (sd(x) == 0 || sd(y) == 0) next
    out[i] <- suppressWarnings(cor(x, y, method = "spearman"))
  }
  out
}

all_res <- list()

for (cancer in CANCERS) {
  cat("[", cancer, "] ...\n")
  rna_scores  <- lapply(METHODS, load_rna, cancer = cancer)
  prot_scores <- lapply(METHODS, load_protein, cancer = cancer)
  names(rna_scores) <- names(prot_scores) <- METHODS

  missing_rna  <- METHODS[vapply(rna_scores, is.null, logical(1))]
  missing_prot <- METHODS[vapply(prot_scores, is.null, logical(1))]
  if (length(missing_rna) > 0 || length(missing_prot) > 0) {
    stop(sprintf("[%s] missing score matrices -- rna: %s, protein: %s", cancer,
                 paste(missing_rna, collapse = ", "), paste(missing_prot, collapse = ", ")))
  }

  for (rna_method in METHODS) {
    for (prot_method in METHODS) {
      rna_mat  <- rna_scores[[rna_method]]
      prot_mat <- prot_scores[[prot_method]]

      common_complexes <- intersect(rownames(rna_mat), rownames(prot_mat))
      common_samples   <- intersect(colnames(rna_mat), colnames(prot_mat))
      if (length(common_samples) < MIN_SAMPLES || length(common_complexes) == 0) {
        cat(sprintf("    skip %s vs %s: samples=%d complexes=%d\n",
                    rna_method, prot_method, length(common_samples), length(common_complexes)))
        next
      }

      rna_sub  <- rna_mat[common_complexes, common_samples, drop = FALSE]
      prot_sub <- prot_mat[common_complexes, common_samples, drop = FALSE]
      r <- spearman_per_row(rna_sub, prot_sub)

      keep <- !is.na(r)
      if (sum(keep) == 0) next

      all_res[[length(all_res) + 1]] <- data.frame(
        cancer = cancer,
        rna_method = rna_method,
        protein_method = prot_method,
        complex = common_complexes[keep],
        n_samples = length(common_samples),
        spearman_r = round(r[keep], 5),
        stringsAsFactors = FALSE
      )
    }
    cat(sprintf("    rna=%-11s done (%d protein references)\n", rna_method, length(METHODS)))
  }
}

res <- do.call(rbind, all_res)
write.csv(res, file.path(OUT_DIR, "cor_complex_level.csv"), row.names = FALSE)

cat("\n=== mean per-complex Spearman r, by (rna_method, protein_method), pooled across cohorts ===\n")
cell_means <- aggregate(spearman_r ~ rna_method + protein_method, data = res, FUN = mean)
print(cell_means[order(cell_means$protein_method, -cell_means$spearman_r), ], row.names = FALSE)

cat("\nOutput ->", file.path(OUT_DIR, "cor_complex_level.csv"), "\n")
