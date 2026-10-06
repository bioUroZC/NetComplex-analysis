rm(list = ls())

library(data.table)

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")

RNA_DIR <- "/proj/c.zihao/work3/00data/cptacT/exprset"
PROT_DIR <- "/proj/c.zihao/work3/00data/cptacT/protein"
OUT_DIR <- "/proj/c.zihao/work3/07multi/00data"
SUMMARY_PATH <- "/proj/c.zihao/work3/07multi/00data/sample_overlap_summary.csv"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()

for (cancer in CANCERS) {
  rna_path <- file.path(RNA_DIR, paste0(cancer, "_exprSet_filtered.csv"))
  prot_path <- file.path(PROT_DIR, paste0(cancer, "_proteomics.csv"))

  rna <- fread(rna_path)
  prot <- fread(prot_path)

  gene_col <- colnames(rna)[1]
  common_genes <- rna[[gene_col]][rna[[gene_col]] %in% prot[[gene_col]]]
  rna <- rna[match(common_genes, get(gene_col))]
  prot <- prot[match(common_genes, get(gene_col))]

  rna_samples <- colnames(rna)[-1]
  prot_samples <- colnames(prot)[-1]
  common_samples <- rna_samples[rna_samples %in% prot_samples]

  keep_cols <- c(colnames(rna)[1], common_samples)
  rna_matched <- rna[, ..keep_cols]
  prot_matched <- prot[, ..keep_cols]

  rna_out <- file.path(OUT_DIR, paste0(cancer, "_RNA.csv"))
  prot_out <- file.path(OUT_DIR, paste0(cancer, "_protein.csv"))

  fwrite(rna_matched, rna_out)
  fwrite(prot_matched, prot_out)

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    cancer = cancer,
    matched_genes = length(common_genes),
    rna_samples = length(rna_samples),
    protein_samples = length(prot_samples),
    matched_samples = length(common_samples),
    status = "done"
  )

  message("[", cancer, "] matched=", length(common_samples),
          " -> ", rna_out, ", ", prot_out)
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, SUMMARY_PATH)

message("[done] matched matrices saved to ", OUT_DIR)
