

rm(list = ls())

library(Matrix)
library(Seurat)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readxl)
library(data.table)
library(R.utils)


path <- "/proj/c.zihao/work2/1pathway/2pertu/GSE226572/"

h5_files <- list.files(
  path = path,
  pattern = "\\.h5$",
  full.names = TRUE
)

print(h5_files)


for (f in h5_files) {
  message("\nProcessing: ", f)
  outdir <- file.path(dirname(f), gsub("\\.h5$", "", basename(f)))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  counts <- Read10X_h5(f)
  class(counts)
  dim(counts)
  mat <- counts
  features <- data.frame(
    gene_id = rownames(mat),
    gene_name = rownames(mat),
    feature_type = "Gene Expression",
    stringsAsFactors = FALSE
  )
  Matrix::writeMM(mat, file.path(outdir, "matrix.mtx"))
  fwrite(
    features,
    file.path(outdir, "features.tsv"),
    sep = "\t",
    col.names = FALSE
  )
  fwrite(
    data.frame(barcode = colnames(mat)),
    file.path(outdir, "barcodes.tsv"),
    sep = "\t",
    col.names = FALSE
  )
  R.utils::gzip(file.path(outdir, "matrix.mtx"), overwrite = TRUE)
  R.utils::gzip(file.path(outdir, "features.tsv"), overwrite = TRUE)
  R.utils::gzip(file.path(outdir, "barcodes.tsv"), overwrite = TRUE)
  message("Saved 10X folder: ", outdir,
          " (", nrow(mat), " genes × ", ncol(mat), " cells)")
}
