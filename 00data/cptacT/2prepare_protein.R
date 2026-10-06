rm(list = ls())

library(data.table)

RAW_DIR <- "/proj/c.zihao/work3/00data/cptacT"

ref <- fread(file.path("/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-BLCA.csv"), select = 1)
tcga_genes <- ref[[1]]
message("Reference genes: ", length(tcga_genes))

cancers <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD",
             "LUSC", "OV", "PDAC", "UCEC", "HNSCC")

out_dir <- "/proj/c.zihao/work3/00data/cptacT/protein"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (name in cancers) {
  path <- file.path(RAW_DIR, name, paste0(name, "_proteomics.csv"))
  if (!file.exists(path)) {
    message("[SKIP] ", name, ": file not found")
    next
  }

  dat <- fread(path)
  setnames(dat, 1, "gene")

  tumor_cols <- c("gene", colnames(dat)[-1][!grepl("\\.N$", colnames(dat)[-1])])
  dat_tumor  <- dat[, ..tumor_cols]
  dat_tumor  <- dat_tumor[gene %in% tcga_genes]

  for (col in names(dat_tumor)[-1]) {
    set(dat_tumor, which(is.na(dat_tumor[[col]])), col, 0)
    set(dat_tumor, j = col, value = round(dat_tumor[[col]], 5))
  }

  out_path <- file.path(out_dir, paste0(name, "_proteomics.csv"))
  fwrite(dat_tumor, out_path)
  message(name, ": ", ncol(dat_tumor) - 1, " tumors | ",
          nrow(dat_tumor), " genes  =>  ", out_path)
}

message("Done.")
