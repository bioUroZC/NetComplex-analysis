

rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)


dataset_paths <- list(
  ACICAM = "/proj/c.zihao/work3/03Survival/CRC/ACICAM/data/exprSet.csv",
  GSE12945 = "/proj/c.zihao/work3/03Survival/CRC/GSE12945/data/exprSet.csv",
  GSE17536 = "/proj/c.zihao/work3/03Survival/CRC/GSE17536/data/exprSet.csv",
  GSE17537 = "/proj/c.zihao/work3/03Survival/CRC/GSE17537/data/exprSet.csv",
  GSE28722 = "/proj/c.zihao/work3/03Survival/CRC/GSE28722/data/exprSet.csv",
  GSE39582 = "/proj/c.zihao/work3/03Survival/CRC/GSE39582/data/exprSet.csv",
  TCGACRC = "/proj/c.zihao/work3/03Survival/CRC/TCGACRC/data/exprSet.csv"
)

links = read.csv("/proj/c.zihao/work3/00data/string/links.csv", row.names = 1)
common_genes <- unique(union(links$protein1, links$protein2))
message("Number of common genes: ", length(common_genes))

names(dataset_paths)




for (name in names(dataset_paths)) {
  original_path <- dataset_paths[[name]]
  expr <- read.csv(original_path, row.names = 1, check.names = FALSE)
  expr_filtered <- expr[which(rownames(expr) %in% common_genes), ]
  expr_filtered <- log2(expr_filtered + 1) 
  print(min(expr_filtered))
  print(max(expr_filtered))
  dir_path <- dirname(original_path)
  base_name <- tools::file_path_sans_ext(basename(original_path))
  save_path <- file.path(dir_path, paste0(base_name, "_filtered.csv"))
  expr_filtered <- round(expr_filtered, 5)
  write.csv(expr_filtered, file = save_path)
  print(dim(expr_filtered))
  message("Saved filtered file to: ", save_path)
}
