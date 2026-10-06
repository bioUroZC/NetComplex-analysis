

rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)



dataset_paths <- list(
  CPTAC = "/proj/c.zihao/work3/03Survival/KIRC/CPTAC/data/exprSet.csv",
  GSE167573 = "/proj/c.zihao/work3/03Survival/KIRC/GSE167573/data/exprSet.csv",
  GSE29609 = "/proj/c.zihao/work3/03Survival/KIRC/GSE29609/data/exprSet.csv",
  MTAB1980 = "/proj/c.zihao/work3/03Survival/KIRC/MTAB1980/data/exprSet.csv",
  TCGAKIRC = "/proj/c.zihao/work3/03Survival/KIRC/TCGAKIRC/data/exprSet.csv"
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
