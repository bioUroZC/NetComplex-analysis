

rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)



dataset_paths <- list(
  GSE26901 = "/proj/c.zihao/work3/03Survival/STAD/GSE26901/data/exprSet.csv",
  GSE29272 = "/proj/c.zihao/work3/03Survival/STAD/GSE29272/data/exprSet.csv",
  GSE57303 = "/proj/c.zihao/work3/03Survival/STAD/GSE57303/data/exprSet.csv",
  GSE62254 = "/proj/c.zihao/work3/03Survival/STAD/GSE62254/data/exprSet.csv",
  TCGASTAD = "/proj/c.zihao/work3/03Survival/STAD/TCGASTAD/data/exprSet.csv"
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
