

rm(list = ls())

library(dplyr)
library(data.table)
library(tidyr)




dataset_paths <- list(
  BLCA = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-BLCA.csv",
  BRCA = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-BRCA.csv",
  CRC  = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-CRC.csv",
  GBM  = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-GBM.csv",
  KIRC = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-KIRC.csv",
  LGG  = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-LGG.csv",
  LUAD = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-LUAD.csv",
  LUSC = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-LUSC.csv",
  MESO = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-MESO.csv",
  OV   = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-OV.csv",
  PAAD = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-PAAD.csv",
  PRAD = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-PRAD.csv",
  SKCM = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-SKCM.csv",
  UVM  = "/proj/c.zihao/work3/00data/TCGAt/files/TCGA-UVM.csv"
)


names(dataset_paths)




for (name in names(dataset_paths)) {
  original_path <- dataset_paths[[name]]
  expr <- read.csv(original_path, row.names = 1)
  print(expr[1:4,1:4])

  expr_filtered <- log2(expr + 1)

  print(min(expr_filtered))
  print(max(expr_filtered))
  dir_path <- dirname(original_path)
  save_path <- paste0(
    "/proj/c.zihao/work3/00data/TCGAt/exprset/", name, "_exprSet_filtered.csv"
  )
  expr_filtered <- round(expr_filtered, 5)
  write.csv(expr_filtered, file = save_path)
  print(dim(expr_filtered))
  message("Saved filtered file to: ", save_path)
}
