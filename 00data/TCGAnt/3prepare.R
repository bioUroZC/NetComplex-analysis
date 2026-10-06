

rm(list = ls())

library(dplyr)




dataset_paths <- list(
  BRCA = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-BRCA.csv",
  KIRC = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-KIRC.csv",
  LUAD = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-LUAD.csv",
  PRAD = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-PRAD.csv",

  LUSC = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-LUSC.csv",
  LIHC = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-LIHC.csv",
  HNSC = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-HNSC.csv",
  CRC  = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-CRC.csv",
  UCEC = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-UCEC.csv",
  KIRP = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-KIRP.csv",
  STAD = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-STAD.csv",
  KICH = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-KICH.csv",
  BLCA = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-BLCA.csv",
  ESCA = "/proj/c.zihao/work3/00data/TCGAnt/files/TCGA-ESCA.csv"
)



names(dataset_paths)




for (name in names(dataset_paths)) {
  original_path <- dataset_paths[[name]]
  expr <- read.csv(original_path, row.names = 1)
  print(expr[1:4,1:4])

  expr_filtered <- log2(expr + 1)
  print(min(expr_filtered))
  print(max(expr_filtered))
  save_path <- paste0(
    "/proj/c.zihao/work3/00data/TCGAnt/exprset/", name, "_exprSet_filtered.csv"
  )
  expr_filtered <- round(expr_filtered, 5)
  write.csv(expr_filtered, file = save_path)
  print(dim(expr_filtered))
  message("Saved filtered file to: ", save_path)
}
