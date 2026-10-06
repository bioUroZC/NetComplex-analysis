

rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)



COAD <- read.csv("/proj/c.zihao/work3/00data/TCGAt/files/TCGA-COAD.csv", header = T, row.names = 1)
READ <- read.csv("/proj/c.zihao/work3/00data/TCGAt/files/TCGA-READ.csv", header = T, row.names = 1)

COAD$genes <- rownames(COAD)
READ$genes <- rownames(READ)

CRC <- merge(COAD, READ, by='genes')
rownames(CRC) <- CRC$genes
CRC$genes <- NULL

setwd('/proj/c.zihao/work3/00data/TCGAt/files/')
write.csv(CRC, file = 'TCGA-CRC.csv')