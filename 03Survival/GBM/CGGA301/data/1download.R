

rm(list = ls())

library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(tibble)
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gdata)
set.seed(1234)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/03survival/GBM/CGGA301/data")



exprSet <- fread("CGGA.mRNA_array_301_gene_level.20200506.txt")
exprSet <- as.data.frame(exprSet)


rownames(exprSet) <- exprSet$Gene_Name
exprSet$Gene_Name <- NULL

exprSet[1:5,1:5]

colnames(exprSet) <- gsub('_', '', colnames(exprSet))
colnames(exprSet) <- gsub('-', '_', colnames(exprSet))
colnames(exprSet) <- paste0("A", colnames(exprSet))

min(exprSet)
max(exprSet)


exprSet <- exprSet + 12

min(exprSet)
max(exprSet)



pd <- fread("CGGA.mRNA_array_301_clinical.20200506.txt")
head(pd)
names(pd)
table(pd$Histology)
table(pd$PRS_type)
table(pd$Grade)


pd <- subset(pd, pd$PRS_type == "Primary")
pd <- subset(pd, pd$Histology == "GBM")

table(pd$Histology)
table(pd$PRS_type)
table(pd$Grade)


pd <- subset(pd, select=c( "CGGA_ID","Age",  "Gender" ,
                           "Censor (alive=0; dead=1)"  ,
                           "OS"))

names(pd) <- c("Sample", "Age", "Gender", "OS", "OS_Time")
pd <- as.data.frame(pd)
head(pd)
names(pd)

str(pd)


str(pd)
table(pd$OS_Time)
pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

str(pd)
table(pd$OS)
mean(pd$OS_Time)

pd$Sample <- gsub('_', '', pd$Sample)
pd$Sample <- gsub('-', '_', pd$Sample)
pd$Sample <- paste0("A", pd$Sample)

cgga325 <- read.csv("/proj/c.zihao/work1/03survival/GBM/CGGA325/data/pd.csv", row.names = 1)
cgga693 <- read.csv("/proj/c.zihao/work1/03survival/GBM/CGGA693/data/pd.csv", row.names = 1)
shared <- c(sub("^B", "A", cgga325$Sample), sub("^C", "A", cgga693$Sample))
table(sub("_[0-9]+$", "", pd$Sample) %in% shared)
pd <- pd[!(sub("_[0-9]+$", "", pd$Sample) %in% shared), ]
nrow(pd)




samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

