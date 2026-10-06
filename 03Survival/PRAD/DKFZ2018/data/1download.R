

rm(list=ls())

library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)
library(biomaRt)
library(readxl)
options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work3/03Survival/PRAD/DKFZ2018/data")



express <- fread('data_mrna_seq_rpkm.txt')
express <- as.data.frame(express)
express$Entrez_Gene_Id <- NULL

colnames(express)
names(express)[1] <- 'gene'
str(express)

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)

exprSet <- as.data.frame(exprSet)
exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

colnames(exprSet) <- sub("(^[^_]+_[^_]+)_(.*)$", "\\1", colnames(exprSet))





pd <- fread("data_clinical_patient.txt")
pd <- as.data.frame(pd)
pd <- pd[-(1:5), ]
head(pd)
names(pd)

pd <- subset(pd, select=c("#Patient Identifier", 
                          "Diagnosis Age",
                          "Stage"  ,
                          "BCR Status", "Time from Surgery to BCR/Last Follow Up" ))

names(pd) <- c("Sample", "Age", "Stage", "OS", "OS_Time")
head(pd)
str(pd)

table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)
str(pd)



table(pd$OS)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)


table(pd$OS_Time)
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
table(pd$OS_Time)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)

table(pd$OS)
mean(pd$OS_Time)




samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

