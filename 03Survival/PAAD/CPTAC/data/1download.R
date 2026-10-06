

rm(list = ls())

library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(tibble)
library(biomaRt)
library(gdata)
set.seed(1234)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/03survival/PAAD/CPTAC/data/")



express <- fread('data_mrna_seq_v2_rsem.txt')
express <- as.data.frame(express)
express[1:5,1:5]
express$Entrez_Gene_Id <- NULL
names(express)[1] <- 'gene'
express[1:5,1:5]


exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)[1:5,1:5]

exprSet <- as.data.frame(exprSet)
exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)




pd <- fread('data_clinical_patient.txt')
head(pd)
names(pd)
pd <- subset(pd, select=c("#Patient ID" , 'Sex',  "Age" , 
                          "Vital Status"    ,                                                    
                          "Follow Up Days"    ))

names(pd) <- c("Sample", "Gender", "Age",  "OS", "OS_Time")
pd <- pd[-c(1:4),]


str(pd)
table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)

str(pd)
table(pd$OS)
pd$OS[pd$OS=="Deceased"] <- 1
pd$OS[pd$OS=="Living"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)




str(pd)
table(pd$OS_Time)
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)



pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

str(pd)
table(pd$OS)
mean(pd$OS_Time)

pd$Sample <- gsub("_", "", pd$Sample)
colnames(exprSet) <- gsub("_", "", colnames(exprSet))

pd$Sample <- gsub("-", "", pd$Sample)
colnames(exprSet) <- gsub("-", "", colnames(exprSet))




samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

