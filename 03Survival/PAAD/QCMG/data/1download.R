

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

setwd("/proj/c.zihao/work1/03survival/PAAD/QCMG/data/")



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

exprSet <- exprSet + 12

min(exprSet)
max(exprSet)




pd <- fread('data_clinical_patient.txt')
head(pd)

pd <- subset(pd, select=c('#Patient Identifier', 'Sex', 'Diagnosis Age', 
                          'Tumor Other Histologic Subtype', "Location"  ,                                                    
                          "Neoplasm Disease Stage American Joint Committee on Cancer Code",
                           "Status"   ,                                                     
                          "Days to Last Followup"    ))

names(pd) <- c("Sample", "Gender", "Age", "Type", "Loc", 'Stage', "OS", "OS_Time")
pd <- pd[-c(1:4),]
table(pd$Type)


pd <- subset(pd, select=c("Sample", "Gender", "Age", 'Stage', "OS", "OS_Time"  ))



str(pd)
table(pd$Gender)
pd$Gender[pd$Gender==''] <- NA
table(pd$Gender)


str(pd)
table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)


table(pd$OS)
pd$OS_binary <- NA
pd$OS_binary[pd$OS %in% c("DOD", "DOC", "Deceased - Of Unknown Cause")] <- 1
pd$OS_binary[pd$OS %in% c("AWD", "NED", "Alive - Disease Status Unknown")] <- 0

table(pd$OS_binary)
pd$OS <- pd$OS_binary
pd$OS_binary <- NULL

str(pd)



str(pd)
table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=="NA"] <- NA
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



samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")



