

rm(list = ls())
set.seed(42)

library(dplyr)
library(tidyr)
library(data.table)

setwd("/proj/c.zihao/work1/03survival/CRC/ACICAM/coad_silu_2022")



express <- fread("data_mrna_seq_expression.txt")

express[1:5, 1:5] 

express <- as.data.frame(express)


exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Hugo_Symbol),
                     FUN = max)
head(exprSet)[1:5, 1:5] 

exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL

min(exprSet)
max(exprSet)



clinicaldata <- fread("data_clinical_patient.txt", skip = 4)
table(clinicaldata$TUMOR_ANATOMIC_LOCATION)
table(clinicaldata$TUMOR_MORPHOLOGY)
clinicaldata[1:4, 1:4]

metadata <- data.frame(Samples=colnames(exprSet))
head(metadata)
metadata$PATIENT_ID <- gsub("-PT-01", "", metadata$Samples)

metadata <- merge(metadata, clinicaldata, by="PATIENT_ID")
head(metadata)
names(metadata)

metadata <- subset(metadata, select = c("Samples", "AGE_AT_DX", "SEX",
                                        "PATH_TUMOR_STAGE" ,            
                                       "PATH_NODES_STAGE",             
                                       "PATH_METASTASIS_STAGE",        
                                        "AJCC_PATH_STAGE",
                                       "OS_STATUS", "OS_MONTHS"))

names(metadata) <- c("Sample", "Age", "Gender", "Tstage", "Nstage", "Mstage", "Stage", "OS", "OS_Time")

str(metadata)


table(metadata$OS)
metadata$OS[metadata$OS=="0:LIVING"] <- 0
metadata$OS[metadata$OS=="1:DECEASED"] <- 1
metadata$OS <- as.numeric(as.character(metadata$OS))
table(metadata$OS)

metadata$OS_Time <- metadata$OS_Time / 12
metadata$OS_Time <- round(metadata$OS_Time, 2)

metadata <- subset(metadata, metadata$OS_Time > 0)
metadata <- metadata[!is.na(metadata$OS) & !is.na(metadata$OS_Time), ]

head(metadata)
pd <- metadata

head(pd)
pd$Sample <- chartr(old = '-', new = '_', x=pd$Sample)
colnames(exprSet)[1:5]
colnames(exprSet) <- chartr(old = '-', new = '_', x=colnames(exprSet))


colnames(exprSet) <- gsub("SER_SILU_CC_", '', colnames(exprSet))
colnames(exprSet) <- gsub("_PT_01", '', colnames(exprSet))
colnames(exprSet)[1:5]


pd$Sample <- gsub("SER_SILU_CC_", '', pd$Sample)
pd$Sample <- gsub("_PT_01", '', pd$Sample)
pd$Sample[1:5]




samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



setwd("/proj/c.zihao/work1/03survival/CRC/ACICAM/data")
exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")



