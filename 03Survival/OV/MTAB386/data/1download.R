

rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/03survival/OV/MTAB386/data/")



exprSet <- read.csv('E.MTAB.386_eset_exprs.csv', header = T, row.names = 1)

min(exprSet)
max(exprSet)

pd <- read.csv('E.MTAB.386_eset_clindata.csv', header = T, row.names = 1)

head(pd)
names(pd)

pd <- subset(pd, select=c("unique_patient_ID",
                          "age_at_initial_pathologic_diagnosis",
                          "tumorstage" ,
                          "vital_status",
                          "days_to_death"))

names(pd) <- c('Sample', "Age", "Stage", "OS", "OS_Time")

colnames(exprSet) <- gsub("\\.", "", colnames(exprSet))
pd$Sample <- gsub("\\.", "", pd$Sample)



table(pd$OS)
pd$OS[pd$OS=="deceased"] <- 1
pd$OS[pd$OS=="living"] <- 0
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

head(pd)





samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")




