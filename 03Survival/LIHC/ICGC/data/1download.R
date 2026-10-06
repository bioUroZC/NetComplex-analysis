

rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(stringr)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/03survival/LIHC/ICGC/data/")



express <- read.table("HCCDB18_mRNA_level3.txt", header = T)

express$Entrez_ID <- NULL

express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
                     FUN = max)
head(exprSet[1:4,1:4])

exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4,1:4]

colnames(exprSet) <- chartr(old = '-', new = '_', x=colnames(exprSet))
colnames(exprSet) <- chartr(old = '.', new = '_', x=colnames(exprSet))


min(exprSet)
max(exprSet)




sam <- read.table("HCCDB18.sample.txt")
sam <- as.data.frame(t(sam))
colnames(sam) <- sam[1,]
sam <- sam[-1,]
table(sam$TYPE)
sam <- subset(sam, sam$TYPE == "HCC")
head(sam)
sam <- sam %>%
  dplyr::arrange(PATIENT, SAMPLE_NAME1)%>%
  distinct(PATIENT, .keep_all = T)


cli <-  read.table("HCCDB18.patient.txt", sep='\t')
cli <- as.data.frame(t(cli))
colnames(cli) <- cli[1,]
cli <- cli[-1,]
head(cli)

cli <- subset(cli, select=c(PATIENT_ID, AGE, GENDER,TNM_STAGE_T, STATUS, SUR))

pd <- merge(cli, sam, by="PATIENT_ID")

head(pd)
pd <- subset(pd, select=c("SAMPLE_ID", 'AGE', 'GENDER', "TNM_STAGE_T",
                          'STATUS', 'SUR'))

names(pd) <- c("Sample", "Age", "Gender", "Stage", "OS", "OS_Time")

pd$Sample <- chartr(old = '-', new = '_', x=pd$Sample)
pd$Sample <- chartr(old = '.', new = '_', x=pd$Sample)

str(pd)

pd$Age <- as.numeric(as.character(pd$Age))

str(pd)

table(pd$Stage)
pd$Stage[pd$Stage=="1"] <-"I"
pd$Stage[pd$Stage=="2"] <-"II"
pd$Stage[pd$Stage=="3"] <- "III"
pd$Stage[pd$Stage=="4"] <- "IV"
table(pd$Stage)

str(pd)


table(pd$OS)
pd$OS[pd$OS=="Alive"] <- 0
pd$OS[pd$OS=="Dead"] <- 1
pd$OS <- as.numeric(as.character(pd$OS))
str(pd)

pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
pd <- subset(pd, pd$OS_Time > 0)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]

table(pd$OS)
mean(pd$OS_Time)


str(pd)



samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")


