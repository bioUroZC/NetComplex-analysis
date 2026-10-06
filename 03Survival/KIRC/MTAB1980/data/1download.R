

rm(list = ls())

library(data.table)
library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(tibble)
library(stringr)
library(gdata)
library(biomaRt)
library(openxlsx)

Sys.setenv("VROOM_CONNECTION_SIZE"=131072*16)




setwd("/proj/c.zihao/work3/03Survival/KIRC/MTAB1980/data/E-MTAB-1980")

exprset <- fread("ccRCC_exp_log_quantile_normalized.txt")
exprset[1:5,1:5]
exprset <- as.data.frame(exprset)
exprset[1:5,1:5]

exprset$`Feature Numbers` <- NULL



setwd("/proj/c.zihao/work3/03Survival/KIRC/MTAB1980/data/")

anno=read.table("refGene.txt", sep="\t")

anno=anno[, c(2, 13)]

names(anno)=c("SystematicName", "symbol")
head(anno)
dim(anno)
anno <- anno %>%
  dplyr::distinct(SystematicName, .keep_all = T)
dim(anno)



data <- merge(anno, exprset, by="SystematicName")
data$SystematicName <- NULL
data$REF <- NULL 

data[1:5,1:5]

exprSet <- aggregate(x = data[, 2:ncol(data)],
                     by = list(data$symbol),
                     FUN = max)

exprSet[1:5,1:5]
exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

colnames(exprSet) <- chartr(old = '-', new = '_', x=colnames(exprSet))





pd <- read.xlsx("41588_2013_BFng2699_MOESM35_ESM.xlsx")
names(pd)
head(pd)

pd <- subset(pd, select = c("sample.ID", "Age", "Sex","Stage.at.diagnosis",
                            "Fuhrman.grade", "outcome",
                            "observation.period.(month)"))

names(pd) <- c("Sample","Age", "Gender", "Stage", "Grade", "OS", "OS_Time")

str(pd)

pd$Age <- as.numeric(as.character(pd$Age))
pd$Age <- round(pd$Age, 2)


table(pd$Age)



table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)

pd$Tstage <- as.numeric(sub(".*T(\\d+).*", "\\1", pd$Stage))
pd$Nstage <- as.numeric(sub(".*N(\\d+).*", "\\1", pd$Stage))
pd$Mstage <- as.numeric(sub(".*M(\\d+).*", "\\1", pd$Stage))

pd$Tstage <- paste0("T", pd$Tstage)
pd$Nstage <- paste0("N", pd$Nstage)
pd$Mstage <- paste0("M", pd$Mstage)


table(pd$Grade)
pd$Grade <- paste0("G", pd$Grade)
pd$Grade[pd$Grade=="GNA"] <- NA
pd$Grade[pd$Grade=="Gundetermined"] <- NA
table(pd$Grade)


pd$Stage <- NULL

table(pd$OS)
pd$OS[pd$OS=="alive"] <- 0
pd$OS[pd$OS=="dead"] <- 1
table(pd$OS)
pd$OS <- as.numeric(as.character(pd$OS))
str(pd)
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)
str(pd)


pd <- pd %>%
  dplyr::select("Sample", "Age", "Gender", "Tstage", "Nstage", "Mstage",
                "Grade",  "OS", "OS_Time")


table(pd$Age)
table(pd$Gender)
table(pd$Grade)
table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)


pd$Sample <- chartr(old = '-', new = '_', x=pd$Sample)



samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)



exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

