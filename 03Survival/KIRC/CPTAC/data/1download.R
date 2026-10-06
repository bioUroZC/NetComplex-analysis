

rm(list = ls())

library(dplyr)
library(tibble)
library(dplyr)
library(data.table)
library(Biobase)
library(tibble)
library(stringr)
library(gdata)
library(biomaRt)
library(openxlsx)
library(tidyr)




setwd("/proj/c.zihao/work3/03Survival/KIRC/CPTAC/data/")

exprsess <- fread("CPTAC-3.star_tpm.tsv")

exprsess[1:5,1:5]

exprsess <- as.data.frame(exprsess)

exprsess[1:5,1:5]

names(exprsess)[1] <- 'gene_id'

exprsess$gene_id <- sub("\\..*", "", exprsess$gene_id)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- exprsess$gene_id
genes_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
                    filters = 'ensembl_gene_id', 
                    values = gene_ids, 
                    mart = mart)
head(genes_info)

genes_info <- subset(genes_info, select=c('ensembl_gene_id', 'external_gene_name'))
names(genes_info) <- c('gene_id', 'gene_name')
data <- merge(genes_info, exprsess, by="gene_id")
data$gene_id <- NULL

data[1:5,1:5]




clinical <- fread("CPTAC-3.clinical.tsv")
clinical <- as.data.frame(clinical)
names(clinical)
table(clinical$primary_site)
table(clinical$sample_type.samples)
table(clinical$tissue_type.samples)


clinical <- subset(clinical, clinical$primary_site == "Kidney")
clinical <- subset(clinical, clinical$sample_type.samples == "Primary Tumor")
clinical <- subset(clinical, clinical$primary_diagnosis.diagnoses == "Renal cell carcinoma, NOS")

table(clinical$sample_type.samples)
table(clinical$tissue_type.samples)

table(clinical$sample_type_id.samples)
table(clinical$tumor_descriptor.samples)
table(clinical$primary_diagnosis.diagnoses)

table(clinical$sample_ordinal.samples)
write.csv(clinical, file = 'clinical.csv')

names(clinical)
head(clinical)


pd <- clinical %>%
 dplyr::select(
    Sample = sample,
    Age = age_at_earliest_diagnosis_in_years.diagnoses.xena_derived,
    Gender = gender.demographic,
    Tstage = ajcc_pathologic_t.diagnoses,
    Nstage = ajcc_pathologic_n.diagnoses,
    Mstage = ajcc_pathologic_m.diagnoses,
    Stage = ajcc_pathologic_stage.diagnoses,
    Grade = tumor_grade.diagnoses,
    OS = vital_status.demographic,
    OS_Time_raw_death = days_to_death.demographic,
    OS_Time_raw_followup = days_to_last_follow_up.diagnoses,
    order = sample_ordinal.samples
  )


table(pd$Age)
str(pd)
pd$Age <- round(pd$Age, 2)



table(pd$Gender)
pd$Gender[pd$Gender=="female"] <- "Female"
pd$Gender[pd$Gender=="male"] <- "Male"
table(pd$Gender)

table(pd$Tstage)
pd$Tstage[pd$Tstage=="TX"] <- NA
pd$Tstage <- sub("^(T[1-4]).*", "\\1", pd$Tstage)
pd$Tstage[pd$Tstage==""] <- NA
table(pd$Tstage)

table(pd$Nstage)
pd$Nstage[pd$Nstage=="NX"] <- NA
pd$Nstage[pd$Nstage==""] <- NA
table(pd$Nstage)

table(pd$Mstage)
pd$Mstage[pd$Mstage=="MX"] <- NA
pd$Mstage[pd$Mstage==""] <- NA
table(pd$Mstage)

head(pd)

table(pd$Stage)
pd$Stage[pd$Stage=="Not Reported"] <- NA
pd$Stage[pd$Stage=="Stage I"] <- "I"
pd$Stage[pd$Stage=="Stage II"] <- "II"
pd$Stage[pd$Stage=="Stage III"] <- "III"
pd$Stage[pd$Stage=="Stage IV"] <- "IV"
table(pd$Stage)

table(pd$Grade)


table(pd$OS)
pd$OS[pd$OS=="Alive"] <- 0
pd$OS[pd$OS=="Dead"] <- 1
pd$OS[pd$OS=="Not Reported"] <- NA
table(pd$OS)

head(pd)


pd$OS_Time <- ifelse(
  pd$OS == 1,
  pd$OS_Time_raw_death,
  pd$OS_Time_raw_followup
)


pd$OS_Time_raw_death <- NULL
pd$OS_Time_raw_followup <- NULL

head(pd)

data[1:5,1:5]
pd$Sample[1:5]


str(pd)
pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)


str(pd)

table(pd$OS)

pd$ID <- substr(x=pd$Sample, start = 1, stop = 9)

pd <- pd %>%
  dplyr::arrange(ID, order)

head(pd)

pd_unique <- pd %>%
  group_by(ID) %>%
  slice_min(order, with_ties = FALSE) %>%  
  ungroup()

pd <- pd_unique

table_IDs <- table(pd$ID)
table_IDs[table_IDs > 1]  


table(pd$Age)
table(pd$Gender)
table(pd$Grade)
table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)


str(pd)



keep_samples <- pd$Sample
data_filtered <- data[, c("gene_name", intersect(colnames(data), keep_samples))]
data_filtered[1:5,1:5]

exprSet <- aggregate(x = data_filtered[, 2:ncol(data_filtered)],
                     by = list(data_filtered$gene_name),
                     FUN = max)

exprSet[1:5,1:5]
exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[-1, ]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]





samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

pd$Sample <- gsub('-', '', pd$Sample)
colnames(exprSet) <- gsub('-', '', colnames(exprSet))

pd$Sample
colnames(exprSet)



exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

