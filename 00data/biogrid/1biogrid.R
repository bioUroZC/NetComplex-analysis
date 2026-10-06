


rm(list = ls())

library(data.table)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(dplyr)




setwd("/proj/c.zihao/work3/00data/biogrid")

PPIdata <- fread('BIOGRID-ORGANISM-Homo_sapiens-4.4.228.tab3.txt')

head(PPIdata)

names(PPIdata)

table(PPIdata$`Experimental System`)
table(PPIdata$`Experimental System Type`)

PPIdata <- subset(PPIdata, `Experimental System Type` == "physical")

table(PPIdata$`Organism Name Interactor A`)
table(PPIdata$`Organism Name Interactor B`)


PPIdata1 <- subset(PPIdata, select = c("Official Symbol Interactor A" ,
                                      "Official Symbol Interactor B" ,
                                       'Organism Name Interactor A',
                                       'Organism Name Interactor B'))

PPIdata1 <- subset(PPIdata1, PPIdata1$'Organism Name Interactor A' == 'Homo sapiens')
PPIdata1 <- subset(PPIdata1, PPIdata1$'Organism Name Interactor B' == 'Homo sapiens')


table(PPIdata1$'Organism Name Interactor A')
table(PPIdata1$'Organism Name Interactor B')
PPIdata1$'Organism Name Interactor A' <- NULL
PPIdata1$'Organism Name Interactor B' <- NULL

PPIdata1 <- as.data.frame(PPIdata1)

names(PPIdata1) <- c("Protein1", "Protein2")

head(PPIdata1)



protein_links <- PPIdata1 %>%
  mutate(
    node1 = pmin(Protein1, Protein2),
    node2 = pmax(Protein1, Protein2)
  ) %>%
  dplyr::select(node1, node2) 

head(protein_links)

names(protein_links) <- c("protein1", "protein2")

head(protein_links)

dim(protein_links)

protein_links <- protein_links %>%
  filter(protein1 != protein2)

dim(protein_links)

protein_links$link <- paste0(protein_links$protein1, '_', protein_links$protein2)

print(dim(protein_links))

protein_links <- protein_links %>%
  dplyr::distinct(link, .keep_all = T)


protein_links$score <- 1

protein_links$link <- NULL

str(protein_links)

print(dim(protein_links))

write.csv(protein_links, file = 'links.csv')