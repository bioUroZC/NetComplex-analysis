

rm(list = ls())

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)



setwd("/proj/c.zihao/work3/00data/string")

protein_links <- fread("9606.protein.links.v12.0.txt", header = TRUE)
protein_links <- as.data.frame(protein_links)
str(protein_links)

min(protein_links$combined_score)

a1 <- dim(protein_links)
print(a1)

protein_links <- protein_links %>%
  dplyr::filter(combined_score > 100)
print(head(protein_links))

a2 <- dim(protein_links)
print(a2)

print(a2[1]/a1[1])

geneinfor <- fread("9606.protein.aliases.v12.0.txt")
geneinfor  <- as.data.frame(geneinfor)
table(geneinfor$source)
hugo_rows <- geneinfor[grep("HUGO", geneinfor$source), ]
print(head(hugo_rows))
names(hugo_rows)[1] <- "protein_id"



id_to_alias <- setNames(hugo_rows$alias, hugo_rows$protein_id)

protein_links$protein1_alias <- id_to_alias[protein_links$protein1]
protein_links$protein2_alias <- id_to_alias[protein_links$protein2]


print(head(protein_links))

protein_links$protein1 <- NULL
protein_links$protein2 <- NULL


head(protein_links)

protein_links <- protein_links %>%
  mutate(
    node1 = pmin(protein1_alias, protein2_alias),
    node2 = pmax(protein1_alias, protein2_alias)
  ) %>%
  dplyr::select(combined_score, node1, node2) 

head(protein_links)

names(protein_links) <- c("score", "protein1", "protein2")

head(protein_links)

dim(protein_links)

protein_links <- protein_links %>%
  filter(protein1 != protein2)

dim(protein_links)

protein_links$link <- paste0(protein_links$protein1, '_', protein_links$protein2)

print(dim(protein_links))

protein_links <- protein_links %>%
  dplyr::distinct(link, .keep_all = T)

protein_links$link <- NULL

str(protein_links)

print(dim(protein_links))

protein_links$score <- protein_links$score / 1000

protein_links <- protein_links %>%
     dplyr::select(protein1, protein2, score)

write.csv(protein_links, file = 'links100.csv')

print(head(protein_links))
