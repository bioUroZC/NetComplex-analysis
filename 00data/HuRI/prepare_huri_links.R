rm(list = ls())

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

input_file <- "/proj/c.zihao/work3/00data/HuRI/HuRI.tsv"
output_file <- "/proj/c.zihao/work3/00data/HuRI/links.csv"

huri_raw <- fread(input_file, header = FALSE, sep = "\t", col.names = c("ensembl1", "ensembl2"))

huri_raw <- huri_raw %>%
  dplyr::mutate(
    ensembl1 = sub("\\..*$", "", ensembl1),
    ensembl2 = sub("\\..*$", "", ensembl2)
  )

all_ensembl <- sort(unique(c(huri_raw$ensembl1, huri_raw$ensembl2)))

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_ensembl,
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)

mapping <- mapping %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(ENSEMBL, SYMBOL) %>%
  dplyr::group_by(ENSEMBL) %>%
  dplyr::summarise(SYMBOL = dplyr::first(SYMBOL), .groups = "drop")

mapping1 <- mapping %>% dplyr::rename(ensembl1 = ENSEMBL, protein1 = SYMBOL)
mapping2 <- mapping %>% dplyr::rename(ensembl2 = ENSEMBL, protein2 = SYMBOL)

huri_links <- huri_raw %>%
  dplyr::left_join(mapping1, by = "ensembl1") %>%
  dplyr::left_join(mapping2, by = "ensembl2") %>%
  dplyr::filter(!is.na(protein1), !is.na(protein2), protein1 != "", protein2 != "") %>%
  dplyr::mutate(
    node1 = pmin(protein1, protein2),
    node2 = pmax(protein1, protein2)
  ) %>%
  dplyr::select(node1, node2) %>%
  rlang::set_names(c("protein1", "protein2")) %>%
  dplyr::filter(protein1 != protein2) %>%
  dplyr::distinct(protein1, protein2) %>%
  dplyr::mutate(score = 1) %>%
  dplyr::select(protein1, protein2, score)

write.csv(huri_links, output_file)

cat("Input pairs:", nrow(huri_raw), "\n")
cat("Mapped genes:", nrow(mapping), "\n")
cat("Output links:", nrow(huri_links), "\n")
cat("Saved to:", output_file, "\n")