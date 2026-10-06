rm(list = ls())

library(data.table)

setwd("/proj/c.zihao/work3/00data/corum")

dat <- fread(
  "/proj/c.zihao/work3/00data/corum/corum_humanComplexes.txt",
  sep = "\t", quote = "", fill = Inf, header = TRUE,
  data.table = FALSE
)

out <- dat[, c("complex_name", "subunits_gene_name", "subunits_stoechiometrie")]

out <- out[!is.na(out$complex_name) & out$complex_name != "" &
           !is.na(out$subunits_gene_name) & out$subunits_gene_name != "", ]

tcga_genes <- rownames(read.csv(
  "/proj/c.zihao/work3/00data/TCGAnt/exprset/BLCA_exprSet_filtered.csv",
  row.names = 1,
  nrows = 0
))

has_overlap <- sapply(out$subunits_gene_name, function(g) {
  genes <- trimws(strsplit(g, ";")[[1]])
  any(genes %in% tcga_genes)
})

out <- out[has_overlap, ]

cat("After gene filter:", nrow(out), "complexes retained.\n")

out$complex_name <- gsub("[^A-Za-z0-9]+", "_", out$complex_name)

out$complex_name <- gsub("^_+|_+$", "", out$complex_name)

names(out) <- c("Complex", "Genes", "Stoichiometry")

is_incomplete <- function(s) {
  if (is.na(s) || !grepl("[0-9]", s)) return(TRUE)
  vals <- trimws(strsplit(s, ";")[[1]])
  any(vals == "") || any(suppressWarnings(is.na(as.numeric(vals))))
}
out$Stoichiometry[sapply(out$Stoichiometry, is_incomplete)] <- NA

gene_key <- vapply(strsplit(out$Genes, ";"), function(g) {
  paste(sort(unique(trimws(g))), collapse = ";")
}, character(1))
out <- out[!duplicated(gene_key), ]

cat("After removing duplicated complexes:", nrow(out), "complexes retained.\n")

out$Complex <- paste0("Comp_", out$Complex)

dup_occurrence <- ave(seq_along(out$Complex), out$Complex, FUN = seq_along)
n_disambiguated <- sum(dup_occurrence > 1)
out$Complex[dup_occurrence > 1] <- paste0(out$Complex[dup_occurrence > 1],
                                          "_v", dup_occurrence[dup_occurrence > 1])
cat("Disambiguated", n_disambiguated, "complex name(s) with distinct gene sets.\n")

has_stoi <- !is.na(out$Stoichiometry)
cat("Complexes with complete stoichiometry:", sum(has_stoi), "/",
    nrow(out), "\n")

print(head(out))

write.csv(out[, c("Complex", "Genes")], file = "complexData.csv", row.names = FALSE)
write.csv(out, file = "complexData_stoich.csv", row.names = FALSE)
