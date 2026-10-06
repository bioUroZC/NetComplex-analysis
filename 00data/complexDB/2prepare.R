rm(list = ls())

library(data.table)

setwd("/proj/c.zihao/work3/00data/complexDB")

dat <- fread(
  "/proj/c.zihao/work3/00data/complexDB/complexportal_parsed.tsv",
  sep = "\t", quote = "", fill = Inf, header = TRUE,
  data.table = FALSE
)

cat("Raw parsed complexes:", nrow(dat), "\n")

out <- dat[, c("complex_name", "genes", "stoichiometry")]

out <- out[!is.na(out$complex_name) & out$complex_name != "" &
           !is.na(out$genes) & out$genes != "", ]

cat("After removing empty names/genes:", nrow(out), "\n")

tcga_genes <- rownames(read.csv(
  "/proj/c.zihao/work3/00data/TCGAnt/exprset/BLCA_exprSet_filtered.csv",
  row.names = 1,
  nrows = 0
))

has_overlap <- sapply(out$genes, function(g) {
  genes <- trimws(strsplit(g, ";")[[1]])
  any(genes %in% tcga_genes)
})

out <- out[has_overlap, ]

cat("After gene filter:", nrow(out), "complexes retained.\n")

out$complex_name <- gsub("[^A-Za-z0-9]+", "_", out$complex_name)
out$complex_name <- gsub("^_+|_+$", "", out$complex_name)

names(out) <- c("Complex", "Genes", "Stoichiometry")

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

is_trivial <- function(s) {
  if (is.na(s) || s == "") return(TRUE)
  vals <- trimws(strsplit(s, ";")[[1]])
  vals <- suppressWarnings(as.numeric(vals))
  if (any(is.na(vals))) return(TRUE)                                    
  if (any(vals == 0)) return(TRUE)                                            
  all(vals == 1)
}
trivial <- sapply(out$Stoichiometry, is_trivial)
out$Stoichiometry[trivial] <- NA

has_stoi <- !is.na(out$Stoichiometry)
cat("Complexes with non-trivial stoichiometry:", sum(has_stoi), "/",
    nrow(out), "\n")

print(head(out))

write.csv(out[, c("Complex", "Genes")], file = "complexData.csv", row.names = FALSE)
cat("[OK] Saved:", file.path(getwd(), "complexData.csv"), "\n")

write.csv(out, file = "complexData_stoich.csv", row.names = FALSE)
cat("[OK] Saved:", file.path(getwd(), "complexData_stoich.csv"), "\n")
