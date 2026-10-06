library(data.table)

DATA_ROOT <- "/proj/c.zihao/work3/06singleCell"
CORUM_FILE <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
RESULT_DIR <- "/proj/c.zihao/work3/10Cellstudy/results"

DATASETS <- c("GSE96583", "GSE226572")
METHODS <- c("netComplex", "ssgsea", "aucell", "ucell", "gsva", "netmean")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

corum <- fread(CORUM_FILE)
gene_of <- setNames(corum$Genes, corum$Complex)

read_auc <- function(dataset, file_name, set_name) {
  auc <- fread(file.path(DATA_ROOT, dataset, "04analysis", file_name))
  auc[, `:=`(dataset = dataset, set = set_name)]
  auc
}

long <- rbindlist(list(
  rbindlist(lapply(DATASETS, read_auc, file_name = "auc.csv", set_name = "positive")),
  rbindlist(lapply(DATASETS, read_auc, file_name = "auc_neg.csv", set_name = "negative"))
))
long <- long[method %in% METHODS]
long[, genes := gene_of[complex]]
long[, n_members := lengths(strsplit(genes, ";"))]
setcolorder(long, c("dataset", "set", "complex", "n_members", "method", "auc", "pr_auc", "genes"))
setorder(long, dataset, set, complex, -auc)
fwrite(long, file.path(RESULT_DIR, "percomplex_auc_long.csv"))

wide <- dcast(long, dataset + set + complex + n_members + genes ~ method,
              value.var = "auc")
setorder(wide, dataset, set, -netComplex)
fwrite(wide, file.path(RESULT_DIR, "percomplex_auc_wide.csv"))

cat("datasets:", paste(DATASETS, collapse = ", "), "\n")
cat("methods :", paste(METHODS, collapse = ", "), "\n\n")
print(long[, .(complexes = uniqueN(complex),
               netComplex = round(mean(auc[method == "netComplex"]), 3),
               ssgsea     = round(mean(auc[method == "ssgsea"]), 3),
               aucell     = round(mean(auc[method == "aucell"]), 3)),
           by = .(dataset, set)], row.names = FALSE)
cat("\nWrote ", file.path(RESULT_DIR, "percomplex_auc_wide.csv"), "\n", sep = "")
