
EXPR_DIR <- "/proj/c.zihao/work3/00data/cptacT/exprset"
CORUM    <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
OUT_DIR  <- "/proj/c.zihao/work3/01protein/01random/01_random_sets"
CANCERS  <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")

N_REPS    <- 50                          
SEED_BASE <- 1                                   

if (dir.exists(OUT_DIR)) unlink(OUT_DIR, recursive = TRUE)
dir.create(OUT_DIR, recursive = TRUE)

all_genes <- sort(rownames(read.csv(file.path(EXPR_DIR, "KIRC_exprSet_filtered.csv"),
                                    row.names = 1, check.names = FALSE)))
for (cancer in CANCERS) {
  f <- file.path(EXPR_DIR, paste0(cancer, "_exprSet_filtered.csv"))
  genes <- sort(rownames(read.csv(f, row.names = 1, check.names = FALSE)))
  if (!identical(genes, all_genes)) stop(cancer, ": gene list differs from KIRC")
}

corum <- read.csv(CORUM, stringsAsFactors = FALSE)
members <- list()
for (i in 1:nrow(corum)) {
  g <- trimws(strsplit(corum$Genes[i], ";")[[1]])
  members[[i]] <- unique(g[g %in% all_genes])
}
corum_genes <- unique(unlist(members))

pool <- setdiff(all_genes, corum_genes)
cat("all genes:", length(all_genes), "| CORUM genes:", length(corum_genes),
    "| genes to draw from:", length(pool), "\n")

for (rep_i in 1:N_REPS) {
  set.seed(SEED_BASE + rep_i)

  random_genes <- character(nrow(corum))
  for (i in 1:nrow(corum)) {
    n <- length(members[[i]])
    random_genes[i] <- paste(sample(pool, n), collapse = ";")                      
  }

  out <- data.frame(Complex = corum$Complex, Genes = random_genes)
  write.csv(out, file.path(OUT_DIR, sprintf("rep%03d.csv", rep_i)), row.names = FALSE)
  cat(sprintf("rep%03d done\n", rep_i))
}

cat("[done]", N_REPS, "random sets ->", OUT_DIR, "\n")
