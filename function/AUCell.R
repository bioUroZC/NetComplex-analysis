
library(AUCell)

.parse_corum <- function(corum_df, gene_set) {
  gs <- lapply(seq_len(nrow(corum_df)), function(i) {
    genes <- trimws(unlist(strsplit(as.character(corum_df[[2L]][i]), ";")))
    genes <- genes[nzchar(genes) & genes %in% gene_set]
    genes
  })
  names(gs) <- as.character(corum_df[[1L]])
  sizes      <- lengths(gs)
  gs[sizes > 0L]
}

run_aucell <- function(expr, corum) {
  mat      <- as.matrix(expr)
  geneset  <- .parse_corum(corum, rownames(mat))
  rankings <- AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE)
  auc_obj  <- AUCell_calcAUC(geneset, rankings, verbose = FALSE)
  as.matrix(getAUC(auc_obj))
}

score_aucell <- function(expr_path, corum_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr  <- as.matrix(read.csv(expr_path, row.names = 1L, check.names = FALSE))
  sep   <- if (grepl("\\.(tsv|txt)$", corum_path)) "\t" else ","
  corum <- read.table(corum_path, sep = sep, header = FALSE,
                      stringsAsFactors = FALSE)

  mat <- run_aucell(expr, corum)
  mat <- round(mat, 5)
  write.csv(mat, file.path(out_dir, "aucell.csv"), quote = FALSE)
  cat(sprintf("Saved aucell: %d x %d\n", nrow(mat), ncol(mat)))
}
