
library(UCell)

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

run_ucell <- function(expr, corum, ncores = 1L) {
  mat     <- as.matrix(expr)
  geneset <- .parse_corum(corum, rownames(mat))

  scores <- ScoreSignatures_UCell(
    matrix     = mat,
    features   = geneset,
    chunk.size = 1000L,
    ncores     = ncores,
    BPPARAM    = BiocParallel::SerialParam()
  )
  out <- t(as.matrix(scores))
  rownames(out) <- sub("_UCell$", "", rownames(out))
  out
}

score_ucell <- function(expr_path, corum_path, out_dir, ncores = 1L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr  <- as.matrix(read.csv(expr_path, row.names = 1L, check.names = FALSE))
  sep   <- if (grepl("\\.(tsv|txt)$", corum_path)) "\t" else ","
  corum <- read.table(corum_path, sep = sep, header = FALSE,
                      stringsAsFactors = FALSE)

  mat <- run_ucell(expr, corum, ncores)
  mat <- round(mat, 5)
  write.csv(mat, file.path(out_dir, "ucell.csv"), quote = FALSE)
  cat(sprintf("Saved ucell: %d x %d\n", nrow(mat), ncol(mat)))
}
