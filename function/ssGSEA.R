
library(GSVA)

.parse_corum <- function(corum_df, gene_set) {
  gs <- lapply(seq_len(nrow(corum_df)), function(i) {
    genes <- trimws(unlist(strsplit(as.character(corum_df[[2L]][i]), ";")))
    genes[nzchar(genes) & genes %in% gene_set]
  })
  names(gs) <- as.character(corum_df[[1L]])
  gs                                        
}

run_ssgsea <- function(expr, corum) {
  mat      <- as.matrix(expr)
  all_gs   <- .parse_corum(corum, rownames(mat))
  nonempty <- all_gs[lengths(all_gs) >= 1L]
  empty    <- names(all_gs)[lengths(all_gs) == 0L]

  result <- if (length(nonempty) > 0L) {
    param <- ssgseaParam(exprData = mat, geneSets = nonempty, normalize = FALSE)
    gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
  } else {
    matrix(nrow = 0L, ncol = ncol(mat), dimnames = list(NULL, colnames(mat)))
  }

  scored        <- if (!is.null(result)) rownames(result) else character(0)
  gs_dropped    <- setdiff(names(nonempty), scored)
  if (length(gs_dropped) > 0L) {
    zero_mat <- matrix(0, nrow = length(gs_dropped), ncol = ncol(mat),
                       dimnames = list(gs_dropped, colnames(mat)))
    result <- rbind(result, zero_mat)
  }

  if (length(empty) > 0L) {
    zero_mat <- matrix(0, nrow = length(empty), ncol = ncol(mat),
                       dimnames = list(empty, colnames(mat)))
    result <- rbind(result, zero_mat)
  }

  result[names(all_gs), , drop = FALSE]
}

score_ssgsea <- function(expr_path, corum_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr  <- as.matrix(read.csv(expr_path, row.names = 1L, check.names = FALSE))
  sep   <- if (grepl("\\.(tsv|txt)$", corum_path)) "\t" else ","
  corum <- read.table(corum_path, sep = sep, header = TRUE,
                      stringsAsFactors = FALSE)

  mat <- run_ssgsea(expr, corum)
  mat <- round(mat, 5)
  write.csv(mat, file.path(out_dir, "ssgsea.csv"), quote = FALSE)
  cat(sprintf("Saved ssgsea: %d x %d\n", nrow(mat), ncol(mat)))
}
