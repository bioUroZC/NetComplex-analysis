
library(GSVA)

.parse_corum <- function(corum_df, gene_set) {
  gs <- lapply(seq_len(nrow(corum_df)), function(i) {
    genes <- trimws(unlist(strsplit(as.character(corum_df[[2L]][i]), ";")))
    genes[nzchar(genes) & genes %in% gene_set]
  })
  names(gs) <- as.character(corum_df[[1L]])
  gs                                        
}

run_plage <- function(expr, corum) {
  mat     <- as.matrix(expr)
  all_gs  <- .parse_corum(corum, rownames(mat))
  geneset <- all_gs[lengths(all_gs) >= 1L]
  empty   <- names(all_gs)[lengths(all_gs) == 0L]

  single  <- geneset[lengths(geneset) == 1L]
  multi   <- geneset[lengths(geneset) >= 2L]

  result_rows <- list()

  if (length(multi) > 0L) {
    param <- plageParam(exprData = mat, geneSets = multi)
    result_rows[["multi"]] <- gsva(param, verbose = FALSE,
                                   BPPARAM = BiocParallel::SerialParam())
  }

  if (length(single) > 0L) {
    single_mat <- do.call(rbind, lapply(names(single), function(cx) {
      mat[single[[cx]], , drop = FALSE]
    }))
    rownames(single_mat) <- names(single)
    result_rows[["single"]] <- single_mat
  }

  result <- do.call(rbind, result_rows)

  scored        <- if (!is.null(result)) rownames(result) else character(0)
  multi_dropped <- setdiff(names(multi), scored)
  if (length(multi_dropped) > 0L) {
    zero_mat <- matrix(0, nrow = length(multi_dropped), ncol = ncol(mat),
                       dimnames = list(multi_dropped, colnames(mat)))
    result <- rbind(result, zero_mat)
  }

  if (length(empty) > 0L) {
    zero_mat <- matrix(0, nrow = length(empty), ncol = ncol(mat),
                       dimnames = list(empty, colnames(mat)))
    result <- rbind(result, zero_mat)
  }

  result[names(all_gs), , drop = FALSE]
}

score_plage <- function(expr_path, corum_path, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr  <- as.matrix(read.csv(expr_path, row.names = 1L, check.names = FALSE))
  sep   <- if (grepl("\\.(tsv|txt)$", corum_path)) "\t" else ","
  corum <- read.table(corum_path, sep = sep, header = TRUE,
                      stringsAsFactors = FALSE)

  mat <- run_plage(expr, corum)
  mat <- round(mat, 5)
  write.csv(mat, file.path(out_dir, "plage.csv"), quote = FALSE)
  cat(sprintf("Saved plage: %d x %d\n", nrow(mat), ncol(mat)))
}
