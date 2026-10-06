
library(vegan)

BASE       <- "/proj/c.zihao/work3/01protein/01random"
PROT_DIR   <- "/proj/c.zihao/work3/00data/cptacT/protein"
SET_DIR    <- file.path(BASE, "01_random_sets")                          
SCORE_DIR  <- file.path(BASE, "02_scores")                                        
OUT_FILE   <- file.path(BASE, "results", "concord.csv")

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHODS <- c("netcomplex", "mean", "min", "zscore", "netmean",
             "perpas_adapted", "gsva", "ssgsea", "plage")

N_PC             <- 3
CV_TOP_FRACTIONS <- c(0.05, 0.10, 0.15, 0.20, 0.25)
N_PERM           <- 50
set.seed(1)

N_REPS <- length(list.files(SET_DIR, pattern = "^rep[0-9]+\\.csv$"))
if (N_REPS == 0) stop("no random sets found in ", SET_DIR, "; run step 01 first")
SETS <- c("real", sprintf("rep%03d", 1:N_REPS))
cat("random sets found:", N_REPS, "\n")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)


score_file <- function(method, cancer, set) {
  if (method == "netcomplex") {
    file.path(SCORE_DIR, set, "netcomplex", "complex_score", cancer,
              paste0(cancer, "_netcomplex_complex_score.csv"))
  } else {
    file.path(SCORE_DIR, set, method, cancer, paste0(method, ".csv"))
  }
}



clean_feature_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  mat <- mat[apply(mat, 1, function(x) sd(x) > 0), , drop = FALSE]
  mat
}

pca_embed <- function(mat, n_pc) {
  mat <- clean_feature_matrix(mat)
  if (nrow(mat) < n_pc || ncol(mat) <= n_pc) return(NULL)
  pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  if (ncol(pc$x) < n_pc) return(NULL)
  pc$x[, 1:n_pc, drop = FALSE]
}

select_top_cv <- function(mat, fraction) {
  mat <- clean_feature_matrix(abs(as.matrix(mat)))
  cvs <- apply(mat, 1, sd) / rowMeans(mat)
  cvs[!is.finite(cvs)] <- NA_real_
  valid <- which(!is.na(cvs))
  if (length(valid) < N_PC) return(NULL)

  n_keep <- max(1L, floor(length(valid) * fraction))
  if (n_keep < N_PC) return(NULL)
  keep <- names(sort(cvs[valid], decreasing = TRUE))[1:n_keep]
  mat[keep, , drop = FALSE]
}

compare_embeddings <- function(prot_embed, score_embed) {
  common <- intersect(rownames(prot_embed), rownames(score_embed))
  X <- prot_embed[common, , drop = FALSE]
  Y <- score_embed[common, , drop = FALSE]
  k <- min(ncol(X), ncol(Y))
  X <- X[, 1:k, drop = FALSE]
  Y <- Y[, 1:k, drop = FALSE]

  d_prot  <- dist(X, method = "euclidean")
  d_score <- dist(Y, method = "euclidean")

  spearman_r <- suppressWarnings(cor(as.vector(d_score), as.vector(d_prot), method = "spearman"))
  mt <- mantel(d_score, d_prot, method = "spearman", permutations = N_PERM)
  pt <- protest(X, Y, symmetric = TRUE, permutations = N_PERM)

  c(spearman_r      = as.numeric(spearman_r),
    mantel_r        = as.numeric(mt$statistic),
    mantel_pval     = as.numeric(mt$signif),
    procrustes_r    = as.numeric(pt$t0),
    procrustes_pval = as.numeric(pt$signif))
}


rows <- list()

for (method in METHODS) {
  cat("\n####", method, "####\n")

  for (cancer in CANCERS) {
    prot <- read.csv(file.path(PROT_DIR, paste0(cancer, "_proteomics.csv")),
                     row.names = 1, check.names = FALSE)
    prot[prot == 0] <- NA

    scores <- list()
    for (set in SETS) {
      f <- score_file(method, cancer, set)
      if (!file.exists(f)) stop("missing score file: ", f)
      scores[[set]] <- read.csv(f, row.names = 1, check.names = FALSE)
    }

    samples   <- colnames(prot)
    complexes <- rownames(scores[["real"]])
    for (set in SETS) {
      samples   <- intersect(samples, colnames(scores[[set]]))
      complexes <- intersect(complexes, rownames(scores[[set]]))
    }
    samples <- sort(samples)
    if (length(samples) < N_PC + 1 || length(complexes) < N_PC) {
      cat("[", cancer, "] skipped: too few samples or complexes\n")
      next
    }

    prot_use <- clean_feature_matrix(prot[, samples, drop = FALSE])
    if (nrow(prot_use) < N_PC) next
    prot_embed <- pca_embed(prot_use, N_PC)
    if (is.null(prot_embed)) next

    cat(sprintf("[%s] samples=%d proteins=%d complexes=%d\n",
                cancer, length(samples), nrow(prot_use), length(complexes)))

    for (cv_fraction in CV_TOP_FRACTIONS) {
      for (set in SETS) {
        mat <- scores[[set]][complexes, samples, drop = FALSE]
        top <- select_top_cv(mat, cv_fraction)
        if (is.null(top)) next
        score_embed <- pca_embed(top, N_PC)
        if (is.null(score_embed)) next

        stat <- compare_embeddings(prot_embed, score_embed)
        rows[[length(rows) + 1]] <- data.frame(
          target_method    = method,
          cancer           = cancer,
          method           = set,
          annotation       = if (set == "real") "real" else "random",
          cv_top_fraction  = cv_fraction,
          n_pc             = N_PC,
          n_samples        = length(samples),
          n_features       = nrow(top),
          n_protein_genes  = nrow(prot_use),
          n_common_complex = length(complexes),
          spearman_r       = round(stat[["spearman_r"]], 5),
          mantel_r         = round(stat[["mantel_r"]], 5),
          mantel_pval      = round(stat[["mantel_pval"]], 5),
          procrustes_r     = round(stat[["procrustes_r"]], 5),
          procrustes_pval  = round(stat[["procrustes_pval"]], 5)
        )
      }
    }
  }
}

res <- do.call(rbind, rows)
write.csv(res, OUT_FILE, row.names = FALSE)
cat("\n[done]", nrow(res), "rows ->", OUT_FILE, "\n")
