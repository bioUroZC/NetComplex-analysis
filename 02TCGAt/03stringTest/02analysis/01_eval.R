rm(list = ls())

library(data.table)
library(cluster)
library(mclust)

SCORE_ROOT <- "/proj/c.zihao/work3/03TCGAt/03stringTest"
OUT_DIR    <- file.path(SCORE_ROOT, "02analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANCERS <- c(
  "BLCA", "BRCA", "CRC", "GBM", "KIRC", "LGG", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PRAD", "SKCM", "UVM"
)

METHODS <- data.frame(
  method = c(
    "netComplex", "mean", "min", "zscore",
    "ssgsea", "gsva", "plage", "netmean", "perpas"
  ),
  label = c(
    "NetComplex", "mean", "min", "z-score",
    "ssGSEA", "GSVA", "PLAGE", "netmean", "PerPAS"
  ),
  path_pat = c(
    "netcomplex/complex_score/%s/%s_netcomplex_complex_score.csv",
    "bench/mean/%s/mean.csv",
    "bench/min/%s/min.csv",
    "bench/zscore/%s/zscore.csv",
    "bench/ssgsea/%s/ssgsea.csv",
    "bench/gsva/%s/gsva.csv",
    "bench/plage/%s/plage.csv",
    "bench/netmean/%s/netmean.csv",
    "bench/perpas_adapted/%s/perpas_adapted.csv"
  ),
  stringsAsFactors = FALSE
)

TOP_FRACS  <- c(0.05, 0.10, 0.15, 0.20, 0.25)
N_PCS      <- 50L
NSTART     <- 50L
SEED       <- 1L

get_path <- function(cancer, path_pat) {
  file.path(SCORE_ROOT, gsub("%s", cancer, path_pat, fixed = TRUE))
}

rows <- list()
rows_silhouette_by_cancer <- list()

for (i in seq_len(nrow(METHODS))) {
  method_label <- METHODS$label[i]
  cat(sprintf("\n=== %s ===\n", method_label))

  mats <- list()
  labels <- character(0)

  for (cancer in CANCERS) {
    path <- get_path(cancer, METHODS$path_pat[i])
    if (!file.exists(path)) {
      stop("Missing score file: ", path)
    }
    mat <- as.matrix(fread(path, header = TRUE), rownames = 1)
    mat <- abs(mat)
    storage.mode(mat) <- "double"
    mat[!is.finite(mat)] <- 0
    mats[[cancer]] <- mat
    labels <- c(labels, rep(cancer, ncol(mat)))
  }

  common_rows <- Reduce(intersect, lapply(mats, rownames))
  pooled <- do.call(cbind, lapply(mats, function(m) m[common_rows, , drop = FALSE]))
  cat(sprintf("  Pooled: %d complexes x %d samples (%d cancers)\n",
              nrow(pooled), ncol(pooled), length(mats)))

  label_int <- as.integer(factor(labels, levels = CANCERS))
  k <- length(unique(label_int))

  means <- rowMeans(pooled)
  sds   <- apply(pooled, 1L, sd)
  cv    <- sds / means
  cv    <- cv[is.finite(cv)]

  for (top_frac in TOP_FRACS) {
    frac_label <- paste0(as.integer(top_frac * 100L), "%")
    top_n <- max(k, floor(length(cv) * top_frac))
    keep  <- names(sort(cv, decreasing = TRUE))[seq_len(min(top_n, length(cv)))]

    X_raw <- t(pooled[keep, , drop = FALSE])
    nonconst <- apply(X_raw, 2L, sd) > 0
    X_raw <- X_raw[, nonconst, drop = FALSE]
    if (ncol(X_raw) == 0L) next

    X <- scale(X_raw)

    n_pc <- min(N_PCS, nrow(X) - 1L, ncol(X) - 1L)
    pca  <- prcomp(X, center = FALSE, scale. = FALSE, rank. = n_pc)

    d <- dist(pca$x)
    sil_obj <- silhouette(label_int, d)
    sil <- mean(sil_obj[, 3L])

    set.seed(SEED)
    km <- suppressWarnings(
      kmeans(pca$x, centers = k, nstart = NSTART, iter.max = 200L)
    )
    if (!is.null(km$ifault) && km$ifault == 4L) {
      cat("    Hartigan-Wong k-means did not converge; retrying with Lloyd k-means\n")
      set.seed(SEED)
      km <- kmeans(pca$x, centers = k, nstart = NSTART,
                   iter.max = 200L, algorithm = "Lloyd")
    }
    ari <- adjustedRandIndex(km$cluster, label_int)

    cat(sprintf("  %s: n_feat=%d  sil=%.3f  ARI=%.3f\n",
                frac_label, length(keep), sil, ari))

    rows[[length(rows) + 1L]] <- data.frame(
      method           = method_label,
      feature_fraction = frac_label,
      n_cancers        = k,
      n_samples        = ncol(pooled),
      n_features       = length(keep),
      silhouette       = round(sil, 4),
      ARI              = round(ari, 4),
      stringsAsFactors = FALSE
    )

    silhouette_by_cancer <- data.table(
      method = method_label,
      feature_fraction = frac_label,
      cancer = labels,
      silhouette = sil_obj[, 3L]
    )[, .(silhouette = mean(silhouette, na.rm = TRUE)),
        by = .(method, feature_fraction, cancer)]
    rows_silhouette_by_cancer[[length(rows_silhouette_by_cancer) + 1L]] <- silhouette_by_cancer

  }
}

result_dt <- rbindlist(rows)
out_path  <- file.path(OUT_DIR, "pancancer_results.csv")
fwrite(result_dt, out_path)
cat("\nSaved:", out_path, "\n")
cat("Total rows:", nrow(result_dt), "\n")

silhouette_by_cancer_dt <- rbindlist(rows_silhouette_by_cancer, use.names = TRUE)
silhouette_out_path <- file.path(OUT_DIR, "pancancer_silhouette_by_cancer.csv")
fwrite(silhouette_by_cancer_dt, silhouette_out_path)
cat("Saved:", silhouette_out_path, "\n")

rank_target <- result_dt[feature_fraction == "5%"][order(-ARI, -silhouette)]
cat("\n=== Global ranking at 5% (by ARI) ===\n")
print(rank_target[, .(method, silhouette, ARI)])
