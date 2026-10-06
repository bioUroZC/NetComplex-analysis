
rm(list = ls())
library(data.table)

OUT_DIR <- "/proj/c.zihao/work3/02TCGAnt/02stringTest/02analysis"
ROOT <- "/proj/c.zihao/work3/02TCGAnt/02stringTest"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANCERS <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KICH", "KIRC", "KIRP",
             "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "UCEC")

METHODS <- data.frame(
  method = c("netcomplex", "mean", "min", "zscore", "ssgsea", "gsva", "plage", "netmean", "perpas"),
  label = c("netcomplex", "mean", "min", "z-score", "ssGSEA", "GSVA", "PLAGE", "netmean", "PerPAS"),
  score_dir = c("netcomplex/complex_score", "bench/mean", "bench/min", "bench/zscore",
                "bench/ssgsea", "bench/gsva", "bench/plage", "bench/netmean",
                "bench/perpas_adapted"),
  filename = c("%s_netcomplex_complex_score.csv", "mean.csv", "min.csv", "zscore.csv",
               "ssgsea.csv", "gsva.csv", "plage.csv", "netmean.csv", "perpas_adapted.csv"),
  stringsAsFactors = FALSE
)

TOP_FRACS <- c(0.05, 0.10, 0.15, 0.20, 0.25)
MIN_FEATURES <- 3L
MIN_SAMPLES <- MIN_FEATURES + 1L
K <- 2
NSTART <- 100
SEED <- 1

roc_auc <- function(labels, scores) {
  pos <- sum(labels == 1L); neg <- sum(labels == 0L)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[labels == 1L]) - pos * (pos + 1L) / 2) / (pos * neg)
}

pr_auc <- function(labels, scores) {
  pos <- sum(labels == 1L)
  y <- labels[order(scores, decreasing = TRUE)]
  precision <- cumsum(y == 1L) / seq_along(y)
  recall <- cumsum(y == 1L) / pos
  sum((recall - c(0, head(recall, -1))) * (precision + c(1, head(precision, -1))) / 2)
}

score_path <- function(cancer, score_dir, filename_pattern) {
  file.path(ROOT, score_dir, cancer, sub("%s", cancer, filename_pattern, fixed = TRUE))
}

load_score <- function(cancer, i) {
  path <- score_path(cancer, METHODS$score_dir[i], METHODS$filename[i])
  if (!file.exists(path)) stop("Missing score matrix: ", path)
  x <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(x) <- "double"
  abs(x)
}

metric_rows <- list()
for (cancer in CANCERS) {
  cat("\n=== ", cancer, " ===\n", sep = "")

  score_list <- lapply(seq_len(nrow(METHODS)), function(i) load_score(cancer, i))
  names(score_list) <- METHODS$label
  reference_score <- score_list[[1]]
  same_samples <- vapply(score_list, function(x) identical(colnames(x), colnames(reference_score)),
                         logical(1))
  same_complexes <- vapply(score_list, function(x) identical(rownames(x), rownames(reference_score)),
                           logical(1))
  if (!all(same_samples) || !all(same_complexes)) {
    mismatched_methods <- names(score_list)[!(same_samples & same_complexes)]
    stop(sprintf("[%s] sample or complex universe differs across methods: %s",
                 cancer, paste(mismatched_methods, collapse = ", ")))
  }
  common_samples <- colnames(reference_score)
  common_complexes <- rownames(reference_score)

  if (length(common_samples) < MIN_SAMPLES || length(common_complexes) < MIN_FEATURES) {
    cat("  skip: insufficient samples or complexes\n")
    next
  }

  for (i in seq_len(nrow(METHODS))) {
    method_label <- METHODS$label[i]
    x <- score_list[[method_label]][common_complexes, common_samples, drop = FALSE]

    valid <- apply(x, 1, function(v) all(is.finite(v)) && sd(v) > 0)
    x <- x[valid, , drop = FALSE]
    feature_cv <- apply(x, 1, sd) / rowMeans(x)
    feature_cv <- feature_cv[is.finite(feature_cv)]
    x <- x[names(feature_cv), , drop = FALSE]
    if (length(feature_cv) < MIN_FEATURES || ncol(x) < MIN_SAMPLES) next

    sample_ids <- colnames(x)
    groups <- ifelse(grepl("11A$", sample_ids), "Normal", "Tumor")
    labels <- as.integer(groups == "Tumor")
    if (length(unique(labels)) != 2L) next
    tumor_ids <- sub("_[^_]+$", "", sample_ids[labels == 1L])
    normal_ids <- sub("_[^_]+$", "", sample_ids[labels == 0L])
    paired_ids <- intersect(tumor_ids, normal_ids)

    for (top_frac in TOP_FRACS) {
      n_top <- max(1L, floor(length(feature_cv) * top_frac))
      if (n_top < MIN_FEATURES) next
      keep <- names(sort(feature_cv, decreasing = TRUE))[seq_len(n_top)]
      feature_matrix <- t(x[keep, , drop = FALSE])
      z <- scale(feature_matrix)

      set.seed(SEED)
      km <- kmeans(z, centers = K, nstart = NSTART)
      agreement <- max(mean(ifelse(km$cluster == 1L, "Normal", "Tumor") == groups),
                       mean(ifelse(km$cluster == 2L, "Normal", "Tumor") == groups))

      pc <- prcomp(z, center = FALSE, scale. = FALSE)
      pc1 <- pc$x[, 1]
      if (mean(pc1[labels == 1L]) < mean(pc1[labels == 0L])) pc1 <- -pc1
      t_idx <- match(paired_ids, tumor_ids)
      n_idx <- match(paired_ids, normal_ids)
      paired_concordance <- if (length(paired_ids) == 0L) NA_real_ else
        mean(pc1[labels == 1L][t_idx] > pc1[labels == 0L][n_idx])

      metric_rows[[length(metric_rows) + 1L]] <- data.frame(
        cancer = cancer, method = method_label,
        feature_fraction = paste0(as.integer(top_frac * 100), "%"),
        n_samples = ncol(x), n_variable_features = length(feature_cv),
        n_selected_features = length(keep),
        cluster_agreement = round(agreement, 5),
        auc = round(roc_auc(labels, pc1), 5),
        aupr = round(pr_auc(labels, pc1), 5),
        paired_concordance = round(paired_concordance, 5),
        pc1_var = round(pc$sdev[1]^2 / sum(pc$sdev^2), 5),
        stringsAsFactors = FALSE
      )
    }
  }
}

metrics <- rbindlist(metric_rows)
fwrite(metrics, file.path(OUT_DIR, "cluster_metrics.csv"))
cat("Saved: ", file.path(OUT_DIR, "cluster_metrics.csv"), "\n", sep = "")
