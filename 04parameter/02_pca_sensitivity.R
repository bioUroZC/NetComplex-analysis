
library(vegan)

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
ALPHAS <- c(0.1, 0.3, 0.5, 0.7, 0.9)
N_PC <- 3L
CV_FRACTIONS <- c(0.05, 0.10, 0.15, 0.20, 0.25)
N_PERM <- 50
set.seed(1)

ROOT <- "/proj/c.zihao/work3/04parameter"
SCORE_ROOT <- file.path(ROOT, "results/complex_score")
PROT_DIR <- "/proj/c.zihao/work3/00data/cptacT/protein"
OUT_DIR <- file.path(ROOT, "results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

clean_feature_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  mat[apply(mat, 1, function(x) sd(x) > 0), , drop = FALSE]
}

select_top_cv <- function(mat, fraction) {
  mat <- abs(as.matrix(mat))
  mat <- clean_feature_matrix(mat)
  means <- rowMeans(mat)
  cvs <- apply(mat, 1, sd) / means
  cvs[!is.finite(cvs)] <- NA_real_
  valid <- which(!is.na(cvs))
  if (length(valid) < N_PC) return(NULL)
  n_keep <- max(1L, floor(length(valid) * fraction))
  if (n_keep < N_PC) return(NULL)
  keep <- names(sort(cvs[valid], decreasing = TRUE))[seq_len(n_keep)]
  mat[keep, , drop = FALSE]
}

pca_embed <- function(mat, n_pc) {
  mat <- clean_feature_matrix(mat)
  if (nrow(mat) < n_pc || ncol(mat) <= n_pc) return(NULL)
  fit <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  if (ncol(fit$x) < n_pc) return(NULL)
  fit$x[, seq_len(n_pc), drop = FALSE]
}

compare_embeddings <- function(prot_embed, score_embed) {
  common <- intersect(rownames(prot_embed), rownames(score_embed))
  X <- prot_embed[common, , drop = FALSE]
  Y <- score_embed[common, , drop = FALSE]
  d_prot <- dist(X)
  d_score <- dist(Y)
  mt <- mantel(d_score, d_prot, method = "spearman", permutations = N_PERM)
  pt <- protest(X, Y, symmetric = TRUE, permutations = N_PERM)
  c(mantel_r = as.numeric(mt$statistic), mantel_pval = as.numeric(mt$signif),
    procrustes_r = as.numeric(pt$t0), procrustes_pval = as.numeric(pt$signif))
}

rows <- list()
for (alpha in ALPHAS) {
  for (cancer in CANCERS) {
    score_file <- file.path(SCORE_ROOT, sprintf("alpha_%.1f", alpha), cancer,
                            sprintf("%s_complex_score.csv", cancer))
    if (!file.exists(score_file)) stop("Missing parameter score matrix: ", score_file)
    scores <- read.csv(score_file, row.names = 1, check.names = FALSE)
    prot <- read.csv(file.path(PROT_DIR, paste0(cancer, "_proteomics.csv")),
                     row.names = 1, check.names = FALSE)
    prot[prot == 0] <- NA
    common_samples <- sort(intersect(colnames(scores), colnames(prot)))
    scores <- scores[, common_samples, drop = FALSE]
    prot <- clean_feature_matrix(prot[, common_samples, drop = FALSE])

    prot_embed <- pca_embed(prot, N_PC)
    if (is.null(prot_embed)) next

    for (cv_fraction in CV_FRACTIONS) {
      selected <- select_top_cv(scores, cv_fraction)
      if (is.null(selected)) next
      score_embed <- pca_embed(selected, N_PC)
      if (is.null(score_embed)) next
      stat <- compare_embeddings(prot_embed, score_embed)
      rows[[length(rows) + 1L]] <- data.frame(
        cancer = cancer, alpha = alpha, n_pc = N_PC,
        cv_top_fraction = cv_fraction,
        n_samples = length(common_samples), n_features = nrow(selected),
        mantel_r = round(stat["mantel_r"], 4), mantel_pval = stat["mantel_pval"],
        procrustes_r = round(stat["procrustes_r"], 4),
        procrustes_pval = stat["procrustes_pval"], stringsAsFactors = FALSE
      )
      cat(sprintf("alpha=%.1f %s PC%d topCV=%g%% Mantel=%.3f Proc=%.3f\n",
                  alpha, cancer, N_PC, 100 * cv_fraction,
                  stat["mantel_r"], stat["procrustes_r"]))
    }
  }
}

res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, "pca_alpha_sensitivity.csv"), row.names = FALSE)
cat("Output ->", file.path(OUT_DIR, "pca_alpha_sensitivity.csv"), "\n")
