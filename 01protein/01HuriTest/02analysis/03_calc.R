
library(vegan)

SCORE_ROOT        <- "/proj/c.zihao/work3/01protein/01HuriTest/bench"
PROT_DIR          <- "/proj/c.zihao/work3/00data/cptacT/protein"
OUT_DIR           <- "/proj/c.zihao/work3/01protein/01HuriTest/02analysis"
NETCOMPLEX_ROOT <- "/proj/c.zihao/work3/01protein/01HuriTest/netcomplex"

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHOD_GROUPS <- list(
  benchmark = c(
    netcomplex = "netcomplex", mean = "mean", min = "min", zscore = "zscore",
    netmean = "netmean", perpas = "perpas_adapted", gsva = "gsva",
    ssgsea = "ssgsea", plage = "plage"
  ),
  ablation = c(
    netcomplex = "netcomplex", noRWR = "noRWR", noRank = "noRank",
    degRandPPI = "degRandPPI", permNodeValues = "permNodeValues"
  )
)

N_PC <- 3
CV_TOP_FRACTIONS <- c(0.05, 0.10, 0.15, 0.20, 0.25)
N_PERM  <- 50
set.seed(1)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

tumor_ids <- function(cols) cols

load_score <- function(cancer, method) {
  if (method == "netcomplex") {
    file <- file.path(NETCOMPLEX_ROOT, "complex_score", cancer,
                      sprintf("%s_netcomplex_complex_score.csv", cancer))
  } else if (method == "noRWR") {
    file <- file.path("/proj/c.zihao/work3/01protein/01HuriTest/ablation", "noRWR", cancer,
                      sprintf("%s_noRWR_complex_score.csv", cancer))
  } else if (method %in% c("noRank", "permNodeValues", "degRandPPI")) {
    file <- file.path("/proj/c.zihao/work3/01protein/01HuriTest/ablation", method,
                      "complex_score", cancer,
                      sprintf("%s_%s_complex_score.csv", cancer, method))
  } else {
    file <- file.path(SCORE_ROOT, method, cancer, sprintf("%s.csv", method))
  }
  if (!file.exists(file)) return(NULL)
  read.csv(file, row.names = 1, check.names = FALSE)
}

clean_feature_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  mat <- mat[apply(mat, 1, function(x) sd(x) > 0), , drop = FALSE]
  mat
}

pca_embed <- function(mat_feat_by_sample, n_pc) {
  mat <- clean_feature_matrix(mat_feat_by_sample)
  if (nrow(mat) < n_pc || ncol(mat) <= n_pc) return(NULL)
  pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  if (ncol(pc$x) < n_pc) return(NULL)
  pc$x[, seq_len(n_pc), drop = FALSE]
}

select_top_cv <- function(mat_feat_by_sample, fraction) {
  mat <- abs(as.matrix(mat_feat_by_sample))
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

compare_embeddings <- function(prot_embed, repr_embed) {
  common <- intersect(rownames(prot_embed), rownames(repr_embed))
  X <- prot_embed[common, , drop = FALSE]
  Y <- repr_embed[common, , drop = FALSE]
  k <- min(ncol(X), ncol(Y))
  X <- X[, seq_len(k), drop = FALSE]
  Y <- Y[, seq_len(k), drop = FALSE]

  d_prot <- dist(X, method = "euclidean")
  d_repr <- dist(Y, method = "euclidean")

  spearman_r <- suppressWarnings(cor(as.vector(d_repr), as.vector(d_prot),
                                     method = "spearman"))
  mt <- mantel(d_repr, d_prot, method = "spearman", permutations = N_PERM)
  pt <- protest(X, Y, symmetric = TRUE, permutations = N_PERM)

  c(spearman_r = as.numeric(spearman_r),
    mantel_r = as.numeric(mt$statistic),
    mantel_pval = as.numeric(mt$signif),
    procrustes_r = as.numeric(pt$t0),
    procrustes_pval = as.numeric(pt$signif))
}

all_res <- list()

for (comparison_group in names(METHOD_GROUPS)) {
  group_sources <- METHOD_GROUPS[[comparison_group]]
  group_methods <- names(group_sources)
  cat("\n===", comparison_group, "===\n")
  for (cancer in CANCERS) {
    cat("[", cancer, "] ...\n")

    prot <- read.csv(file.path(PROT_DIR, paste0(cancer, "_proteomics.csv")),
                     row.names = 1, check.names = FALSE)
    prot[prot == 0] <- NA
    score_list <- lapply(group_sources, function(source) load_score(cancer, source))
    names(score_list) <- group_methods
    missing_methods <- names(score_list)[vapply(score_list, is.null, logical(1))]
    if (length(missing_methods) > 0) {
      stop(sprintf("[%s/%s] missing score matrices: %s", comparison_group, cancer,
                   paste(missing_methods, collapse = ", ")))
    }

    common_samples <- Reduce(intersect, c(
      list(tumor_ids(colnames(prot))),
      lapply(score_list, function(x) tumor_ids(colnames(x)))
    ))
    common_samples <- sort(common_samples)
    if (length(common_samples) < N_PC + 1) {
      cat(sprintf("    skip: only %d common tumor samples\n", length(common_samples)))
      next
    }

    common_complexes <- Reduce(intersect, lapply(score_list, rownames))
    if (length(common_complexes) < N_PC) {
      cat(sprintf("    skip: only %d common complexes\n", length(common_complexes)))
      next
    }

    prot_use <- clean_feature_matrix(prot[, common_samples, drop = FALSE])
    if (nrow(prot_use) < N_PC) next

    reps <- lapply(score_list, function(x) x[common_complexes, common_samples, drop = FALSE])

    cat(sprintf("    samples=%d, protein_genes=%d, common_complex=%d\n",
                length(common_samples), nrow(prot_use), length(common_complexes)))

    prot_embed <- pca_embed(prot_use, N_PC)
    if (is.null(prot_embed)) next

    for (cv_fraction in CV_TOP_FRACTIONS) {
      for (m in names(reps)) {
        repr_mat <- select_top_cv(reps[[m]], cv_fraction)
        if (is.null(repr_mat)) next
        repr_embed <- pca_embed(repr_mat, N_PC)
        if (is.null(repr_embed)) next

        stat <- compare_embeddings(prot_embed, repr_embed)
        all_res[[length(all_res) + 1]] <- data.frame(
          comparison_group = comparison_group,
          cancer = cancer,
          method = m,
          cv_top_fraction = cv_fraction,
          n_pc = N_PC,
          n_samples = length(common_samples),
          n_features = nrow(repr_mat),
          n_protein_genes = nrow(prot_use),
          n_common_complex = length(common_complexes),
          spearman_r = round(stat["spearman_r"], 5),
          mantel_r = round(stat["mantel_r"], 5),
          mantel_pval = round(stat["mantel_pval"], 5),
          procrustes_r = round(stat["procrustes_r"], 5),
          procrustes_pval = round(stat["procrustes_pval"], 5),
          stringsAsFactors = FALSE
        )
        cat(sprintf("      topCV=%-3g%% %-15s n=%d spearman=%.3f mantel=%.3f proc=%.3f\n",
                    100 * cv_fraction, m, nrow(repr_mat), stat["spearman_r"],
                    stat["mantel_r"], stat["procrustes_r"]))
      }
    }
  }
}

res <- do.call(rbind, all_res)
write.csv(res, file.path(OUT_DIR, "concord.csv"), row.names = FALSE)

for (comparison_group in names(METHOD_GROUPS)) {
  group_res <- res[res$comparison_group == comparison_group, ]
  cat("\n===", comparison_group, ": PC3 cross-cancer median Mantel r by top-CV fraction ===\n")
  print(with(group_res, tapply(mantel_r, list(cv_top_fraction, method), median, na.rm = TRUE)))
  cat("\n===", comparison_group, ": PC3 cross-cancer median Procrustes r by top-CV fraction ===\n")
  print(with(group_res, tapply(procrustes_r, list(cv_top_fraction, method), median, na.rm = TRUE)))
}

cat("\nOutput ->", file.path(OUT_DIR, "concord.csv"), "\n")
