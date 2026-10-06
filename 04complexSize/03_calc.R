
library(vegan)

SCORE_ROOT      <- "/proj/c.zihao/work3/01protein/01stringTest/bench"
PROT_DIR        <- "/proj/c.zihao/work3/00data/cptacT/protein"
OUT_DIR         <- "/proj/c.zihao/work3/04complexSize"
NETCOMPLEX_ROOT <- "/proj/c.zihao/work3/01protein/01stringTest/netcomplex"
COMPLEX_DB      <- "/proj/c.zihao/work3/00data/corum/complexData.csv"

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHOD_GROUPS <- list(
  benchmark = c(
    netcomplex = "netcomplex", mean = "mean", min = "min", zscore = "zscore",
    netmean = "netmean", perpas = "perpas_adapted", gsva = "gsva",
    ssgsea = "ssgsea", plage = "plage"
  )
)

SIZE_BIN_LEVELS <- c("1", "2", "3", "4", "5-7", "8-10", "11+")

N_PC <- 3
CV_TOP_FRACTIONS <- c(0.05, 0.10, 0.15, 0.20, 0.25)
N_PERM  <- 50
set.seed(1)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cdb <- read.csv(COMPLEX_DB, stringsAsFactors = FALSE)
complex_size <- vapply(strsplit(cdb$Genes, ";"),
                       function(x) sum(nzchar(trimws(x))), integer(1))
size_bin <- ifelse(complex_size <= 4L, as.character(complex_size),
            ifelse(complex_size <= 7L, "5-7",
            ifelse(complex_size <= 10L, "8-10", "11+")))
BIN_MEMBERS <- split(cdb$Complex, factor(size_bin, levels = SIZE_BIN_LEVELS))

cat("=== complex size bins ===\n")
for (lv in SIZE_BIN_LEVELS) {
  cat(sprintf("  size %-5s %d complexes\n", lv, length(BIN_MEMBERS[[lv]])))
}

tumor_ids <- function(cols) cols

load_score <- function(cancer, method) {
  if (method == "netcomplex") {
    file <- file.path(NETCOMPLEX_ROOT, "complex_score", cancer,
                      sprintf("%s_netcomplex_complex_score.csv", cancer))
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

    prot_use <- clean_feature_matrix(prot[, common_samples, drop = FALSE])
    if (nrow(prot_use) < N_PC) next

    prot_embed <- pca_embed(prot_use, N_PC)
    if (is.null(prot_embed)) next

    cat(sprintf("    samples=%d, protein_genes=%d, complexes_scored=%d\n",
                length(common_samples), nrow(prot_use), length(common_complexes)))

    for (bin in SIZE_BIN_LEVELS) {
      bin_complexes <- intersect(common_complexes, BIN_MEMBERS[[bin]])
      if (length(bin_complexes) < N_PC) {
        cat(sprintf("    size %-5s skip: only %d complexes\n", bin, length(bin_complexes)))
        next
      }

      reps <- lapply(score_list, function(x) x[bin_complexes, common_samples, drop = FALSE])

      for (cv_fraction in CV_TOP_FRACTIONS) {
        for (m in names(reps)) {
          repr_mat <- select_top_cv(reps[[m]], cv_fraction)
          if (is.null(repr_mat)) next
          repr_embed <- pca_embed(repr_mat, N_PC)
          if (is.null(repr_embed)) next

          stat <- compare_embeddings(prot_embed, repr_embed)
          all_res[[length(all_res) + 1]] <- data.frame(
            comparison_group = comparison_group,
            size_bin = bin,
            cancer = cancer,
            method = m,
            cv_top_fraction = cv_fraction,
            n_pc = N_PC,
            n_samples = length(common_samples),
            n_features = nrow(repr_mat),
            n_protein_genes = nrow(prot_use),
            n_common_complex = length(bin_complexes),
            spearman_r = round(stat["spearman_r"], 5),
            mantel_r = round(stat["mantel_r"], 5),
            mantel_pval = round(stat["mantel_pval"], 5),
            procrustes_r = round(stat["procrustes_r"], 5),
            procrustes_pval = round(stat["procrustes_pval"], 5),
            stringsAsFactors = FALSE
          )
        }
      }
      cat(sprintf("    size %-5s complexes=%4d  topCV=%d%%..%d%% done\n", bin,
                  length(bin_complexes), as.integer(100 * min(CV_TOP_FRACTIONS)),
                  as.integer(100 * max(CV_TOP_FRACTIONS))))
    }
  }
}

res <- do.call(rbind, all_res)
res$size_bin <- factor(res$size_bin, levels = SIZE_BIN_LEVELS)
res <- res[order(res$size_bin, res$cancer, res$cv_top_fraction), ]
write.csv(res, file.path(OUT_DIR, "concord.csv"), row.names = FALSE)

cat("\n=== PC3 cross-cancer mean Mantel r by size bin (averaged over CV fractions) ===\n")
print(round(with(res, tapply(mantel_r, list(size_bin, method), mean, na.rm = TRUE)), 3))
cat("\n=== PC3 cross-cancer mean Procrustes r by size bin (averaged over CV fractions) ===\n")
print(round(with(res, tapply(procrustes_r, list(size_bin, method), mean, na.rm = TRUE)), 3))

cat("\nOutput ->", file.path(OUT_DIR, "concord.csv"), "\n")
