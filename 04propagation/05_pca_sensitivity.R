
library(vegan)
library(data.table)

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC",
             "HNSCC")
K_VALUES <- c(0, 1, 2, 3, 5, 8, 12, 20, 50)
SEEDS <- 1:20
N_PC <- 3L
CV_FRACTIONS <- c(0.05, 0.10, 0.15, 0.20, 0.25)
N_PERM <- 0

ROOT <- "/proj/c.zihao/work3/04propagation"
PROT_DIR <- "/proj/c.zihao/work3/00data/cptacT/protein"
OUT_DIR <- file.path(ROOT, "results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

score_path <- function(network, seed, k, cancer) {
  if (network == "netcomplex") {
    file.path(ROOT, "netcomplex", sprintf("K_%d", k), "complex_score", cancer,
              sprintf("%s_netcomplex_complex_score.csv", cancer))
  } else if (network == "noRestart") {
    file.path(ROOT, "ablation/noRestart", sprintf("K_%d", k), "complex_score",
              cancer, sprintf("%s_noRestart_complex_score.csv", cancer))
  } else {
    file.path(ROOT, "ablation/degRandPPI", sprintf("seed_%03d", seed),
              sprintf("K_%d", k), "complex_score", cancer,
              sprintf("%s_degRandPPI_complex_score.csv", cancer))
  }
}

read_scores <- function(file) {
  dt <- fread(file, data.table = FALSE)
  mat <- as.matrix(dt[, -1, drop = FALSE])
  rownames(mat) <- dt[[1]]
  storage.mode(mat) <- "double"
  mat
}

clean_feature_matrix <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  mat[apply(mat, 1, function(x) sd(x) > 0), , drop = FALSE]
}

cv_ordered <- function(mat) {
  mat <- abs(as.matrix(mat))
  mat <- clean_feature_matrix(mat)
  if (nrow(mat) == 0) return(NULL)
  cvs <- apply(mat, 1, sd) / rowMeans(mat)
  cvs[!is.finite(cvs)] <- NA_real_
  valid <- which(!is.na(cvs))
  if (length(valid) < N_PC) return(NULL)
  keep <- names(sort(cvs[valid], decreasing = TRUE))
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
  mt <- mantel(dist(Y), dist(X), method = "spearman", permutations = N_PERM)
  pt <- protest(X, Y, symmetric = TRUE, permutations = N_PERM)
  c(mantel_r = as.numeric(mt$statistic), procrustes_r = as.numeric(pt$t0))
}

prot_raw <- list()
for (cancer in CANCERS) {
  prot <- read.csv(file.path(PROT_DIR, paste0(cancer, "_proteomics.csv")),
                   row.names = 1, check.names = FALSE)
  prot[prot == 0] <- NA
  prot_raw[[cancer]] <- prot
}

evaluate_one <- function(network, seed, k, cancer, cache) {
  file <- score_path(network, seed, k, cancer)
  if (!file.exists(file)) stop("Missing score matrix: ", file)
  scores <- read_scores(file)

  prot <- prot_raw[[cancer]]
  common_samples <- sort(intersect(colnames(scores), colnames(prot)))
  scores <- scores[, common_samples, drop = FALSE]

  key <- paste(cancer, length(common_samples), sep = "_")
  if (is.null(cache[[key]])) {
    cache[[key]] <- pca_embed(
      clean_feature_matrix(as.matrix(prot[, common_samples, drop = FALSE])), N_PC)
  }
  prot_embed <- cache[[key]]
  if (is.null(prot_embed)) return(list(rows = NULL, cache = cache))

  ordered <- cv_ordered(scores)
  if (is.null(ordered)) return(list(rows = NULL, cache = cache))

  rows <- list()
  for (cv_fraction in CV_FRACTIONS) {
    n_keep <- max(1L, floor(nrow(ordered) * cv_fraction))
    if (n_keep < N_PC) next
    score_embed <- pca_embed(ordered[seq_len(n_keep), , drop = FALSE], N_PC)
    if (is.null(score_embed)) next
    stat <- compare_embeddings(prot_embed, score_embed)
    rows[[length(rows) + 1L]] <- data.frame(
      network = network, seed = seed, cancer = cancer, k = k, n_pc = N_PC,
      cv_top_fraction = cv_fraction, n_samples = length(common_samples),
      n_features = n_keep,
      mantel_r = round(stat["mantel_r"], 5),
      procrustes_r = round(stat["procrustes_r"], 5),
      stringsAsFactors = FALSE
    )
  }
  list(rows = rows, cache = cache)
}

cache <- list()
out <- list()

cat("=== real network ===\n")
for (k in K_VALUES) {
  for (cancer in CANCERS) {
    res <- evaluate_one("netcomplex", 0L, k, cancer, cache)
    cache <- res$cache
    out <- c(out, res$rows)
  }
  cat(sprintf("  K=%-2d done (%d rows)\n", k, length(out)))
}

cat("=== noRestart control (alpha = 0, pure diffusion) ===\n")
for (k in K_VALUES) {
  for (cancer in CANCERS) {
    res <- evaluate_one("noRestart", 0L, k, cancer, cache)
    cache <- res$cache
    out <- c(out, res$rows)
  }
  cat(sprintf("  K=%-2d done (%d rows)\n", k, length(out)))
}

cat(sprintf("=== degRandPPI null, %d seeds ===\n", length(SEEDS)))
for (seed in SEEDS) {
  for (k in K_VALUES) {
    for (cancer in CANCERS) {
      res <- evaluate_one("degRandPPI", seed, k, cancer, cache)
      cache <- res$cache
      out <- c(out, res$rows)
    }
  }
  cat(sprintf("  seed=%-3d done (%d rows total)\n", seed, length(out)))
}

res <- do.call(rbind, out)
fwrite(res, file.path(OUT_DIR, "pca_step_sensitivity.csv"))

cat("\n=== Cross-cancer mean Mantel r by K ===\n")
agg <- aggregate(mantel_r ~ k + network, data = res, FUN = mean, na.rm = TRUE)
print(reshape(agg, idvar = "k", timevar = "network", direction = "wide"), digits = 4)

cat("\n=== Cross-cancer mean Procrustes r by K ===\n")
agg <- aggregate(procrustes_r ~ k + network, data = res, FUN = mean, na.rm = TRUE)
print(reshape(agg, idvar = "k", timevar = "network", direction = "wide"), digits = 4)

cat("\nOutput ->", file.path(OUT_DIR, "pca_step_sensitivity.csv"), "\n")
