library(data.table)
library(survival)

CANCERS       <- c("GBM", "PDAC", "LUAD", "LUSC", "KIRC")
METHODS       <- c("NetComplex_rna", "NetComplex_prot",
                   "NetComplex_multi", "NetComplex_multilayer")
SCORE_BASE    <- "/proj/c.zihao/work3/07multi/netcomplex"
CLINICAL_BASE <- "/proj/c.zihao/work3/00data/cptacT"
OUT_DIR       <- "/proj/c.zihao/work3/07multi/02survival"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CV_CUTOFFS <- c("5%"=0.05, "10%"=0.10, "15%"=0.15, "20%"=0.20, "25%"=0.25)
N_STARTS   <- 50
set.seed(42)


load_clinical <- function(cancer) {
  path <- file.path(CLINICAL_BASE, cancer, paste0(cancer, "_clinical.csv"))
  dt <- fread(path, check.names = FALSE)
  id_col <- if ("Patient_ID" %in% names(dt)) "Patient_ID" else names(dt)[1]
  setnames(dt, id_col, "sample_id")
  dt[, sample_id := as.character(sample_id)]
  os_col <- grep("Overall survival.*days", names(dt), value = TRUE)
  os_col <- os_col[!grepl("collection", os_col)]
  st_col <- grep("Survival status", names(dt), value = TRUE)
  if (!length(os_col) || !length(st_col))
    stop(paste(cancer, ": OS columns not found"))
  dt[, os_days   := suppressWarnings(as.numeric(get(os_col[1])))]
  dt[, os_status := suppressWarnings(as.numeric(get(st_col[1])))]
  dt <- dt[!is.na(sample_id) & !is.na(os_days) & !is.na(os_status) & os_days > 0]
  dt[, .(sample_id, os_days, os_status)]
}


load_score_matrix <- function(cancer, method) {
  path <- file.path(SCORE_BASE, cancer, paste0(method, ".csv"))
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt$Complex
  mat
}


top_complexes_by_cv <- function(mat, pct) {
  means <- rowMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 1, sd, na.rm = TRUE)
  cv    <- ifelse(abs(means) > 1e-9, sds / abs(means), NA_real_)
  thresh <- quantile(cv, 1 - pct, na.rm = TRUE)
  names(cv)[!is.na(cv) & cv >= thresh]
}


cluster_logrank <- function(mat, clin) {
  common <- intersect(colnames(mat), clin$sample_id)
  if (length(common) < 10L) return(NULL)

  mat  <- mat[, common, drop = FALSE]
  clin <- clin[sample_id %in% common][match(common, sample_id)]

  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.nan(mat_scaled)] <- 0

  km <- kmeans(t(mat_scaled), centers = 2, nstart = N_STARTS, iter.max = 100)
  labels <- km$cluster

  mean1 <- mean(mat_scaled[, km$cluster == 1], na.rm = TRUE)
  mean2 <- mean(mat_scaled[, km$cluster == 2], na.rm = TRUE)
  if (mean1 > mean2) labels <- ifelse(labels == 1, 2, 1)

  df <- merge(
    data.table(sample_id = names(labels), cluster = labels),
    clin, by = "sample_id"
  )

  if (sum(df$cluster == 1) < 3L || sum(df$cluster == 2) < 3L) return(NULL)

  lr <- tryCatch(
    survdiff(Surv(os_days, os_status) ~ cluster, data = df),
    error = function(e) NULL
  )
  if (is.null(lr)) return(NULL)

  list(
    n_samples  = nrow(df),
    n_events   = sum(df$os_status),
    n_cluster1 = sum(df$cluster == 1),
    n_cluster2 = sum(df$cluster == 2),
    p_logrank  = pchisq(lr$chisq, df = 1, lower.tail = FALSE),
    cluster_df = df
  )
}


all_rows    <- list()
cluster_dfs <- list()

for (cancer in CANCERS) {
  clin <- tryCatch(load_clinical(cancer),
                   error = function(e) { message(e); NULL })
  if (is.null(clin)) next

  for (method in METHODS) {
    mat <- load_score_matrix(cancer, method)
    if (is.null(mat)) {
      cat("[", cancer, "] missing:", method, "\n"); next
    }

    for (cv_label in names(CV_CUTOFFS)) {
      pct    <- CV_CUTOFFS[cv_label]
      top_cx <- top_complexes_by_cv(mat, pct)
      mat_sub <- mat[top_cx, , drop = FALSE]

      cat("[", cancer, "][", method, "][", cv_label, "]",
          " n_cx=", nrow(mat_sub), " clustering ...\n")

      res <- cluster_logrank(mat_sub, clin)
      if (is.null(res)) {
        cat("[", cancer, "][", method, "][", cv_label, "] skipped\n"); next
      }

      cat("  p=", signif(res$p_logrank, 3), "\n")

      all_rows[[length(all_rows) + 1L]] <- data.table(
        cancer     = cancer,
        method     = method,
        cv_filter  = cv_label,
        n_cx_used  = nrow(mat_sub),
        n_samples  = res$n_samples,
        n_events   = res$n_events,
        n_cluster1 = res$n_cluster1,
        n_cluster2 = res$n_cluster2,
        p_logrank  = round(res$p_logrank, 6)
      )

      cdf <- res$cluster_df
      cdf[, cancer := cancer][, method := method][, cv_filter := cv_label]
      cluster_dfs[[length(cluster_dfs) + 1L]] <- cdf
    }
  }
}

result_dt <- rbindlist(all_rows, fill = TRUE)
fwrite(result_dt, file.path(OUT_DIR, "cluster_logrank.csv"))
cat("[saved] cluster_logrank.csv\n")

cluster_all <- rbindlist(cluster_dfs, fill = TRUE)
fwrite(cluster_all, file.path(OUT_DIR, "cluster_assignments.csv"))
cat("[saved] cluster_assignments.csv\n")

for (.cv in names(CV_CUTOFFS)) {
  cat("\n[cluster log-rank p-values |", .cv, "]\n")
  wide <- dcast(result_dt[cv_filter == .cv],
                cancer ~ method, value.var = "p_logrank")
  print(wide)
}

cat("[done]\n")
