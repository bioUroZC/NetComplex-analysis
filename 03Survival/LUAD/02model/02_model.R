
library(data.table)
library(survival)

function_dir <- "/proj/c.zihao/work3/03Survival/function"
source(file.path(function_dir, "lodo_03_feature_select.R"))
source(file.path(function_dir, "lodo_04_model.R"))
source(file.path(function_dir, "lodo_05_eval.R"))

cancer    <- "LUAD"
model_dir <- "/proj/c.zihao/work3/03Survival/LUAD/02model"
data_dir  <- file.path(model_dir, "data")

feature_sets <- c("raw_expr", "netcomplex", "ssgsea", "mean")
n_stage1 <- 10000
n_stage2 <- 2000
n_stage3 <- 200
seed     <- 1


build_model_frame <- function(mat, samples, clin_sub) {
  df <- as.data.frame(t(mat[, samples, drop = FALSE]))
  df$Sample <- rownames(df)
  df <- merge(df, clin_sub[, c("Sample", "OS", "OS_Time")], by = "Sample")
  rownames(df) <- df$Sample
  df$Sample <- NULL
  df
}


clin <- as.data.frame(fread(file.path(data_dir, "OS.csv")))
datasets <- unique(clin$dataset)
cat(cancer, ": clinical n =", nrow(clin), ", events =", sum(clin$OS == 1), "\n")

results <- list()

for (fs in feature_sets) {
  path <- file.path(data_dir, paste0(fs, ".csv"))
  if (!file.exists(path)) stop("file not found: ", path)
  dt <- fread(path)
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt[[1]]
  mat <- mat[, colnames(mat) %in% clin$Sample, drop = FALSE]
  cat("\n", fs, ":", nrow(mat), "features x", ncol(mat), "samples\n")

  for (ds in datasets) {
    train_clin <- clin[clin$dataset != ds, ]
    test_clin  <- clin[clin$dataset == ds, ]
    train_samples <- colnames(mat)[colnames(mat) %in% train_clin$Sample]
    test_samples  <- colnames(mat)[colnames(mat) %in% test_clin$Sample]
    cat("  held out:", ds, "| train =", length(train_samples),
        "| test =", length(test_samples), "\n")

    if (length(train_samples) < 20 || length(test_samples) < 5) {
      cat("  skipped: too few samples\n")
      results[[length(results) + 1]] <- data.frame(
        cancer = cancer, feature_set = fs, dataset = ds,
        n_train = length(train_samples), n_test = length(test_samples),
        n_events_test = NA, n_features_used = NA, cindex = NA, mean_tAUC = NA
      )
      next
    }

    top1 <- get_top_features(mat[, train_samples, drop = FALSE], n_stage1)
    mat1 <- mat[top1, , drop = FALSE]

    train_df  <- build_model_frame(mat1, train_samples, train_clin)
    feat_cols <- setdiff(colnames(train_df), c("OS", "OS_Time"))

    top2 <- top_correlated_features(as.matrix(train_df[, feat_cols, drop = FALSE]),
                                    train_df$OS_Time, n_stage2)
    top3 <- cox_screen_topN(train_df, top2, min(n_stage3, length(top2)))

    x_train <- as.matrix(train_df[, top3, drop = FALSE])
    train_var <- apply(x_train, 2, var, na.rm = TRUE)
    keep_feat <- names(train_var)[is.finite(train_var) & train_var > 0]
    x_train <- x_train[, keep_feat, drop = FALSE]
    y_train <- Surv(train_df$OS_Time, train_df$OS)

    test_df <- build_model_frame(mat1[keep_feat, , drop = FALSE], test_samples, test_clin)
    x_test  <- as.matrix(test_df[, keep_feat, drop = FALSE])

    risk <- fit_lasso_cox_risk(x_train, y_train, x_test, train_df, test_df, keep_feat, seed)

    ci  <- safe_concordance_index(risk, test_df$OS_Time, test_df$OS)
    auc <- safe_mean_auc(risk, test_df$OS_Time, test_df$OS)
    cat("  features =", length(keep_feat), "| C-index =", round(ci, 4),
        "| mean tAUC =", round(auc, 4), "\n")

    results[[length(results) + 1]] <- data.frame(
      cancer = cancer, feature_set = fs, dataset = ds,
      n_train = nrow(train_df), n_test = nrow(test_df),
      n_events_test = sum(test_df$OS == 1), n_features_used = length(keep_feat),
      cindex = round(ci, 4), mean_tAUC = round(auc, 4)
    )
  }
}

results <- do.call(rbind, results)
fwrite(results, file.path(model_dir, "lodo_results.csv"))

cat("\nMedian C-index per feature set:\n")
print(aggregate(cindex ~ feature_set, data = results, FUN = median))
cat("done\n")
