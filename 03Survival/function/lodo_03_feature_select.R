
library(survival)

get_top_features <- function(mat, top_k) {
  means <- rowMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 1, sd, na.rm = TRUE)
  cv    <- sds / means
  cv[!is.finite(cv)] <- NA_real_
  ord <- order(cv, decreasing = TRUE, na.last = TRUE)
  rownames(mat)[ord][seq_len(min(top_k, nrow(mat)))]
}

top_correlated_features <- function(x, y, top_k) {
  cor_vec <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs")[, 1])
  cor_vec[!is.finite(cor_vec)] <- 0
  names(sort(abs(cor_vec), decreasing = TRUE))[seq_len(min(top_k, length(cor_vec)))]
}

cox_screen_topN <- function(train_df, feature_cols, topN) {
  pvals <- vapply(feature_cols, function(g) {
    fml <- as.formula(paste0("Surv(OS_Time, OS) ~ `", g, "`"))
    tryCatch({
      fit <- coxph(fml, data = train_df, ties = "efron")
      s <- summary(fit)
      if (nrow(s$coefficients) == 0) return(1.0)
      p <- s$coefficients[1, "Pr(>|z|)"]
      if (!is.finite(p)) 1.0 else p
    }, error = function(e) 1.0)
  }, numeric(1))
  names(pvals) <- feature_cols
  names(sort(pvals))[seq_len(min(topN, length(pvals)))]
}
