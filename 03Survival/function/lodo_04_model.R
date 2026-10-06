
library(survival)
library(glmnet)

fit_lasso_cox_risk <- function(x_train, y_train, x_test, train_df, test_df,
                                keep_feat, seed = 1) {
  risk <- rep(0, nrow(x_test))
  if (ncol(x_train) == 0L) return(risk)

  nfolds <- max(3, min(5, nrow(x_train)))
  fit <- tryCatch({
    set.seed(seed)
    foldid <- sample(rep(seq_len(nfolds), length.out = nrow(x_train)))
    cvfit <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, foldid = foldid)
    glmnet(x_train, y_train, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
  }, error = function(e) NULL)

  if (!is.null(fit)) {
    return(as.numeric(predict(fit, newx = x_test, type = "link")))
  }

  if (length(keep_feat) > 0L) {
    fb_fml <- as.formula(paste0("Surv(OS_Time, OS) ~ `", keep_feat[1], "`"))
    fb_fit <- tryCatch(coxph(fb_fml, data = train_df, ties = "efron"), error = function(e) NULL)
    if (!is.null(fb_fit)) {
      fb_pred <- tryCatch(predict(fb_fit, newdata = test_df, type = "lp"), error = function(e) NULL)
      if (!is.null(fb_pred)) return(as.numeric(fb_pred))
    }
  }
  risk
}
