
library(survcomp)
library(timeROC)

safe_concordance_index <- function(risk, time, event) {
  keep <- is.finite(risk) & is.finite(time) & is.finite(event)
  risk <- risk[keep]; time <- time[keep]; event <- event[keep]
  if (length(unique(risk)) <= 1) return(0.5)
  ci <- tryCatch(
    concordance.index(x = risk, surv.time = time, surv.event = event)$c.index,
    error = function(e) NA
  )
  if (is.na(ci)) 0.5 else ci
}

safe_mean_auc <- function(risk, time, event) {
  keep <- is.finite(risk) & is.finite(time) & is.finite(event)
  risk <- risk[keep]; time <- time[keep]; event <- event[keep]
  if (length(unique(risk)) <= 1) return(0.5)
  tp <- unique(as.numeric(quantile(time, c(0.05, 0.25, 0.5, 0.75, 0.95))))
  tp <- tp[is.finite(tp) & tp > 0]
  if (length(tp) == 0) return(0.5)
  roc <- tryCatch(
    timeROC(T = time, delta = event, marker = risk, cause = 1, times = tp),
    error = function(e) NULL
  )
  if (is.null(roc)) return(0.5)
  mean(roc$AUC, na.rm = TRUE)
}
