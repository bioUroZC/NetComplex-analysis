
RES_DIR <- "/proj/c.zihao/work3/01protein/01stringTest/02analysis"
IN_FILE <- file.path(RES_DIR, "concord.csv")
OUT_FILE <- file.path(RES_DIR, "netcomplex_vs_ablations.csv")

MAIN_PC <- 3L
FULL_METHOD <- "netcomplex"
ABLATIONS <- c("noRWR", "noRank", "degRandPPI", "permNodeValues")
METRICS <- c(mantel_r = "Mantel r",
             procrustes_r = "Symmetric Procrustes r")

if (!file.exists(IN_FILE)) stop("Missing concordance table: ", IN_FILE)
res <- read.csv(IN_FILE, stringsAsFactors = FALSE)

required <- c("comparison_group", "cancer", "method", "n_pc", names(METRICS))
missing <- setdiff(required, names(res))
if (length(missing)) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

res <- res[
  res$comparison_group == "ablation" &
    res$n_pc == MAIN_PC &
    res$method %in% c(FULL_METHOD, ABLATIONS),
  , drop = FALSE
]

compare_metric <- function(metric) {
  cohort_mean <- aggregate(res[[metric]],
                           by = list(cancer = res$cancer, method = res$method),
                           FUN = mean, na.rm = TRUE)
  names(cohort_mean)[3] <- "value"

  full <- cohort_mean[cohort_mean$method == FULL_METHOD, c("cancer", "value")]
  names(full)[2] <- "netcomplex_value"

  rows <- lapply(ABLATIONS, function(ablation) {
    ablated <- cohort_mean[cohort_mean$method == ablation, c("cancer", "value")]
    names(ablated)[2] <- "ablation_value"
    paired <- merge(full, ablated, by = "cancer", all = FALSE)
    if (nrow(paired) < 3L) stop("Too few shared cohorts for ", metric, " vs ", ablation)

    test <- wilcox.test(paired$netcomplex_value, paired$ablation_value,
                        paired = TRUE, alternative = "greater", exact = TRUE)
    difference <- paired$netcomplex_value - paired$ablation_value
    data.frame(
      metric = unname(METRICS[[metric]]),
      ablation = ablation,
      n_cancers = nrow(paired),
      mean_netcomplex = mean(paired$netcomplex_value),
      mean_ablation = mean(paired$ablation_value),
      mean_difference = mean(difference),
      median_difference = median(difference),
      cancers_netcomplex_higher = sum(difference > 0),
      p_value = unname(test$p.value),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$p_value_BH <- p.adjust(out$p_value, method = "BH")
  out
}

results <- do.call(rbind, lapply(names(METRICS), compare_metric))
results$ablation <- factor(results$ablation, levels = ABLATIONS)
results <- results[order(results$metric, results$ablation), ]
results$ablation <- as.character(results$ablation)

write.csv(results, OUT_FILE, row.names = FALSE)
print(results, row.names = FALSE, digits = 5)
cat("Wrote:", OUT_FILE, "\n")
