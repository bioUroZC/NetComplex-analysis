
library(ggplot2)

RES_DIR <- "/proj/c.zihao/work3/01protein/01stringTest/03SD"
OUT_DIR <- "/proj/c.zihao/work3/01protein/01stringTest/03SD"
MAIN_PC <- 3
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

res <- read.csv(file.path(RES_DIR, "concord.csv"),
                stringsAsFactors = FALSE)
required_cols <- c("comparison_group", "cancer", "method", "n_pc", "mantel_r", "procrustes_r")
missing_cols <- setdiff(required_cols, names(res))
if (length(missing_cols) > 0) {
  stop("Results file does not use the grouped-comparison format; rerun 03_calc.R. Missing: ",
       paste(missing_cols, collapse = ", "))
}

CANCER_ORDER <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHODS_BY_GROUP <- list(
  benchmark = c("netcomplex", "mean", "min", "zscore", "netmean", "perpas", "gsva", "ssgsea", "plage"),
  ablation = c("netcomplex", "noRWR", "noRank", "degRandPPI", "permNodeValues")
)
METHOD_LABELS <- c(
  netcomplex = "NetComplex", mean = "Mean", min = "Min", zscore = "Z-score",
  netmean = "NetMean", perpas = "PerPAS", gsva = "GSVA", ssgsea = "ssGSEA",
  plage = "PLAGE", noRWR = "NoRWR", noRank = "NoRank",
  degRandPPI = "RewiredPPI", permNodeValues = "PermutedRanks"
)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}

plot_group_heatmap <- function(group) {
  expected_methods <- METHODS_BY_GROUP[[group]]
  main <- res[res$comparison_group == group & res$n_pc == MAIN_PC &
              res$method %in% expected_methods, , drop = FALSE]
  missing_methods <- setdiff(expected_methods, unique(main$method))
  if (length(missing_methods) > 0) {
    stop(sprintf("%s results missing method(s): %s", group,
                 paste(missing_methods, collapse = ", ")))
  }

  mantel <- aggregate(mantel_r ~ cancer + method, data = main, FUN = mean, na.rm = TRUE)
  mantel$metric <- "Mantel r"
  names(mantel)[3] <- "value"
  procrustes <- aggregate(procrustes_r ~ cancer + method, data = main, FUN = mean, na.rm = TRUE)
  procrustes$metric <- "Symmetric Procrustes r"
  names(procrustes)[3] <- "value"
  heat <- rbind(mantel, procrustes)

  avg <- aggregate(value ~ method + metric, data = heat, FUN = mean, na.rm = TRUE)
  avg$cancer <- "Average"
  heat <- rbind(heat, avg[, names(heat)])

  heat$cancer <- factor(heat$cancer, levels = c(CANCER_ORDER, "Average"))
  heat$method_label <- factor(METHOD_LABELS[heat$method],
                              levels = rev(METHOD_LABELS[expected_methods]))

  title <- sprintf("%s comparison: mean concordance across SD thresholds",
                   tools::toTitleCase(group))
  p <- ggplot(heat, aes(x = cancer, y = method_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.65) +
    geom_text(aes(label = sprintf("%.3f", value)), colour = "#1F2937", size = 2.8) +
    scale_fill_gradientn(
      colours = c("#F7F7F7", "#DCEFEA", "#8ACFC0", "#4C91BD", "#356FA5"),
      values = scales::rescale(c(-0.05, 0.15, 0.40, 0.65, 0.90)),
      limits = c(-0.05, 0.90),
      breaks = c(0.0, 0.2, 0.4, 0.6, 0.8),
      name = "Mean r"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    facet_wrap(~ metric, nrow = 1) +
    theme_pub(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      panel.background = element_rect(fill = "#FAFBFC", colour = NA),
      strip.background = element_rect(fill = "#E8EEF5", colour = "#A8B7C8", linewidth = 0.4),
      strip.text = element_text(face = "bold", colour = "#23374D"),
      plot.margin = margin(6, 10, 6, 38)
    )

  filename <- sprintf("heat_%s.pdf", c(benchmark = "bench", ablation = "abl")[group])
  ggsave(file.path(OUT_DIR, filename), p,
         width = 14.0, height = ifelse(group == "benchmark", 5.8, 4.8), device = cairo_pdf)
}

for (group in names(METHODS_BY_GROUP)) {
  plot_group_heatmap(group)
}

cat("Plots ->", OUT_DIR, "\n")
