rm(list = ls())

library(data.table)
library(ggplot2)

OUT_DIR <- "/proj/c.zihao/work3/03TCGAt/03stringTest/02analysis"
result_dt <- fread(file.path(OUT_DIR, "pancancer_results.csv"))
silhouette_by_cancer_dt <- fread(file.path(OUT_DIR, "pancancer_silhouette_by_cancer.csv"))

METHOD_ORDER <- c("NetComplex", "mean", "min", "z-score",
                  "ssGSEA", "GSVA", "PLAGE", "netmean", "PerPAS")
METHOD_ORDER <- METHOD_ORDER[METHOD_ORDER %in% unique(result_dt$method)]
FRAC_ORDER <- c("5%", "10%", "15%", "20%", "25%")
FRAC_ORDER <- FRAC_ORDER[FRAC_ORDER %in% unique(result_dt$feature_fraction)]

method_label <- c(
  "NetComplex" = "NetComplex", "mean" = "Mean", "min" = "Min",
  "z-score" = "ZScore", "ssGSEA" = "ssGSEA", "GSVA" = "GSVA",
  "PLAGE" = "PLAGE", "netmean" = "NetMean", "PerPAS" = "PerPAS"
)
method_palette <- c(
  "NetComplex" = "#000000", "mean" = "#D55E00", "min" = "#009E73",
  "z-score" = "#CC79A7", "ssGSEA" = "#0072B2", "GSVA" = "#56B4E9",
  "PLAGE" = "#E69F00", "netmean" = "#6A3D9A", "PerPAS" = "#666666"
)

plot_dt <- melt(
  result_dt,
  id.vars = c("method", "feature_fraction"),
  measure.vars = c("ARI", "silhouette"),
  variable.name = "metric", value.name = "score"
)
plot_dt[, method := factor(method, levels = METHOD_ORDER)]
plot_dt[, method_label := factor(as.character(method), levels = METHOD_ORDER,
                                 labels = unname(method_label[METHOD_ORDER]))]
plot_dt[, feature_fraction := factor(feature_fraction, levels = FRAC_ORDER)]
plot_dt[, metric := factor(metric, levels = c("ARI", "silhouette"),
                           labels = c("Adjusted Rand Index", "Mean silhouette"))]

theme_pub <- function() {
  theme_classic(base_size = 10, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "#E8EEF5", color = "#A8B7C8", linewidth = 0.4),
      strip.text = element_text(face = "bold", color = "#23374D"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey88"),
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

make_curve <- function(metric_name, y_label, output_file, add_zero_line = FALSE) {
  d <- plot_dt[metric == metric_name]
  p <- ggplot(d, aes(feature_fraction, score, group = method_label, color = method_label)) +
    geom_line(aes(linewidth = method == "NetComplex"),
              show.legend = c(colour = TRUE, linewidth = FALSE)) +
    geom_point(aes(size = method == "NetComplex"), show.legend = FALSE) +
    scale_color_manual(values = unname(method_palette[METHOD_ORDER]), name = NULL) +
    scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 1.1)) +
    scale_size_manual(values = c("FALSE" = 1.1, "TRUE" = 2.0)) +
    guides(linewidth = "none", size = "none") +
    labs(
      title = "Pan-cancer structure preservation across feature fractions",
      x = "features retained (top CV fraction)", y = y_label
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(face = "bold"), legend.text = element_text(size = 8.5))

  if (add_zero_line) {
    p <- p + geom_hline(yintercept = 0, color = "grey60", linewidth = 0.35, linetype = "dashed")
  }
  ggsave(file.path(OUT_DIR, output_file), p, width = 5, height = 5)
}

make_curve("Adjusted Rand Index", "Adjusted Rand Index", "curve_ari_bench.pdf")
make_curve("Mean silhouette", "Mean silhouette", "3C_curve_silhouette_bench.pdf", add_zero_line = TRUE)

CANCER_ORDER <- c(
  "BLCA", "BRCA", "CRC", "GBM", "KIRC", "LGG", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PRAD", "SKCM", "UVM"
)
silhouette_by_cancer_dt[, method := factor(method, levels = METHOD_ORDER)]
silhouette_by_cancer_dt[, method_label := factor(as.character(method), levels = METHOD_ORDER,
                                                  labels = unname(method_label[METHOD_ORDER]))]
silhouette_by_cancer_dt[, cancer := factor(cancer, levels = CANCER_ORDER)]
silhouette_by_cancer_dt[, feature_fraction := factor(feature_fraction, levels = FRAC_ORDER)]

heat_dt <- silhouette_by_cancer_dt[, .(silhouette = mean(silhouette, na.rm = TRUE)),
                                   by = .(method, cancer)]
heat_dt <- heat_dt[cancer %in% CANCER_ORDER]
heat_average_dt <- heat_dt[, .(silhouette = mean(silhouette, na.rm = TRUE)), by = method]
heat_average_dt[, cancer := "Average"]
heat_dt <- rbind(heat_dt, heat_average_dt, use.names = TRUE)
heat_dt[, method := factor(method, levels = METHOD_ORDER)]
heat_dt[, method_label := factor(as.character(method), levels = METHOD_ORDER,
                                 labels = unname(method_label[METHOD_ORDER]))]
heat_dt[, cancer := factor(cancer, levels = c(CANCER_ORDER, "Average"))]

p_heat <- ggplot(heat_dt,
                 aes(x = cancer, y = method_label, fill = silhouette)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", silhouette)), color = "#1F2937", size = 1.7) +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#DCEFEA", "#8ACFC0", "#4C91BD", "#356FA5"),
    limits = range(heat_dt$silhouette, na.rm = TRUE),
    name = "Mean\nsilhouette"
  ) +
  labs(
    title = "Cancer-specific structure preservation across feature fractions",
    subtitle = "Fourteen TCGA tumour-only cancer types; mean silhouette averaged across five feature fractions",
    x = NULL, y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 7),
    axis.text.y = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(OUT_DIR, "SupFigure3_heat_bench.pdf"), p_heat, width = 9, height = 5)

ranking_dt <- result_dt[feature_fraction == "5%"][order(-ARI, -silhouette)]
cat("\n=== Global ranking at 5% (by ARI) ===\n")
print(ranking_dt[, .(method, ARI, silhouette)])
