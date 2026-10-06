rm(list = ls())

library(ggplot2)
library(ggrepel)
library(scales)

RESULT_DIR <- "/proj/c.zihao/work3/02TCGAnt/02stringTest/02analysis"
PLOT_DIR <- "/proj/c.zihao/work3/02TCGAnt/02stringTest/02analysis"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
metrics_dt <- read.csv(file.path(RESULT_DIR, "cluster_metrics.csv"))

METHOD_ORDER <- c(
  "netcomplex", "mean", "min", "z-score", "ssGSEA", "GSVA", "PLAGE", "netmean", "PerPAS"
)
METHOD_ORDER <- METHOD_ORDER[METHOD_ORDER %in% unique(metrics_dt$method)]
if (length(METHOD_ORDER) != 9L) {
  stop("Benchmark metrics are incomplete. Run 01_cal.R before plotting.")
}

FRAC_ORDER <- c("5%", "10%", "15%", "20%", "25%")
CANCER_ORDER <- c("BLCA", "BRCA", "CRC", "ESCA", "HNSC", "KICH", "KIRC",
                  "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD",
                  "UCEC")
method_palette <- c(
  "netcomplex" = "#000000", "mean" = "#D55E00", "min" = "#009E73",
  "z-score" = "#59a14f", "ssGSEA" = "#9c755f", "GSVA" = "#76b7b2",
  "PLAGE" = "#edc948", "netmean" = "#af7aa1", "PerPAS" = "#6b7280"
)
method_label <- c(
  "netcomplex" = "NetComplex", "mean" = "Mean", "min" = "Min",
  "z-score" = "ZScore", "ssGSEA" = "ssGSEA", "GSVA" = "GSVA",
  "PLAGE" = "PLAGE", "netmean" = "NetMean", "PerPAS" = "PerPAS"
)

metrics_dt <- metrics_dt[metrics_dt$method %in% METHOD_ORDER, ]
metrics_dt$method <- factor(metrics_dt$method, levels = METHOD_ORDER)
metrics_dt$feature_fraction <- factor(metrics_dt$feature_fraction, levels = FRAC_ORDER)
metrics_dt$cancer <- factor(metrics_dt$cancer, levels = CANCER_ORDER)

format_method <- function(x) {
  factor(as.character(x), levels = METHOD_ORDER, labels = unname(method_label[METHOD_ORDER]))
}

theme_pub <- function() {
  theme_classic(base_size = 10, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.35, color = "black"),
          axis.ticks = element_line(linewidth = 0.35, color = "black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major.y = element_line(linewidth = 0.2, color = "grey88")
          )
}

metric_map <- c(auc = "AUC", cluster_agreement = "Clustering accuracy")
line_dt <- do.call(rbind, lapply(names(metric_map), function(metric) {
  d <- aggregate(metrics_dt[[metric]],
                 by = list(method = metrics_dt$method,
                           feature_fraction = metrics_dt$feature_fraction),
                 FUN = mean, na.rm = TRUE)
  names(d)[3] <- "mean_metric"
  d$metric <- metric_map[[metric]]
  d
}))
line_dt$method <- factor(line_dt$method, levels = METHOD_ORDER)
line_dt$feature_fraction <- factor(line_dt$feature_fraction, levels = FRAC_ORDER)
line_dt$method_label <- format_method(line_dt$method)
line_dt$metric <- factor(line_dt$metric, levels = unname(metric_map))
end_dt <- line_dt[line_dt$feature_fraction == "25%", ]

make_curve <- function(metric_name, y_label, output_file) {
  d <- line_dt[line_dt$metric == metric_name, ]
  d_end <- end_dt[end_dt$metric == metric_name, ]
  p <- ggplot(d, aes(feature_fraction, mean_metric, group = method_label, color = method_label)) +
    geom_line(aes(linewidth = method == "netcomplex"), show.legend = FALSE) +
    geom_point(aes(size = method == "netcomplex"), show.legend = FALSE) +
    geom_text_repel(data = d_end, aes(label = method_label), colour = "grey15",
                    segment.colour = unname(method_palette[as.character(d_end$method)]),
                    size = 2.7, hjust = 0, direction = "y", nudge_x = 0.18,
                    segment.size = 0.35, min.segment.length = 0,
                    box.padding = 0.12, max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = unname(method_palette[METHOD_ORDER])) +
    scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 1.1)) +
    scale_size_manual(values = c("FALSE" = 1.1, "TRUE" = 2.0)) +
    scale_y_continuous(limits = c(0, 1.02), expand = expansion(mult = c(0.01, 0.02))) +
    labs(title = "Benchmark: tumor-normal discrimination (14 cancers)",
         x = "features retained (top CV fraction)", y = y_label) +
    theme_pub() +
    theme(axis.text.x = element_text(face = "bold"), legend.position = "none") +
    scale_x_discrete(expand = expansion(add = c(0.08, 0.8)))
  ggsave(file.path(PLOT_DIR, output_file), p, width = 5, height = 5)
}

make_curve("AUC", "Mean AUC", "3B_curve_auc_bench.pdf")
make_curve("Clustering accuracy", "Mean clustering accuracy", "3A_curve_accuracy_bench.pdf")

heat_dt <- rbind(
  transform(aggregate(auc ~ cancer + method, data = metrics_dt, FUN = mean, na.rm = TRUE),
            metric = "AUC", value = auc)[, c("cancer", "method", "metric", "value")],
  transform(aggregate(cluster_agreement ~ cancer + method, data = metrics_dt, FUN = mean, na.rm = TRUE),
            metric = "Clustering accuracy", value = cluster_agreement)[, c("cancer", "method", "metric", "value")]
)
heat_average_dt <- aggregate(value ~ method + metric, data = heat_dt, FUN = mean, na.rm = TRUE)
heat_average_dt$cancer <- "Average"
heat_dt <- rbind(heat_dt, heat_average_dt[, c("cancer", "method", "metric", "value")])
heat_dt$method <- factor(heat_dt$method, levels = METHOD_ORDER)
heat_dt$method_label <- format_method(heat_dt$method)
heat_dt$cancer <- factor(heat_dt$cancer, levels = c(CANCER_ORDER, "Average"))
heat_dt$metric <- factor(heat_dt$metric, levels = c("AUC", "Clustering accuracy"))

p_heat <- ggplot(heat_dt, aes(cancer, method_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.4,
            fontface = "bold", color = "#1F2937") +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#DCEFEA", "#8ACFC0", "#4C91BD", "#356FA5"),
    limits = c(0.4, 1), oob = squish, name = "Score"
  ) +
  facet_wrap(~ metric, nrow = 1) +
  labs(title = "Cancer-level tumor-normal discrimination across CV thresholds",
       subtitle = "Values are averaged across the 5%, 10%, 15%, 20%, and 25% feature fractions",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, face = "bold", size = 8),
    axis.text.y = element_text(face = "bold"),
    strip.background = element_rect(fill = "#E8EEF5", color = "#A8B7C8", linewidth = 0.4),
    strip.text = element_text(face = "bold", color = "#23374D"),
    legend.position = "right"
  )
ggsave(file.path(PLOT_DIR, "SupFigure1_heat_bench.pdf"), p_heat, width = 13, height = 5)

ranking_dt <- aggregate(cbind(cluster_agreement, auc) ~ method,
                        data = metrics_dt[metrics_dt$feature_fraction == "5%", ],
                        FUN = function(x) round(mean(x, na.rm = TRUE), 4))
ranking_dt <- ranking_dt[order(-ranking_dt$auc, -ranking_dt$cluster_agreement), ]
cat("\n=== NetComplex vs benchmark ranking at 5% ===\n")
print(ranking_dt)
