rm(list = ls())

library(ggplot2)

RESULT_DIR <- "/proj/c.zihao/work3/02TCGAnt/02stringTest/02analysis"
PLOT_DIR <- "/proj/c.zihao/work3/02TCGAnt/02stringTest/02analysis"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

metrics_dt <- read.csv(file.path(RESULT_DIR, "cluster_metrics.csv"))
theme_pub <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#243B53", fill = NA, linewidth = 0.7),
      axis.text.y = element_text(face = "bold"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

plot_dumbbell <- function(metric, metric_label) {
  top5_dt <- metrics_dt[
    metrics_dt$feature_fraction == "5%" &
      metrics_dt$method %in% c("netcomplex", "mean", "ssGSEA"),
    c("cancer", "method", metric)
  ]
  names(top5_dt)[3] <- "value"

  if (length(unique(top5_dt$cancer)) != 14L || length(unique(top5_dt$method)) != 3L) {
    stop("The 5% results for NetComplex, Mean, and ssGSEA are incomplete. Run 01_cal.R first.")
  }

  wide_dt <- reshape(top5_dt, idvar = "cancer", timevar = "method", direction = "wide")
  names(wide_dt) <- sub("value\\.", "", names(wide_dt))
  names(wide_dt)[names(wide_dt) == "netcomplex"] <- "NetComplex"
  names(wide_dt)[names(wide_dt) == "mean"] <- "Mean"
  names(wide_dt)[names(wide_dt) == "ssGSEA"] <- "ssGSEA"
  wide_dt$difference <- wide_dt$NetComplex - pmax(wide_dt$Mean, wide_dt$ssGSEA)
  wide_dt$cancer <- factor(wide_dt$cancer, levels = wide_dt$cancer[order(wide_dt$difference)])

  point_dt <- rbind(
    data.frame(cancer = wide_dt$cancer, method = "NetComplex", value = wide_dt$NetComplex),
    data.frame(cancer = wide_dt$cancer, method = "Mean", value = wide_dt$Mean),
    data.frame(cancer = wide_dt$cancer, method = "ssGSEA", value = wide_dt$ssGSEA)
  )
  point_dt$method <- factor(point_dt$method, levels = c("NetComplex", "Mean", "ssGSEA"))
  range_dt <- aggregate(value ~ cancer, data = point_dt,
                        FUN = function(x) c(min = min(x), max = max(x)))
  range_dt <- data.frame(
    cancer = range_dt$cancer,
    xmin = range_dt$value[, "min"],
    xmax = range_dt$value[, "max"]
  )

  p <- ggplot(wide_dt, aes(y = cancer)) +
    geom_segment(data = range_dt, aes(x = xmin, xend = xmax, y = cancer, yend = cancer),
                 color = "grey70", linewidth = 0.7) +
    geom_point(data = point_dt, aes(x = value, color = method, shape = method), size = 2.8) +
    scale_color_manual(values = c("NetComplex" = "#000000", "Mean" = "#D55E00", "ssGSEA" = "#0072B2"), name = NULL) +
    scale_shape_manual(values = c("NetComplex" = 16, "Mean" = 17, "ssGSEA" = 15), name = NULL) +
    scale_x_continuous(limits = c(0.4, 1.02), breaks = seq(0.4, 1.0, 0.1)) +
    labs(
      title = "Tumor-normal discrimination by cancer type",
      subtitle = paste0(metric_label, " at the top 5% CV-ranked complexes"),
      x = metric_label, y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "bottom", legend.direction = "horizontal")

  p
}

p_accuracy <- plot_dumbbell("cluster_agreement", "Accuracy")
p_auc <- plot_dumbbell("auc", "AUC")

grDevices::pdf(file.path(PLOT_DIR, "SupFigure2_dumbbell_bench.pdf"), width = 12.0, height = 5)
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
print(p_accuracy, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_auc, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
grDevices::dev.off()
