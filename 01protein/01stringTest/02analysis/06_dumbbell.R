
library(ggplot2)

RES_DIR <- "/proj/c.zihao/work3/01protein/01stringTest/02analysis"
RES_FILE <- file.path(RES_DIR, "concord.csv")
OUT_FILES <- c(
  mantel_r = file.path(RES_DIR, "2C_dumbbell_bench.pdf"),
  procrustes_r = file.path(RES_DIR, "2D_dumbbell_bench_procrustes.pdf")
)

MAIN_PC <- 3
METHODS <- c("netcomplex", "mean", "ssgsea")
CANCER_ORDER <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD",
                  "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHOD_LABELS <- c(netcomplex = "NetComplex", mean = "Mean", ssgsea = "ssGSEA")
METHOD_COLORS <- c(netcomplex = "#000000", mean = "#D55E00", ssgsea = "#0072B2")

if (!file.exists(RES_FILE)) {
  stop("Missing concordance table: ", RES_FILE)
}

res <- read.csv(RES_FILE, stringsAsFactors = FALSE)
required <- c("comparison_group", "cancer", "method", "n_pc",
              "cv_top_fraction", names(OUT_FILES))
missing <- setdiff(required, names(res))
if (length(missing) > 0) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

plot_data <- res[
  res$comparison_group == "benchmark" &
    res$n_pc == MAIN_PC &
    res$method %in% METHODS,
  c("cancer", "method", "cv_top_fraction", names(OUT_FILES))
]

make_dumbbell <- function(metric, metric_label, out_file) {
  summary_data <- aggregate(
    plot_data[[metric]], by = list(cancer = plot_data$cancer, method = plot_data$method),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(summary_data)[3] <- "value"

  expected <- expand.grid(cancer = CANCER_ORDER, method = METHODS,
                          stringsAsFactors = FALSE)
  missing_pairs <- merge(expected, summary_data,
                         by = c("cancer", "method"), all.x = TRUE)
  if (anyNA(missing_pairs$value)) {
    absent <- missing_pairs[is.na(missing_pairs$value), c("cancer", "method")]
    stop("Missing cancer-method result(s): ",
         paste(paste(absent$cancer, absent$method, sep = "/"), collapse = ", "))
  }

  summary_data$cancer <- factor(summary_data$cancer, levels = rev(CANCER_ORDER))
  summary_data$method <- factor(summary_data$method, levels = METHODS)
  ranges <- aggregate(value ~ cancer, data = summary_data,
                      FUN = function(x) c(min = min(x), max = max(x)))
  ranges <- data.frame(cancer = ranges$cancer,
                       xmin = ranges$value[, "min"], xmax = ranges$value[, "max"])

  p <- ggplot() +
    geom_segment(data = ranges,
                 aes(x = xmin, xend = xmax, y = cancer, yend = cancer),
                 colour = "grey70", linewidth = 0.7) +
    geom_point(data = summary_data,
               aes(x = value, y = cancer, colour = method, shape = method),
               size = 3.0) +
    scale_colour_manual(values = METHOD_COLORS, labels = METHOD_LABELS, name = NULL) +
    scale_shape_manual(values = c(netcomplex = 16, mean = 17, ssgsea = 15),
                       labels = METHOD_LABELS, name = NULL) +
    labs(title = "Cross-omics concordance by cancer type",
         subtitle = paste("Mean", metric_label, "across top-CV fractions"),
         x = paste("Mean", metric_label), y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal")

  ggsave(out_file, p, width = 7.0, height = 4.6)
  cat("Wrote:", out_file, "\n")
}

make_dumbbell("mantel_r", "Mantel r", OUT_FILES[["mantel_r"]])
make_dumbbell("procrustes_r", "symmetric Procrustes r", OUT_FILES[["procrustes_r"]])
