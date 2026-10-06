
library(ggplot2)

BR        <- "/proj/c.zihao/work3/04ComplexPortal"
OUT_DIR   <- file.path(BR, "02analysis")
RES_FILE  <- file.path(OUT_DIR, "concord.csv")
OUT_BENCH <- file.path(OUT_DIR, "Figure4C_boxplot_bench.pdf")
OUT_PROCRUSTES <- file.path(OUT_DIR, "Figure4D_boxplot_bench_procrustes.pdf")
MAIN_PC   <- 3

METHODS <- c("netcomplex", "mean", "min", "zscore", "netmean", "perpas",
             "gsva", "ssgsea", "plage")
METHOD_LABELS <- c(
  netcomplex = "NetComplex", mean = "Mean", min = "Min", zscore = "Z-score",
  netmean = "NetMean", perpas = "PerPAS", gsva = "GSVA", ssgsea = "ssGSEA",
  plage = "PLAGE"
)
OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#6A3D9A", "#666666")
PALETTE <- c(
  netcomplex = "#000000",
  setNames(OKABE, c("ssgsea", "mean", "min", "zscore", "netmean", "perpas", "gsva", "plage"))
)
FRAC_ORDER <- c(0.05, 0.10, 0.15, 0.20, 0.25)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey88"),
      strip.background = element_rect(fill = "#E8EEF5", color = "#A8B7C8", linewidth = 0.4),
      strip.text = element_text(face = "bold", color = "#23374D"),
      legend.position = "none"
    )
}

if (!file.exists(RES_FILE)) stop("no concordance table; run 03_calc.R first: ", RES_FILE)
res <- read.csv(RES_FILE, stringsAsFactors = FALSE)

plot_boxplot <- function(metric, metric_label, output_file) {
  d <- res[
    res$n_pc == MAIN_PC &
      res$comparison_group == "benchmark" &
      res$method %in% METHODS,
    c("cancer", "method", "cv_top_fraction", metric)
  ]
  names(d)[4] <- "value"
  if (length(unique(d$cancer)) != 10L) stop("Expected results from 10 cancer cohorts.")
  d$method <- factor(d$method, levels = METHODS,
                     labels = unname(METHOD_LABELS[METHODS]))

  p <- ggplot(d, aes(x = method, y = value, fill = method, colour = method)) +
    geom_boxplot(width = 0.68, alpha = 0.45, outlier.shape = NA,
                 linewidth = 0.45) +
    geom_jitter(width = 0.12, height = 0, size = 1.35, alpha = 0.82) +
    scale_fill_manual(values = unname(PALETTE[METHODS])) +
    scale_colour_manual(values = unname(PALETTE[METHODS])) +
    labs(
      x = NULL,
      y = metric_label
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))

  ggsave(output_file, p, width = 6, height = 6)
  cat("Plot ->", output_file, "\n")
}

plot_boxplot("mantel_r", "Mantel r", OUT_BENCH)
plot_boxplot("procrustes_r", "Symmetric Procrustes r", OUT_PROCRUSTES)
