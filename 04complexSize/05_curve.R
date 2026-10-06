
library(ggplot2)
library(ggrepel)

BR       <- "/proj/c.zihao/work3/04complexSize"
OUT_DIR  <- BR
RES_FILE <- file.path(OUT_DIR, "concord.csv")
OUT_BENCH <- file.path(OUT_DIR, "Figure4F_curve_bench.pdf")
OUT_PROCRUSTES <- file.path(OUT_DIR, "Figure4G_curve_bench_procrustes.pdf")
MAIN_PC  <- 3

SIZE_BIN_LEVELS <- c("1", "2", "3", "4", "5-7", "8-10", "11+")
METHODS_BY_GROUP <- list(
  benchmark = c("netcomplex", "mean", "min", "zscore", "netmean", "perpas",
                "gsva", "ssgsea", "plage")
)
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

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey88"),
      legend.position = "none",
      plot.margin = margin(6, 8, 6, 8)
    )
}

if (!file.exists(RES_FILE)) stop("no concordance table; run 03_calc.R first: ", RES_FILE)
res <- read.csv(RES_FILE, stringsAsFactors = FALSE)
res <- res[res$n_pc == MAIN_PC, , drop = FALSE]

make_aggregate <- function(metric) {
  agg <- aggregate(res[[metric]],
                   by = list(comparison_group = res$comparison_group,
                             method = res$method,
                             size_bin = res$size_bin),
                   FUN = mean, na.rm = TRUE)
  names(agg)[4] <- "value"
  agg$label <- METHOD_LABELS[agg$method]
  agg$size_bin <- factor(agg$size_bin, levels = SIZE_BIN_LEVELS)
  agg$x <- as.integer(agg$size_bin)
  agg
}

panel <- function(agg, group, title, y_label) {
  d <- agg[agg$comparison_group == group &
           agg$method %in% METHODS_BY_GROUP[[group]], , drop = FALSE]
  ends <- d[d$x == max(d$x), , drop = FALSE]

  ggplot(d, aes(x, value, colour = method, group = method)) +
    geom_line(aes(linewidth = method == "netcomplex"), show.legend = FALSE) +
    geom_point(aes(size = method == "netcomplex"), show.legend = FALSE) +
    geom_text_repel(data = ends, aes(label = label),
                    colour = "grey15", segment.colour = PALETTE[ends$method],
                    size = 2.7, hjust = 0,
                    direction = "y", nudge_x = 0.12, segment.size = 0.4,
                    min.segment.length = 0,
                    box.padding = 0.12, max.overlaps = Inf, show.legend = FALSE) +
    scale_colour_manual(values = PALETTE) +
    scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 1.1)) +
    scale_size_manual(values = c("FALSE" = 1.1, "TRUE" = 2.0)) +
    scale_x_continuous(breaks = seq_along(SIZE_BIN_LEVELS),
                       labels = SIZE_BIN_LEVELS,
                       expand = expansion(mult = c(0.04, 0.26))) +
    labs(x = "complex size (members)",
         y = y_label,
         title = title) +
    theme_pub()
}

p1 <- panel(make_aggregate("mantel_r"), "benchmark",
            "Benchmark: concordance by complex size (10 cancers, PC3)",
            "Mean Mantel r")
p2 <- panel(make_aggregate("procrustes_r"), "benchmark",
            "Benchmark: concordance by complex size (10 cancers, PC3)",
            "Mean symmetric Procrustes r")

ggsave(OUT_BENCH, p1, width = 6.4, height = 4.4)
ggsave(OUT_PROCRUSTES, p2, width = 6.4, height = 4.4)

cat("Plot ->", OUT_BENCH, "\n")
cat("Plot ->", OUT_PROCRUSTES, "\n")
