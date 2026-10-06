
library(ggplot2)
library(ggrepel)
library(ggbreak)

BR       <- "/proj/c.zihao/work3/01protein/01stringTest"
OUT_DIR  <- file.path(BR, "03SD")
RES_FILE <- file.path(OUT_DIR, "concord.csv")
OUT_BENCH <- file.path(OUT_DIR, "SupplementaryFigure7.pdf")
OUT_ABL <- file.path(OUT_DIR, "curve_abl.pdf")
MAIN_PC  <- 3

METHODS_BY_GROUP <- list(
  benchmark = c("netcomplex", "mean", "min", "zscore", "netmean", "perpas",
                "gsva", "ssgsea", "plage"),
  ablation  = c("netcomplex", "noRWR", "noRank", "degRandPPI", "permNodeValues")
)
METHOD_LABELS <- c(
  netcomplex = "NetComplex", mean = "Mean", min = "Min", zscore = "Z-score",
  netmean = "NetMean", perpas = "PerPAS", gsva = "GSVA", ssgsea = "ssGSEA",
  plage = "PLAGE", noRWR = "NoRWR", noRank = "NoRank",
  degRandPPI = "RewiredPPI", permNodeValues = "PermutedRanks"
)

OKABE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#6A3D9A", "#666666")
PALETTE <- c(
  netcomplex = "#000000",
  setNames(OKABE, c("ssgsea", "mean", "min", "zscore", "netmean", "perpas", "gsva", "plage")),
  setNames(OKABE[1:4], c("noRWR", "noRank", "degRandPPI", "permNodeValues"))
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

agg <- aggregate(mantel_r ~ comparison_group + method + cv_top_fraction,
                 data = res, FUN = mean, na.rm = TRUE)
agg$label <- METHOD_LABELS[agg$method]

panel <- function(group, title) {
  d <- agg[agg$comparison_group == group &
           agg$method %in% METHODS_BY_GROUP[[group]], , drop = FALSE]
  ends <- d[d$cv_top_fraction == max(d$cv_top_fraction), , drop = FALSE]

  if (group == "benchmark") {
    return(
      ggplot(d, aes(cv_top_fraction, mantel_r, colour = method, group = method)) +
        geom_line(aes(linewidth = method == "netcomplex"), show.legend = FALSE) +
        geom_point(aes(size = method == "netcomplex"), show.legend = FALSE) +
        geom_text_repel(data = ends, aes(label = label),
                        colour = "grey15", segment.colour = PALETTE[ends$method],
                        size = 2.7, hjust = 0,
                        direction = "y", nudge_x = 0.006, segment.size = 0.4,
                        min.segment.length = 0,
                        box.padding = 0.12, max.overlaps = Inf, show.legend = FALSE) +
        scale_colour_manual(values = PALETTE) +
        scale_linewidth_manual(values = c("FALSE" = 0.5, "TRUE" = 1.1)) +
        scale_size_manual(values = c("FALSE" = 1.1, "TRUE" = 2.0)) +
        scale_x_continuous(breaks = sort(unique(d$cv_top_fraction)),
                           labels = function(v) paste0(100 * v, "%"),
                           expand = expansion(mult = c(0.04, 0.26))) +
        labs(x = "features retained (top SD fraction)",
             y = "Mean Mantel r",
             title = title) +
        theme_pub()
    )
  }

  p <- ggplot(d, aes(cv_top_fraction, mantel_r, colour = method, group = method)) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = PALETTE, labels = METHOD_LABELS, name = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(breaks = sort(unique(d$cv_top_fraction)),
                       labels = function(v) paste0(100 * v, "%"),
                       expand = expansion(mult = c(0.04, 0.04))) +
    labs(x = NULL, y = "Mean Mantel r", title = title) +
    theme_pub() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 8.5), axis.title.x = element_blank())

  p + scale_y_break(c(0.05, 0.42), scales = 0.6, space = 0.06)
}

p1 <- panel("benchmark", "Benchmark: concordance with matched proteome (10 cancers, PC3)")
p2 <- panel("ablation",  "Ablation: concordance with matched proteome (10 cancers, PC3)")

ggsave(OUT_BENCH, p1, width = 6.4, height = 4.4)
ggsave(OUT_ABL, p2, width = 6.8, height = 5.0)

cat("Plot ->", OUT_BENCH, "\n")
cat("Plot ->", OUT_ABL, "\n")
