library(data.table)
library(ggplot2)

IN_DIR  <- "/proj/c.zihao/work3/07multi/02survival"
OUT_DIR <- "/proj/c.zihao/work3/07multi/02survival"

dt <- fread(file.path(IN_DIR, "cluster_logrank.csv"))

METHOD_ORDER  <- c("NetComplex_rna", "NetComplex_prot", "NetComplex_multi", "NetComplex_multilayer")
METHOD_LABELS <- c("NetComplex_rna", "NetComplex_prot", "NetComplex_multi", "NetComplex_multilayer")
CV_ORDER      <- c("5%", "10%", "15%", "20%", "25%")
METHOD_COLORS <- c(
  NetComplex_rna        = "#888888",
  NetComplex_prot       = "#AAAAAA",
  NetComplex_multi      = "#74B9E7",
  NetComplex_multilayer = "#2A9D8F"
)
MAIN_OUTPUT_FILES <- c(
  PDAC = "Figure6E_PDAC_survival_cluster.pdf",
  KIRC = "Figure6F_KIRC_survival_cluster.pdf"
)
SUPPLEMENTARY_CANCERS <- c("GBM", "LUAD", "LUSC")
SUPPLEMENTARY_OUTPUT <- "SupplementaryFigure9.pdf"

dt[, method    := factor(method,    levels = METHOD_ORDER, labels = METHOD_LABELS)]
dt[, cv_filter := factor(cv_filter, levels = CV_ORDER)]
dt[, neg_log10_p := -log10(p_logrank)]

plot_cancer <- function(cancer_sel) {
  sub <- dt[cancer == cancer_sel]

  ggplot(sub, aes(x = cv_filter, y = neg_log10_p,
                  group = method, color = method)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = METHOD_COLORS, name = "Method") +
    scale_x_discrete(labels = CV_ORDER) +
    annotate("text", x = 0.6, y = -log10(0.05) + 0.07,
             label = "p = 0.05", size = 2.8, color = "grey50", hjust = 0) +
    labs(
      x     = "Top % complexes by CV",
      y     = expression(-log[10]("P-value")),
      title = cancer_sel
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", hjust = 0.5)
    )
}

for (cancer in names(MAIN_OUTPUT_FILES)) {
  p <- plot_cancer(cancer)
  fname <- file.path(OUT_DIR, MAIN_OUTPUT_FILES[[cancer]])
  ggsave(fname, p, width = 6, height = 4)
  cat("[saved]", fname, "\n")
}

p_supp <- ggplot(dt[cancer %in% SUPPLEMENTARY_CANCERS],
                 aes(x = cv_filter, y = neg_log10_p,
                     group = method, color = method)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_text(
    data = data.table(
      cancer = SUPPLEMENTARY_CANCERS,
      cv_filter = "5%",
      neg_log10_p = -log10(0.05) + 0.04,
      label = "P = 0.05"
    ),
    aes(x = cv_filter, y = neg_log10_p, label = label),
    inherit.aes = FALSE, color = "grey50", size = 2.6, hjust = 0
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_color_manual(values = METHOD_COLORS, name = "Method") +
  scale_x_discrete(labels = CV_ORDER) +
  facet_wrap(~ cancer, nrow = 1) +
  labs(
    x = "Top % complexes by CV",
    y = expression(-log[10]("P-value"))
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

supp_fname <- file.path(OUT_DIR, SUPPLEMENTARY_OUTPUT)
ggsave(supp_fname, p_supp, width = 9.2, height = 3.2)
cat("[saved]", supp_fname, "\n")

cat("[plot done]\n")
