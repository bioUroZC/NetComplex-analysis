library(ggplot2)

COUNT_DIR <- "/proj/c.zihao/work3/01protein/count"

NET_ORDER  <- c("HuRI", "string", "biogrid")
NET_LABELS <- c(HuRI = "HuRI", string = "STRING", biogrid = "BioGRID")
NET_COLORS <- c(HuRI = "#D55E00", string = "#0072B2", biogrid = "#009E73")

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0, face = "bold", size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, size = base_size - 1.5, colour = "#4B5563"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.size = unit(9, "pt"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

as_net_factor <- function(x) factor(NET_LABELS[x], levels = NET_LABELS[NET_ORDER])

deg <- read.csv(file.path(COUNT_DIR, "network_degradation.csv"), stringsAsFactors = FALSE)

CONSTANT_COLS <- c(
  "frac_nonisolated",
  "mean_degree_all_genes",
  "mean_degree_nonisolated",
  "frac_corum_nonisolated",
  "frac_complexes_all_isolated",
  "mean_frac_isolated_members"
)
for (col in CONSTANT_COLS) {
  spread <- tapply(deg[[col]], deg$network, function(v) diff(range(v)))
  if (any(spread > 1e-9)) {
    stop("Column varies across cohorts, panel would hide that: ", col)
  }
}

u <- deg[!duplicated(deg$network), ]
rownames(u) <- u$network
u <- u[NET_ORDER, ]

cplxC <- data.frame(
  network = rep(NET_ORDER, times = 2),
  metric = factor(
    rep(c(
      "Complexes with every member isolated",
      "Isolated members per complex (mean)"
    ), each = 3),
    levels = c(
      "Complexes with every member isolated",
      "Isolated members per complex (mean)"
    )
  ),
  value = c(u$frac_complexes_all_isolated, u$mean_frac_isolated_members),
  n = c(u$n_complexes_all_isolated, rep(NA_integer_, 3))
)
cplxC$net <- as_net_factor(cplxC$network)
cplxC$label <- ifelse(
  is.na(cplxC$n),
  sprintf("%.1f%%", 100 * cplxC$value),
  sprintf("%.1f%%\n(%d / 4,923)", 100 * cplxC$value, cplxC$n)
)

pC <- ggplot(cplxC, aes(x = net, y = value, fill = network)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = label), vjust = -0.25, size = 2.6, colour = "#1F2937", lineheight = 0.95) +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_manual(values = NET_COLORS, guide = "none") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.72),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Fraction",
    title = "PPI coverage of CORUM complex members",
    subtitle = "Isolated nodes retain their rank-normalized scores during propagation"
  ) +
  theme_pub() +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8.5))

ggsave(
  file.path(COUNT_DIR, "Figure_complex_level_consequence.pdf"),
  pC,
  width = 6.4,
  height = 4.6,
  device = cairo_pdf
)
