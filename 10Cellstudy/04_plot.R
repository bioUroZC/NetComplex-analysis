library(data.table)
library(ggplot2)

RESULT_DIR <- "/proj/c.zihao/work3/10Cellstudy/results"
PLOT_DIR <- "/proj/c.zihao/work3/10Cellstudy/plots"
SCORE_ROOT <- "/proj/c.zihao/work3/06singleCell"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

auc_long <- fread(file.path(RESULT_DIR, "percomplex_auc_long.csv"))
method_order <- c("ssgsea", "gsva", "netmean", "aucell", "ucell", "netComplex")
method_labels <- c("ssGSEA", "GSVA", "NetMean", "AUCell", "UCell", "NetComplex")

auc_long <- auc_long[set == "positive" & method %in% method_order]
auc_long[, method := factor(method, levels = method_order, labels = method_labels)]
auc_long[, label := sub("_complex$", "", sub("^Comp_", "", complex))]

complex_order <- auc_long[method == "NetComplex",
  .(mean_auc = mean(auc)), by = label]
setorder(complex_order, mean_auc)
auc_long[, label := factor(label, levels = complex_order$label)]

auc_plot <- ggplot(auc_long, aes(x = method, y = label, fill = auc)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(. ~ dataset) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0.5, limits = c(0, 1), name = "AUROC"
  ) +
  scale_y_discrete(labels = function(x) parse(text = sprintf("italic('%s')", x))) +
  labs(
    x = NULL,
    y = NULL,
    title = "Interferon-responsive complexes, scored separately"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

cases <- fread(file.path(RESULT_DIR, "case_complexes.csv"))
case_names <- unique(cases$complex)
case_labels <- setNames(
  sub("_complex$", "", sub("^Comp_", "", case_names)),
  case_names
)

conditions <- data.table(
  dataset = c(rep("GSE96583", 2), rep("GSE226572", 9)),
  sample = c(
    "GSM2560248", "GSM2560249",
    "GSM7079694", "GSM7079695", "GSM7079696",
    "GSM7079702", "GSM7079703", "GSM7079705",
    "GSM7079710", "GSM7079711", "GSM7079713"
  ),
  condition = c(
    "control", "stim",
    "control", "control", "stim",
    "control", "control", "stim",
    "control", "control", "stim"
  )
)
conditions[, condition := factor(condition, levels = c("control", "stim"))]

score_parts <- vector("list", nrow(conditions))
for (i in seq_len(nrow(conditions))) {
  dataset <- conditions$dataset[i]
  sample_id <- conditions$sample[i]
  condition <- conditions$condition[i]

  score_file <- file.path(
    SCORE_ROOT, dataset, "netcomplex", "bulk", sample_id, "netComplex.csv"
  )
  score_table <- fread(score_file)
  setnames(score_table, 1L, "complex")
  score_table <- score_table[complex %in% case_names]

  score_parts[[i]] <- melt(
    score_table,
    id.vars = "complex",
    variable.name = "cell",
    value.name = "netComplex_score"
  )[, `:=`(dataset = dataset, sample = sample_id, condition = condition)]
}

score_data <- rbindlist(score_parts)
score_data[, case := factor(
  case_labels[complex],
  levels = rev(case_labels[case_names])
)]

score_plot <- ggplot(score_data, aes(condition, netComplex_score, fill = condition)) +
  geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.4) +
  facet_grid(case ~ dataset, scales = "free_y") +
  scale_fill_manual(values = c(control = "#9ecae1", stim = "#e34a33")) +
  scale_x_discrete(labels = c(control = "Control cells", stim = "Stimulated cells")) +
  labs(
    x = NULL,
    y = "NetComplex score",
    title = "Selected-complex scores by condition"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", colour = NA)
  )

grDevices::pdf(
  file.path(PLOT_DIR, "Figure7CD.pdf"),
  width = 17,
  height = 10.5,
  useDingbats = FALSE
)
grid::grid.newpage()
print(auc_plot, vp = grid::viewport(x = 0.25, y = 0.47, width = 0.49, height = 0.91))
print(score_plot, vp = grid::viewport(x = 0.75, y = 0.47, width = 0.49, height = 0.91))
grid::grid.text("C", x = 0.012, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grid::grid.text("D", x = 0.512, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grDevices::dev.off()

member_signals <- fread(file.path(RESULT_DIR, "member_signal_summary.csv"))
member_signals <- member_signals[complex %in% case_names]
member_signals[, case := factor(
  case_labels[complex],
  levels = rev(case_labels[case_names])
)]

detection_data <- melt(
  member_signals,
  id.vars = c("dataset", "complex", "case", "gene"),
  measure.vars = c("detection_control", "detection_stim"),
  variable.name = "condition",
  value.name = "detection"
)
detection_data[, condition := factor(
  sub("detection_", "", condition),
  levels = c("control", "stim")
)]

detection_plot <- ggplot(
  detection_data,
  aes(condition, detection, group = gene, colour = gene)
) +
  geom_line(alpha = 0.65) +
  geom_point(size = 1.6) +
  facet_grid(case ~ dataset) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "Fraction of cells with member detected",
    colour = "Member",
    title = "Member detection is measured directly"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92", colour = NA)
  )

grDevices::pdf(
  file.path(PLOT_DIR, "SupplementaryFigure12.pdf"),
  width = 8.5,
  height = 5.1,
  useDingbats = FALSE
)
grid::grid.newpage()
print(detection_plot, vp = grid::viewport(x = 0.5, y = 0.47, width = 0.96, height = 0.91))
grid::grid.text("A", x = 0.012, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grDevices::dev.off()

cat("Wrote final single-cell PDFs to ", PLOT_DIR, "\n", sep = "")
