library(data.table)
library(ggplot2)

ROOT <- "/proj/c.zihao/work3"
STUDY_DIR <- file.path(ROOT, "09TCGAstudy")
RESULT_DIR <- file.path(STUDY_DIR, "results")
PLOT_DIR <- file.path(STUDY_DIR, "plots")
SCORE_ROOT <- file.path(ROOT, "02TCGAnt", "02stringTest")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

selected <- fread(file.path(RESULT_DIR, "kirc_selected_case.csv"))
members <- fread(file.path(RESULT_DIR, "kirc_selected_member_genes.csv"))

net_raw <- fread(file.path(SCORE_ROOT, "netcomplex", "complex_score", "KIRC",
                           "KIRC_netcomplex_complex_score.csv"), check.names = FALSE)
net <- as.matrix(net_raw[, -1, with = FALSE])
storage.mode(net) <- "numeric"
rownames(net) <- net_raw[[1L]]

mean_raw <- fread(file.path(SCORE_ROOT, "bench", "mean", "KIRC", "mean.csv"), check.names = FALSE)
mean_score <- as.matrix(mean_raw[, -1, with = FALSE])
storage.mode(mean_score) <- "numeric"
rownames(mean_score) <- mean_raw[[1L]]

sample_ids <- colnames(net)
patient <- sub("_(01A|11A)$", "", sample_ids)
tumour <- grepl("_01A$", sample_ids)
normal <- grepl("_11A$", sample_ids)
patient_id <- sort(intersect(patient[tumour], patient[normal]))
pairs <- data.table(
  patient_id,
  tumour_id = unname(setNames(sample_ids[tumour], patient[tumour])[patient_id]),
  normal_id = unname(setNames(sample_ids[normal], patient[normal])[patient_id])
)

net_long <- data.table(
  patient = rep(pairs$patient_id, each = 2L),
  condition = factor(rep(c("Normal", "Tumour"), length(pairs$patient_id)),
                     levels = c("Normal", "Tumour")),
  value = as.numeric(c(rbind(net[selected$complex, pairs$normal_id],
                             net[selected$complex, pairs$tumour_id]))),
  panel = "NetComplex"
)
mean_long <- data.table(
  patient = rep(pairs$patient_id, each = 2L),
  condition = factor(rep(c("Normal", "Tumour"), length(pairs$patient_id)),
                     levels = c("Normal", "Tumour")),
  value = as.numeric(c(rbind(mean_score[selected$complex, pairs$normal_id],
                             mean_score[selected$complex, pairs$tumour_id]))),
  panel = "Member mean"
)
score_long <- rbind(net_long, mean_long)

case_label <- sub("^Comp_", "", selected$complex)
case_label <- gsub("_", " ", case_label, fixed = TRUE)
score_labels <- merge(data.table(
  panel = c("NetComplex", "Member mean"),
  label = c(
    sprintf("d[z] = %.2f\nFDR = %.2g", selected$net_dz, selected$net_fdr),
    sprintf("d[z] = %.2f\nFDR = %.2g", selected$mean_dz, selected$mean_fdr)
  ),
  condition = factor("Normal", levels = c("Normal", "Tumour"))
), score_long[, .(value = max(value) + 0.035 * diff(range(value))), by = panel], by = "panel", sort = FALSE)
fig_a <- ggplot(score_long, aes(condition, value, group = patient)) +
  geom_line(colour = "grey55", alpha = 0.28, linewidth = 0.35) +
  geom_point(aes(colour = condition), alpha = 0.62, size = 1.35) +
  geom_boxplot(aes(group = condition), width = 0.48, outlier.shape = NA,
               fill = NA, colour = "black", linewidth = 0.45) +
  geom_text(data = score_labels, aes(x = condition, y = value, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1.15, size = 3.3) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = c(Normal = "#277DA1", Tumour = "#D1495B"), guide = "none") +
  labs(
    title = paste("KIRC case:", case_label),
    subtitle = sprintf("KIRC paired tumour-normal samples (n = %d)", selected$n_pairs),
    x = NULL, y = "Score"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

decomposition <- fread(file.path(RESULT_DIR, "kirc_signal_decomposition.csv"))
drivers <- fread(file.path(RESULT_DIR, "kirc_signal_contributions.csv"))
external_input <- decomposition[component == "External PPI-neighbour input", score_delta]
top_drivers <- drivers[is_complex_member == FALSE][order(-rwr_contribution_to_score_delta)][1:5]
neighbour_label <- paste(
  "External PPI neighbours\n(top RWR contributors)",
  paste(sprintf("%s  (RNA delta %+.2f)", top_drivers$gene, top_drivers$raw_expression_delta), collapse = "\n"),
  sep = "\n"
)
member_label <- paste(
  "Complex members: direct RNA", 
  paste(sprintf("%s  (delta %+.2f)", members$gene, members$mean_delta), collapse = "\n"),
  sep = "\n"
)
diagram_nodes <- data.table(
  node = c("neighbours", "members", "rwr", "mean", "net"),
  x = c(1.2, 4.0, 4.0, 7.1, 7.1),
  y = c(4.0, 4.0, 1.7, 4.0, 1.7),
  label = c(
    neighbour_label,
    member_label,
    sprintf("RWR on STRING PPI\nadds neighbour context\nexternal input = %+.5f", external_input),
    sprintf("Member mean: weak\nd[z] = %.2f; FDR = %.2g", selected$mean_dz, selected$mean_fdr),
    sprintf("NetComplex: strong\nd[z] = %.2f; FDR = %.2g", selected$net_dz, selected$net_fdr)
  ),
  fill = c("#FADBD8", "#D6EAF8", "#FCF3CF", "#F2F3F4", "#D5F5E3")
)
diagram_edges <- data.table(
  x = c(1.95, 4.0, 4.75, 4.75), y = c(3.5, 3.35, 4.0, 1.7),
  xend = c(3.35, 4.0, 6.35, 6.35), yend = c(1.95, 2.35, 4.0, 1.7)
)
fig_c <- ggplot() +
  geom_curve(data = diagram_edges, aes(x = x, y = y, xend = xend, yend = yend),
             curvature = 0.12, linewidth = 0.75, colour = "grey35",
             arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed")) +
  geom_label(data = diagram_nodes, aes(x = x, y = y, label = label, fill = fill),
             linewidth = 0.4, size = 3.35, lineheight = 1.08, colour = "black",
             show.legend = FALSE) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.1, 8.25), ylim = c(0.45, 5.05), clip = "off") +
  labs(
    title = "Why can NetComplex detect this KIRC complex?",
    subtitle = "Member expression is weak, but positive RNA signals in connected PPI neighbours are propagated into the complex score."
  ) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0),
        plot.margin = margin(12, 10, 10, 12))

grDevices::pdf(file.path(PLOT_DIR, "Figure7AB.pdf"), width = 17, height = 6.1, useDingbats = FALSE)
grid::grid.newpage()
print(fig_a, vp = grid::viewport(x = 0.25, y = 0.47, width = 0.49, height = 0.91))
print(fig_c, vp = grid::viewport(x = 0.75, y = 0.47, width = 0.49, height = 0.91))
grid::grid.text("A", x = 0.012, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grid::grid.text("B", x = 0.512, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grDevices::dev.off()

cat("Wrote KIRC figures to: ", PLOT_DIR, "\n", sep = "")
