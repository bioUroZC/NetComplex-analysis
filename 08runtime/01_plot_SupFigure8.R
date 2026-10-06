

library(dplyr)
library(ggplot2)
library(patchwork)

ROOT <- "/proj/c.zihao/work3/08runtime"
N_REPEATS <- 3L
OUTFILE <- file.path(ROOT, "SupplementaryFigure10.pdf")

summarise_runtime <- function(result_dir) {
  files <- list.files(result_dir, pattern = "^runtime_.*_rep[0-9]+\\.csv$",
                      full.names = TRUE)
  if (!length(files)) {
    stop("No replicated runtime files found in ", result_dir)
  }

  raw <- bind_rows(lapply(files, function(path) {
    d <- read.csv(path, stringsAsFactors = FALSE)
    stem <- sub("^runtime_", "", basename(path))
    d$method <- sub("_rep[0-9]+\\.csv$", "", stem)
    d
  }))

  summary <- raw %>%
    group_by(method, sweep, value) %>%
    summarise(
      median_seconds = median(runtime_seconds),
      q25_seconds = quantile(runtime_seconds, 0.25),
      q75_seconds = quantile(runtime_seconds, 0.75),
      n_repeats = n(),
      .groups = "drop"
    )
  write.csv(summary, file.path(result_dir, "runtime_summary.csv"), row.names = FALSE)
  summary
}

bulk <- summarise_runtime(file.path(ROOT, "bulk/results"))
cell <- summarise_runtime(file.path(ROOT, "cell/results"))
required_columns <- c("method", "sweep", "value", "median_seconds",
                      "q25_seconds", "q75_seconds", "n_repeats")
for (d in list(bulk, cell)) {
  if (!all(required_columns %in% names(d))) {
    stop("A runtime summary does not contain the required columns.")
  }
  if (any(d$n_repeats < N_REPEATS)) {
    stop("Supplementary Figure 8 requires ", N_REPEATS,
         " timing repeats per point.")
  }
}

bulk_colors <- c(
  netcomplex = "#E63946", gsva = "#457B9D", ssgsea = "#1D3557",
  plage = "#2A9D8F", netmean = "#6A4C93", perpas_adapted = "#B5838D",
  mean = "#F4A261", min = "#E76F51", zscore = "#A8DADC"
)
bulk_labels <- c(
  netcomplex = "NetComplex", gsva = "GSVA", ssgsea = "ssGSEA",
  plage = "PLAGE", netmean = "NetMean", perpas_adapted = "PerPAS",
  mean = "Mean", min = "Min", zscore = "z-score"
)
cell_colors <- c(
  netcomplex = "#E63946", netmean = "#6A4C93", perpas_adapted = "#B5838D",
  gsva = "#457B9D", ssgsea = "#1D3557", aucell = "#F4A261", ucell = "#2A9D8F"
)
cell_labels <- c(
  netcomplex = "NetComplex", netmean = "NetMean", perpas_adapted = "PerPAS",
  gsva = "GSVA", ssgsea = "ssGSEA", aucell = "AUCell", ucell = "UCell"
)

make_plot <- function(data, sweep_name, x_lab, title, colors, labels) {
  d <- data %>% filter(.data$sweep == sweep_name)
  ggplot(d, aes(x = value, y = median_seconds, colour = method,
                fill = method, group = method)) +
    geom_ribbon(aes(ymin = q25_seconds, ymax = q75_seconds),
                alpha = 0.14, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.1) +
    scale_x_log10(breaks = sort(unique(d$value))) +
    scale_y_log10(labels = function(x) {
      vapply(x, format, character(1), scientific = FALSE, trim = TRUE)
    }) +
    scale_colour_manual(values = colors, labels = labels, name = NULL) +
    scale_fill_manual(values = colors, guide = "none") +
    labs(x = x_lab, y = "Runtime (s; median and IQR)", title = title) +
    theme_classic(base_size = 10) +
    theme(
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "top",
      legend.text = element_text(size = 8)
    )
}

p_a <- make_plot(bulk, "sample_size", "Number of samples", "Bulk sample-size scaling",
                 bulk_colors, bulk_labels)
p_b <- make_plot(bulk, "complex_count", "Number of complexes", "Bulk complex-count scaling",
                 bulk_colors, bulk_labels)
p_c <- make_plot(cell, "cell_count", "Number of cells", "Single-cell count scaling",
                 cell_colors, cell_labels)
p_d <- make_plot(cell, "complex_count", "Number of complexes", "Single-cell complex-count scaling",
                 cell_colors, cell_labels)

bulk_row <- ((p_a + p_b) + plot_layout(guides = "collect")) &
  theme(legend.position = "top")
cell_row <- ((p_c + p_d) + plot_layout(guides = "collect")) &
  theme(legend.position = "top")
figure <- (bulk_row / cell_row) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(OUTFILE, figure, width = 8, height = 8, bg = "white")
message("Wrote ", OUTFILE)
