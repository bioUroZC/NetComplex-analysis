
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

BASE_DIR    <- "/proj/c.zihao/work3/05noise"
RESULT_DIR  <- file.path(BASE_DIR, "results")
SCAN_DIR    <- file.path(RESULT_DIR, "scan")
PLOT_DIR    <- file.path(BASE_DIR, "plots")
METHOD_ORDER <- c("netcomplex", "RewiredPPI", "noRWR")
METHOD_LABELS <- c(
  netcomplex = "NetComplex",
  RewiredPPI = "RewiredPPI",
  noRWR      = "noRWR"
)
METHOD_COLOURS <- c(
  netcomplex = "#B2182B",
  noRWR      = "#4D4D4D",
  RewiredPPI = "#009E73"
)
METHOD_LINETYPES <- c(
  netcomplex = "solid",
  noRWR      = "solid",
  RewiredPPI = "dashed"
)
NETWORK_LABELS <- c(string = "STRING", biogrid = "BioGRID", HuRI = "HuRI")
NOISE_LEVELS <- c(0, 0.25, 0.5, 1, 2, 4)

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line        = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.45, colour = "black"),
      axis.text        = element_text(colour = "black", size = rel(0.95)),
      axis.title       = element_text(colour = "black", size = rel(1.0)),
      legend.position  = "top",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      legend.text      = element_text(size = rel(0.95)),
      legend.key       = element_blank(),
      legend.key.width = unit(1.15, "lines"),
      legend.margin    = margin(0, 0, 2, 0),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = rel(1.0),
                                      margin = margin(0, 0, 4, 0)),
      plot.title       = element_blank(),
      plot.margin      = margin(4, 8, 4, 4),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.minor   = element_blank()
    )
}

files <- list.files(SCAN_DIR, pattern = "_noise_scan\\.csv$",
                    full.names = TRUE)
if (length(files) == 0L) {
  stop("No scan results under ", SCAN_DIR,
       " -- run 02_*.py scoring scripts then 03_compute_recovery.py first.")
}
raw <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
if (!"RewiredPPI" %in% unique(raw$method)) {
  stop("No RewiredPPI rows found -- re-run 02_RewiredPPI.py.")
}
raw <- raw %>% filter(method != "randPPI")
present_methods <- intersect(METHOD_ORDER, unique(raw$method))
raw$method <- factor(raw$method, levels = present_methods)
raw$network_label <- factor(
  dplyr::recode(as.character(raw$network), !!!NETWORK_LABELS),
  levels = unname(NETWORK_LABELS[intersect(names(NETWORK_LABELS),
                                           unique(raw$network))])
)
cat(sprintf("Loaded %d rows from %d files (%d networks, %d cancers)\n",
            nrow(raw), length(files),
            length(unique(raw$network)), length(unique(raw$cancer))))

cohort_tbl <- raw %>%
  group_by(network, network_label, cancer, noise_level, method) %>%
  summarise(complex_recovery  = mean(complex_recovery,  na.rm = TRUE),
            geometry_recovery = mean(geometry_recovery, na.rm = TRUE),
            .groups = "drop")

summary_tbl <- cohort_tbl %>%
  group_by(network, network_label, noise_level, method) %>%
  summarise(
    n = n(),
    complex_recovery_se  = sd(complex_recovery,    na.rm = TRUE) / sqrt(n),
    geometry_recovery_se = sd(geometry_recovery,   na.rm = TRUE) / sqrt(n),
    complex_recovery     = mean(complex_recovery,  na.rm = TRUE),
    geometry_recovery    = mean(geometry_recovery, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_tbl %>% select(-network_label),
          file.path(RESULT_DIR, "noise_summary.csv"), row.names = FALSE)

cat("\n-- complex_recovery (mean +/- SE over cancers) --\n")
print(as.data.frame(
  summary_tbl %>%
    transmute(network, noise_level, method,
              complex_recovery = sprintf("%.4f+/-%.4f",
                                         complex_recovery, complex_recovery_se)) %>%
    pivot_wider(names_from = method, values_from = complex_recovery)
), row.names = FALSE)

wide <- raw %>%
  select(network, network_label, cancer, noise_level, seed, method,
         complex_recovery, geometry_recovery) %>%
  pivot_wider(names_from = method,
              values_from = c(complex_recovery, geometry_recovery))

contrast_pairs <- wide %>%
  mutate(
    vs_noRWR_complex       = complex_recovery_netcomplex  - complex_recovery_noRWR,
    vs_RewiredPPI_complex  = complex_recovery_netcomplex  - complex_recovery_RewiredPPI,
    vs_noRWR_geometry      = geometry_recovery_netcomplex - geometry_recovery_noRWR,
    vs_RewiredPPI_geometry = geometry_recovery_netcomplex - geometry_recovery_RewiredPPI
  )

contrast_cohort <- contrast_pairs %>%
  group_by(network, network_label, cancer, noise_level) %>%
  summarise(across(starts_with("vs_"), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

contrast <- contrast_cohort %>%
  group_by(network, network_label, noise_level) %>%
  summarise(
    across(starts_with("vs_"), list(mean = ~ mean(.x, na.rm = TRUE),
                                    se   = ~ sd(.x, na.rm = TRUE) / sqrt(n())),
           .names = "{.col}_{.fn}"),
    n_pairs = n(),
    .groups = "drop"
  ) %>%
  rename_with(~ sub("_mean$", "", .x), ends_with("_mean"))

write.csv(
  contrast %>%
    select(network, noise_level,
           vs_noRWR_complex, vs_noRWR_complex_se,
           vs_RewiredPPI_complex, vs_RewiredPPI_complex_se,
           vs_noRWR_geometry, vs_noRWR_geometry_se,
           vs_RewiredPPI_geometry, vs_RewiredPPI_geometry_se,
           n_pairs),
  file.path(RESULT_DIR, "noise_contrasts.csv"), row.names = FALSE
)

cat("\n-- NetComplex minus control (positive = NetComplex retains more) --\n")
print(as.data.frame(
  contrast %>%
    select(network, noise_level, vs_noRWR_complex, vs_RewiredPPI_complex,
           vs_noRWR_geometry, vs_RewiredPPI_geometry, n_pairs)
), row.names = FALSE)

n_networks <- nlevels(summary_tbl$network_label)
method_cols <- METHOD_COLOURS[as.character(present_methods)]
method_lty  <- METHOD_LINETYPES[as.character(present_methods)]
method_labs <- METHOD_LABELS[as.character(present_methods)]

plot_recovery <- function(mean_col, se_col, ylab) {
  plot_tbl <- summary_tbl %>%
    mutate(noise_x = match(noise_level, NOISE_LEVELS))
  p <- ggplot(plot_tbl,
              aes(x = noise_x, y = .data[[mean_col]],
                  colour = method, linetype = method,
                  fill = method, group = method)) +
    geom_ribbon(aes(ymin = .data[[mean_col]] - .data[[se_col]],
                    ymax = .data[[mean_col]] + .data[[se_col]]),
                alpha = 0.16, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.1, shape = 16) +
    scale_colour_manual(values = method_cols, labels = method_labs,
                        name = NULL) +
    scale_fill_manual(values = method_cols, labels = method_labs,
                      name = NULL) +
    scale_linetype_manual(values = method_lty, labels = method_labs,
                          name = NULL) +
    scale_x_continuous(breaks = seq_along(NOISE_LEVELS),
                       labels = NOISE_LEVELS) +
    scale_y_continuous(limits = c(0, 1.02),
                       breaks = seq(0, 1, by = 0.2),
                       expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "Noise level (\u00d7 per-gene SD)", y = ylab) +
    theme_pub() +
    guides(
      colour = guide_legend(nrow = 1,
                            override.aes = list(fill = NA, linewidth = 0.95, size = 2.2)),
      linetype = guide_legend(nrow = 1),
      fill = "none"
    )
  if (n_networks > 1L) p <- p + facet_wrap(~ network_label, nrow = 1)
  p
}

p_complex <- plot_recovery("complex_recovery", "complex_recovery_se",
                           "Complex recovery")
p_geom    <- plot_recovery("geometry_recovery", "geometry_recovery_se",
                           "Geometry recovery")

panel_w <- if (n_networks > 1L) 2.4 * n_networks + 1.2 else 3.6
panel_h <- 3.2

save_pdf <- function(path, plot, width, height) {
  ggsave(path, plot, width = width, height = height, bg = "white")
}

save_pdf(file.path(PLOT_DIR, "Figure5E_noise_complex_recovery.pdf"),
         p_complex, panel_w, panel_h)
save_pdf(file.path(PLOT_DIR, "Figure5F_noise_geometry_recovery.pdf"),
         p_geom, panel_w, panel_h)

cat("\nPlots ->", PLOT_DIR, "\n")
cat("=== Done ===\n")
