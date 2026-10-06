
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

BASE_DIR   <- "/proj/c.zihao/work3/05independence"
RESULT_DIR <- file.path(BASE_DIR, "results")
PLOT_DIR   <- file.path(BASE_DIR, "plots")
BENCH_DIR  <- file.path(RESULT_DIR, "bench")
NC_DIR     <- file.path(RESULT_DIR, "netcomplex", "complex_score")

LABELS <- c("frac0.20", "frac0.40", "frac0.60", "frac0.80", "frac1.00")
REF    <- "frac1.00"

METHOD_FILE <- c(
  mean           = "mean.csv",
  min            = "min.csv",
  zscore         = "zscore.csv",
  netmean        = "netmean.csv",
  perpas_adapted = "perpas_adapted.csv",
  gsva           = "gsva.csv",
  ssgsea         = "ssgsea.csv",
  plage          = "plage.csv"
)

METHOD_ORDER <- c("netcomplex", "mean", "min", "ssgsea",
                  "zscore", "gsva", "plage", "netmean", "perpas_adapted")
METHOD_LABELS <- c(
  netcomplex     = "NetComplex",
  mean           = "mean",
  min            = "min",
  ssgsea         = "ssGSEA",
  zscore         = "z-score",
  gsva           = "GSVA",
  plage          = "PLAGE",
  netmean        = "netmean",
  perpas_adapted = "PerPAS-adapted"
)
METHOD_COLOURS <- c(
  netcomplex     = "#B2182B",
  mean           = "#000000",
  min            = "#7F7F7F",
  ssgsea         = "#56B4E9",
  zscore         = "#E69F00",
  gsva           = "#0072B2",
  plage          = "#CC79A7",
  netmean        = "#009E73",
  perpas_adapted = "#D55E00"
)
INDEPENDENT <- c("netcomplex", "mean", "min", "ssgsea")

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line          = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks         = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length  = unit(1.4, "mm"),
      axis.text          = element_text(colour = "black", size = rel(1.0)),
      axis.title         = element_text(colour = "black", size = rel(1.05)),
      axis.title.x       = element_text(margin = margin(t = 4)),
      axis.title.y       = element_text(margin = margin(r = 4)),
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.title       = element_blank(),
      legend.text        = element_text(size = rel(0.95), colour = "black"),
      legend.key         = element_blank(),
      legend.key.width   = unit(1.35, "lines"),
      legend.key.height  = unit(0.7, "lines"),
      legend.margin      = margin(0, 0, 1, 0),
      legend.spacing.x   = unit(0.15, "cm"),
      plot.title         = element_blank(),
      plot.margin        = margin(3, 6, 3, 3),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

read_score <- function(method, label) {
  path <- if (identical(method, "netcomplex")) {
    file.path(NC_DIR, label, paste0(label, "_netcomplex_complex_score.csv"))
  } else {
    file.path(BENCH_DIR, method, label, METHOD_FILE[[method]])
  }
  if (!file.exists(path)) return(NULL)
  as.matrix(read.csv(path, row.names = 1L, check.names = FALSE))
}

spearman_flat <- function(a, b) {
  cols <- intersect(colnames(a), colnames(b))
  rows <- intersect(rownames(a), rownames(b))
  if (length(cols) < 1L || length(rows) < 5L) return(NA_real_)
  x <- as.numeric(a[rows, cols, drop = FALSE])
  y <- as.numeric(b[rows, cols, drop = FALSE])
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

mae_shared <- function(a, b) {
  cols <- intersect(colnames(a), colnames(b))
  rows <- intersect(rownames(a), rownames(b))
  if (length(cols) < 1L || length(rows) < 1L) return(NA_real_)
  x <- as.numeric(a[rows, cols, drop = FALSE])
  y <- as.numeric(b[rows, cols, drop = FALSE])
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  mean(abs(x[ok] - y[ok]))
}

methods_present <- c(
  if (dir.exists(NC_DIR)) "netcomplex" else character(0),
  names(METHOD_FILE)[dir.exists(file.path(BENCH_DIR, names(METHOD_FILE)))]
)
if (length(methods_present) == 0L) {
  stop("No scored methods under ", RESULT_DIR, " -- run 02cal/run_array.sh first.")
}

rows <- list()
for (method in methods_present) {
  ref <- read_score(method, REF)
  if (is.null(ref)) {
    warning("Missing reference for ", method); next
  }
  for (label in LABELS) {
    mat <- read_score(method, label)
    if (is.null(mat)) next
    rows[[length(rows) + 1L]] <- data.frame(
      method = method,
      label = label,
      fraction = as.numeric(sub("^frac", "", label)),
      n_samples = ncol(mat),
      n_shared = length(intersect(colnames(ref), colnames(mat))),
      spearman_flat = spearman_flat(ref, mat),
      mae = mae_shared(ref, mat),
      independence_class = if (method %in% INDEPENDENT) "independent" else "dependent",
      stringsAsFactors = FALSE
    )
  }
}

raw <- bind_rows(rows)
raw$method <- factor(raw$method, levels = intersect(METHOD_ORDER, unique(raw$method)))
raw$independence_class <- factor(raw$independence_class,
                                 levels = c("independent", "dependent"))
write.csv(raw, file.path(RESULT_DIR, "independence_summary.csv"), row.names = FALSE)

curve <- raw %>% filter(label != REF)
cat(sprintf("Loaded %d method x fraction rows (%d methods)\n",
            nrow(raw), n_distinct(raw$method)))
print(as.data.frame(
  curve %>%
    transmute(fraction, method,
              spearman = sprintf("%.4f", spearman_flat)) %>%
    pivot_wider(names_from = method, values_from = spearman) %>%
    arrange(fraction)
), row.names = FALSE)

present <- levels(droplevels(curve$method))
method_cols <- METHOD_COLOURS[present]

heat_order <- curve %>%
  filter(fraction == min(fraction)) %>%
  arrange(independence_class, desc(spearman_flat), method) %>%
  mutate(method_label = METHOD_LABELS[as.character(method)]) %>%
  pull(method_label)
heat <- curve %>%
  mutate(
    method_label = factor(METHOD_LABELS[as.character(method)],
                          levels = rev(heat_order)),
    class_label = factor(
      ifelse(independence_class == "independent",
             "sample-independent", "cohort-dependent"),
      levels = c("sample-independent", "cohort-dependent")
    ),
    frac_label = percent(fraction, accuracy = 1),
    txt_col = ifelse(spearman_flat >= 0.97, "white", "#1F2937")
  )
heat$frac_label <- factor(heat$frac_label,
                          levels = percent(sort(unique(heat$fraction)), accuracy = 1))

p_heat <- ggplot(heat, aes(x = frac_label, y = method_label, fill = spearman_flat)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", spearman_flat), colour = txt_col),
            size = 2.5, family = "sans") +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#F7F7F7", "#DCEFEA", "#8ACFC0", "#4C91BD", "#356FA5"),
    values = rescale(c(0.70, 0.82, 0.92, 0.97, 1.00)),
    limits = c(0.70, 1.00),
    breaks = c(0.70, 0.80, 0.90, 1.00),
    name = "Spearman correlation"
  ) +
  labs(x = "Cohort size (fraction of samples)", y = NULL) +
  facet_grid(class_label ~ ., scales = "free_y", space = "free_y") +
  theme_pub() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = rel(0.9), colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "#FAFBFC", colour = NA),
    strip.background = element_rect(fill = "#E8EEF5", colour = "#A8B7C8",
                                    linewidth = 0.4),
    strip.text.y = element_text(size = rel(0.85), colour = "#23374D")
  )

snap <- curve %>%
  filter(fraction == min(fraction)) %>%
  arrange(spearman_flat, desc(as.character(method))) %>%
  mutate(
    method_label = factor(METHOD_LABELS[as.character(method)],
                          levels = METHOD_LABELS[as.character(method)]),
    lbl = sprintf("%.3f", spearman_flat),
    txt_col = ifelse(spearman_flat >= 0.92, "inside", "outside")
  )

p_bar <- ggplot(snap,
                aes(x = method_label, y = spearman_flat, fill = method)) +
  geom_col(width = 0.70, colour = NA) +
  geom_col(width = 0.70, fill = NA, colour = "grey20", linewidth = 0.28) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey50",
             linewidth = 0.35) +
  geom_text(data = dplyr::filter(snap, txt_col == "inside"),
            aes(label = lbl),
            hjust = 1.15, size = 2.4, colour = "white", family = "sans") +
  geom_text(data = dplyr::filter(snap, txt_col == "outside"),
            aes(label = lbl),
            hjust = -0.12, size = 2.4, colour = "black", family = "sans") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL,
       y = sprintf("Spearman concordance at %s cohort",
                   percent(min(snap$fraction), accuracy = 1))) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = rel(1.05), colour = "black"),
    plot.margin = margin(3, 12, 3, 3)
  )

save_pdf <- function(path, plot, width, height) {
  ggsave(path, plot, width = width, height = height, units = "in", bg = "white")
}

save_pdf(file.path(PLOT_DIR, "Figure5D_independence_heatmap.pdf"), p_heat, 4.6, 4.0)
save_pdf(file.path(PLOT_DIR, "independence_bar.pdf"), p_bar, 3.4, 3.1)

cat("Wrote plots to ", PLOT_DIR, "\n", sep = "")
