
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)

BASE_DIR <- "/proj/c.zihao/work3/06singleCell"
DATASETS <- c("GSE96583", "GSE226572")
OUTPUT_FILES <- list(
  GSE96583 = c(benchmark = "Figure6A_GSE96583_benchmark.pdf",
               ablation = "Figure6B_GSE96583_ablation_difference.pdf"),
  GSE226572 = c(benchmark = "Figure6C_GSE226572_benchmark.pdf",
                ablation = "Figure6D_GSE226572_ablation_difference.pdf")
)

METHOD_ORDER <- c(
  "netComplex", "ssgsea", "aucell", "ucell",
  "gsva", "netmean", "perpas_adapted"
)
METHOD_LABELS <- c(
  netComplex = "NetComplex",
  ssgsea = "ssGSEA",
  aucell = "AUCell",
  ucell = "UCell",
  gsva = "GSVA",
  netmean = "netmean",
  perpas_adapted = "PerPAS-adapted"
)

ABL_MAP <- tibble::tribble(
  ~method,           ~variant,
  "netComplex",     "NetComplex",
  "noRWR_bulk",     "noRWR",
  "degRandPPI_bulk", "RewiredPPI"
)

theme_pub <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(4, 8, 4, 4),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    strip.text = element_text(face = "bold")
  )

benchmark_panel <- function(pos, neg) {
  tab <- bind_rows(
    pos %>% mutate(set = "Stimulus-responsive complexes"),
    neg %>% mutate(set = "Negative-control complexes")
  ) %>%
    filter(method %in% METHOD_ORDER) %>%
    mutate(
      method_label = METHOD_LABELS[method],
      method_label = factor(method_label,
                            levels = METHOD_LABELS[METHOD_ORDER]),
      set = factor(set, levels = c("Stimulus-responsive complexes",
                                   "Negative-control complexes"))
    )

  dodge <- position_dodge(width = 0.70)
  ggplot(tab, aes(x = method_label, y = auc, fill = set, colour = set)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_boxplot(aes(group = interaction(method_label, set)),
                 position = dodge, width = 0.62, outlier.shape = NA,
                 linewidth = 0.45, alpha = 0.35) +
    geom_point(aes(group = set),
               position = position_jitterdodge(jitter.width = 0.10,
                                                dodge.width = 0.70),
               size = 1.15, alpha = 0.65) +
    scale_colour_manual(values = c(
      "Stimulus-responsive complexes" = "#B2182B",
      "Negative-control complexes" = "#4D4D4D"
    )) +
    scale_fill_manual(values = c(
      "Stimulus-responsive complexes" = "#B2182B",
      "Negative-control complexes" = "#4D4D4D"
    )) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Complex-level AUROC",
         title = "Curated positive and negative-control complexes") +
    theme_pub +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)
    )
}

ablation_differences <- function(df, set_label) {
  wide <- df %>%
    filter(method %in% ABL_MAP$method) %>%
    inner_join(ABL_MAP, by = "method") %>%
    select(variant, complex, auc) %>%
    pivot_wider(names_from = variant, values_from = auc)

  wide %>%
    transmute(
      complex,
      `RWR contribution` = NetComplex - noRWR,
      `PPI-topology contribution` = NetComplex - RewiredPPI
    ) %>%
    pivot_longer(-complex, names_to = "comparison", values_to = "difference") %>%
    mutate(set = set_label)
}

ablation_panel <- function(pos, neg) {
  delta <- bind_rows(
    ablation_differences(pos, "Stimulus-responsive complexes"),
    ablation_differences(neg, "Negative-control complexes")
  ) %>%
    filter(is.finite(difference)) %>%
    mutate(
      comparison = factor(comparison,
                          levels = c("RWR contribution", "PPI-topology contribution")),
      set = factor(set, levels = c("Stimulus-responsive complexes",
                                   "Negative-control complexes")),
      set_short = factor(
        recode(as.character(set),
               "Stimulus-responsive complexes" = "Stimulus-responsive",
               "Negative-control complexes" = "Negative control"),
        levels = c("Stimulus-responsive", "Negative control")
      )
    )

  ggplot(delta, aes(x = comparison, y = difference)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_boxplot(aes(x = set_short, group = set, fill = set, colour = set),
                 width = 0.62, outlier.shape = NA, linewidth = 0.45,
                 alpha = 0.35) +
    geom_point(aes(x = set_short, colour = set),
               position = position_jitter(width = 0.10, height = 0),
               size = 1.15, alpha = 0.65) +
    facet_wrap(~ comparison, nrow = 1) +
    scale_colour_manual(values = c(
      "Stimulus-responsive complexes" = "#B2182B",
      "Negative-control complexes" = "#4D4D4D"
    )) +
    scale_fill_manual(values = c(
      "Stimulus-responsive complexes" = "#B2182B",
      "Negative-control complexes" = "#4D4D4D"
    )) +
    coord_cartesian(ylim = c(-0.35, 0.85)) +
    labs(x = NULL, y = "AUROC difference (NetComplex - ablation)",
         title = "Bulk NetComplex ablation effects") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
}

for (gse in DATASETS) {
  pos_file <- file.path(BASE_DIR, gse, "04analysis/auc.csv")
  neg_file <- file.path(BASE_DIR, gse, "04analysis/auc_neg.csv")
  if (!file.exists(pos_file) || !file.exists(neg_file)) {
    warning("Skipping ", gse, ": missing AUROC input files.")
    next
  }

  pos <- read.csv(pos_file, stringsAsFactors = FALSE)
  neg <- read.csv(neg_file, stringsAsFactors = FALSE)
  p_benchmark <- benchmark_panel(pos, neg)
  p_ablation <- ablation_panel(pos, neg)

  ggsave(file.path(BASE_DIR, OUTPUT_FILES[[gse]][["benchmark"]]),
         p_benchmark, width = 5, height = 5, bg = "white")
  ggsave(file.path(BASE_DIR, OUTPUT_FILES[[gse]][["ablation"]]),
         p_ablation, width = 5, height = 5, bg = "white")
  cat("Saved ", OUTPUT_FILES[[gse]][["benchmark"]], " and ",
      OUTPUT_FILES[[gse]][["ablation"]], "\n", sep = "")
}

cat("\n=== Done ===\n")
