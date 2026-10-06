
library(ggplot2)

BASE     <- "/proj/c.zihao/work3/01protein/01random"
IN_FILE  <- file.path(BASE, "results", "concord.csv")
OUT_CSV  <- file.path(BASE, "results", "summary.csv")

CANCERS <- c("BRCA", "COAD", "GBM", "KIRC", "LUAD", "LUSC", "OV", "PDAC", "UCEC", "HNSCC")
METHODS <- c("netcomplex", "mean", "min", "zscore", "netmean",
             "perpas_adapted", "gsva", "ssgsea", "plage")
LABELS  <- c(netcomplex = "NetComplex", mean = "Mean", min = "Min", zscore = "Z-score",
             netmean = "NetMean", perpas_adapted = "PerPAS", gsva = "GSVA",
             ssgsea = "ssGSEA", plage = "PLAGE")
PALETTE <- c(
  netcomplex = "#000000", ssgsea = "#0072B2", mean = "#D55E00", min = "#009E73",
  zscore = "#CC79A7", netmean = "#E69F00", perpas_adapted = "#56B4E9",
  gsva = "#6A3D9A", plage = "#666666"
)
METRICS <- c(mantel_r = "Mantel r", procrustes_r = "Symmetric Procrustes r")

res <- read.csv(IN_FILE, stringsAsFactors = FALSE)
res <- res[res$n_pc == 3, ]
n_sets  <- length(unique(res$method[res$annotation == "random"]))
cancers <- CANCERS[CANCERS %in% res$cancer]
methods <- METHODS[METHODS %in% res$target_method]


per_cancer <- aggregate(cbind(mantel_r, procrustes_r) ~ target_method + cancer + method + annotation,
                        data = res, FUN = mean, na.rm = TRUE)

overall <- aggregate(cbind(mantel_r, procrustes_r) ~ target_method + method + annotation,
                     data = per_cancer, FUN = mean, na.rm = TRUE)
overall$cancer <- "overall"

values <- rbind(per_cancer, overall)
scopes <- c(cancers, "overall")


rows <- list()
for (scope in scopes) {
  for (m in methods) {
    for (metric in names(METRICS)) {
      d <- values[values$cancer == scope & values$target_method == m, ]
      real <- d[[metric]][d$annotation == "real"]
      null <- d[[metric]][d$annotation == "random"]
      if (length(real) != 1 || length(null) < 2) next

      null_sd <- sd(null)
      if (is.na(null_sd) || null_sd == 0) z <- NA_real_ else z <- (real - mean(null)) / null_sd

      row <- data.frame(
        scope         = scope,
        target_method = m,
        metric        = metric,
        real          = round(real, 5),
        null_mean     = round(mean(null), 5),
        null_sd       = round(null_sd, 5),
        null_min      = round(min(null), 5),
        null_max      = round(max(null), 5),
        null_q2.5     = round(unname(quantile(null, 0.025)), 5),
        null_q97.5    = round(unname(quantile(null, 0.975)), 5),
        n_null        = length(null),
        z             = round(z, 3),
        emp_p         = round((1 + sum(null >= real)) / (1 + length(null)), 4)
      )

      if (row$real > row$null_max) {
        row$verdict <- "WIN"
      } else if (row$real < row$null_min) {
        row$verdict <- "LOSE"
      } else if (row$real >= row$null_q2.5 && row$real <= row$null_q97.5) {
        row$verdict <- "TIE"
      } else {
        row$verdict <- "AMBIGUOUS"
      }
      rows[[length(rows) + 1]] <- row
    }
  }
}
summary_tab <- do.call(rbind, rows)
write.csv(summary_tab, OUT_CSV, row.names = FALSE)

cat(sprintf("\n=== overall, %d random sets, CV fractions averaged ===\n", n_sets))
print(summary_tab[summary_tab$scope == "overall",
                  c("target_method", "metric", "real", "null_mean", "z", "emp_p", "verdict")],
      row.names = FALSE)


curve_data <- aggregate(cbind(mantel_r, procrustes_r) ~ target_method + method + annotation + cv_top_fraction,
                        data = res, FUN = mean, na.rm = TRUE)
curve_data$label <- factor(LABELS[curve_data$target_method], levels = LABELS[methods])

theme_pub <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(2.2, "pt"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black", size = base_size),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
      panel.spacing.x = unit(12, "pt"),
      panel.spacing.y = unit(10, "pt"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.key.width = unit(26, "pt"),
      legend.key.height = unit(10, "pt"),
      legend.text = element_text(size = base_size - 0.5),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-6, 0, 0, 0),
      legend.box = "horizontal",
      plot.margin = margin(6, 8, 2, 4)
    )
}

plot_curves <- function(metric, keep = methods) {
  d <- curve_data[curve_data$target_method %in% keep, ]
  d$label <- factor(LABELS[d$target_method], levels = LABELS[keep])
  d$value <- d[[metric]]
  real <- d[d$annotation == "real", ]
  rand <- d[d$annotation == "random", ]

  band <- aggregate(value ~ label + cv_top_fraction, data = rand, FUN = mean)
  names(band)[3] <- "mid"
  band$lo <- aggregate(value ~ label + cv_top_fraction, data = rand, FUN = quantile, probs = 0.025)$value
  band$hi <- aggregate(value ~ label + cv_top_fraction, data = rand, FUN = quantile, probs = 0.975)$value
  ggplot() +
    geom_ribbon(data = band, aes(cv_top_fraction, ymin = lo, ymax = hi, fill = "2.5% to 97.5% interval"),
                colour = NA, show.legend = c(fill = TRUE, colour = FALSE, linetype = FALSE)) +
    geom_point(data = rand, aes(cv_top_fraction, value),
               position = position_jitter(width = 0.004, height = 0, seed = 1),
               size = 0.55, colour = "#8A8A8A", alpha = 0.45, stroke = 0) +
    geom_line(data = band, aes(cv_top_fraction, mid, linetype = "Mean of random complexes"),
              colour = "#4D4D4D", linewidth = 0.45) +
    geom_line(data = real, aes(cv_top_fraction, value, colour = target_method, group = label),
              linewidth = 0.9) +
    geom_point(data = real, aes(cv_top_fraction, value, colour = target_method),
               size = 1.7, show.legend = FALSE) +
    facet_wrap(~label, ncol = 3) +
    scale_x_continuous(breaks = sort(unique(d$cv_top_fraction)),
                       labels = paste0(100 * sort(unique(d$cv_top_fraction)), "%"),
                       expand = expansion(mult = 0.06)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    scale_fill_manual(values = c("2.5% to 97.5% interval" = "#E3E3E3"), name = NULL) +
    scale_colour_manual(values = PALETTE, breaks = keep, labels = LABELS[keep], name = NULL) +
    scale_linetype_manual(values = c("Mean of random complexes" = "dashed"), name = NULL) +
    guides(fill = guide_legend(order = 3),
           linetype = guide_legend(order = 2,
                                   override.aes = list(colour = "#4D4D4D", linewidth = 0.55,
                                                       linetype = "dashed", fill = NA)),
           colour = guide_legend(order = 1,
                                 override.aes = list(linewidth = 0.9, linetype = "solid",
                                                     fill = NA, shape = NA))) +
    labs(x = "Top-CV fraction", y = METRICS[[metric]]) +
    theme_pub()
}

MANTEL_MAIN  <- c("netcomplex", "mean", "ssgsea")
MANTEL_OTHER <- c("min", "zscore", "netmean", "perpas_adapted", "gsva", "plage")
OUT_MANTEL_MAIN  <- file.path(BASE, "results", "figure2c.pdf")
OUT_MANTEL_OTHER <- file.path(BASE, "results", "supfigure2.pdf")

ggsave(OUT_MANTEL_MAIN, plot_curves("mantel_r", MANTEL_MAIN),
       width = 7.2, height = 3.35)
ggsave(OUT_MANTEL_OTHER, plot_curves("mantel_r", MANTEL_OTHER),
       width = 7.2, height = 5.7)

cat("\nSummary ->", OUT_CSV, "\n")
cat("Figure  ->", OUT_MANTEL_MAIN, "\n")
cat("Figure  ->", OUT_MANTEL_OTHER, "\n")
