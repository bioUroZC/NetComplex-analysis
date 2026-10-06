
library(ggplot2)

ROOT <- "/proj/c.zihao/work3/04propagation"
res <- read.csv(file.path(ROOT, "results/pca_step_sensitivity.csv"))
out_dir <- file.path(ROOT, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

METRICS <- c(mantel_r = "Mantel r", procrustes_r = "Procrustes r")
K_LEVELS <- sort(unique(res$k))
as_k <- function(x) factor(as.character(x), levels = as.character(K_LEVELS))

base_theme <- theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "grey95", color = "grey80"),
        legend.position = "top")

by_cohort <- aggregate(cbind(mantel_r, procrustes_r) ~ network + seed + cancer + k,
                       data = res, FUN = mean, na.rm = TRUE)
pooled <- aggregate(cbind(mantel_r, procrustes_r) ~ network + seed + k,
                    data = by_cohort, FUN = mean, na.rm = TRUE)

summarise_metric <- function(column) {
  rows <- lapply(K_LEVELS, function(k) {
    real <- pooled[pooled$network == "netcomplex" & pooled$k == k, column]
    null <- pooled[pooled$network == "degRandPPI" & pooled$k == k, column]
    if (length(real) != 1L || length(null) == 0L) return(NULL)
    null_mean <- mean(null)
    q <- quantile(null, c(0.025, 0.975), names = FALSE)
    data.frame(k = k, metric = unname(METRICS[column]),
               real = real, null_mean = null_mean, null_sd = sd(null),
               n_null = length(null),
               delta = real - null_mean,
               band_lo = q[1] - null_mean, band_hi = q[2] - null_mean,
               null_percentile = mean(null < real))
  })
  do.call(rbind, rows)
}

summary_all <- do.call(rbind, lapply(names(METRICS), summarise_metric))
summary_all$metric <- factor(summary_all$metric, levels = unname(METRICS))
write.csv(summary_all, file.path(ROOT, "results/step_effect_vs_null.csv"),
          row.names = FALSE)

decomp <- do.call(rbind, lapply(unname(METRICS), function(m) {
  s <- summary_all[summary_all$metric == m, ]
  base <- s$real[s$k == 0]
  rbind(
    data.frame(k = s$k, metric = m, component = "generic smoothing",
               value = s$null_mean - base),
    data.frame(k = s$k, metric = m, component = "real topology",
               value = s$real - s$null_mean)
  )
}))
decomp$component <- factor(decomp$component,
                           levels = c("generic smoothing", "real topology"))
decomp$metric <- factor(decomp$metric, levels = unname(METRICS))

decomp_totals <- aggregate(value ~ metric + k, data = decomp, FUN = sum)
decomp$k <- as_k(decomp$k)
decomp_totals$k <- as_k(decomp_totals$k)

p_decomp <- ggplot(decomp, aes(x = k, y = value, fill = component)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_col(width = 0.72) +
  geom_point(data = decomp_totals, aes(x = k, y = value), inherit.aes = FALSE,
             shape = 18, size = 2.4, color = "grey20") +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("generic smoothing" = "#BDBDBD",
                               "real topology" = "#08519C")) +
  labs(x = "Propagation steps K (alpha = 0.3 fixed)",
       y = "Gain over K = 0 (no propagation)", fill = NULL) +
  base_theme
ggsave(file.path(out_dir, "Figure5C_step_contribution_split.pdf"), p_decomp,
       width = 8.5, height = 4.6)

summary_all$k <- as_k(summary_all$k)

p_delta <- ggplot(summary_all, aes(x = k, y = delta, group = 1)) +
  geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = "grey80", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.8, color = "#08519C") +
  geom_point(size = 2.2, color = "#08519C") +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  labs(x = "Propagation steps K (alpha = 0.3 fixed)",
       y = "Topology-specific concordance gain") +
  base_theme
ggsave(file.path(out_dir, "Figure5B_step_delta_vs_null.pdf"), p_delta,
       width = 8.5, height = 4.4)

abs_long <- do.call(rbind, lapply(names(METRICS), function(column) {
  agg <- aggregate(pooled[[column]], by = list(k = pooled$k, network = pooled$network),
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))
  data.frame(k = agg$k, network = agg$network, mean = agg$x[, "mean"],
             sd = agg$x[, "sd"], metric = unname(METRICS[column]))
}))
abs_long$sd[is.na(abs_long$sd)] <- 0
abs_long$metric <- factor(abs_long$metric, levels = unname(METRICS))
abs_long$k <- as_k(abs_long$k)
abs_long$network <- factor(
  abs_long$network,
  levels = c("netcomplex", "degRandPPI", "noRestart"),
  labels = c("real STRING (alpha = 0.3)", "rewired null (mean of 20)",
             "no restart (alpha = 0)"))
abs_long <- abs_long[!is.na(abs_long$network), ]

p_abs <- ggplot(abs_long, aes(x = k, y = mean, colour = network, group = network)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("real STRING (alpha = 0.3)" = "#08519C",
                                 "rewired null (mean of 20)" = "#A63603",
                                 "no restart (alpha = 0)" = "#6A51A3")) +
  labs(x = "Propagation steps K",
       y = "Concordance pooled over cohorts (PC = 3)", colour = NULL) +
  base_theme
ggsave(file.path(out_dir, "SupplementaryFigure8.pdf"), p_abs,
       width = 9, height = 4.8)

cat("=== real vs rewired null, pooled over cohorts ===\n")
print(summary_all[, c("k", "metric", "real", "null_mean", "delta",
                      "band_lo", "band_hi", "null_percentile")],
      digits = 4, row.names = FALSE)

cat("\n=== contribution split at the converged solution (K =",
    max(K_LEVELS), ") ===\n")
split_tbl <- do.call(rbind, lapply(unname(METRICS), function(m) {
  d <- decomp[decomp$metric == m & decomp$k == as.character(max(K_LEVELS)), ]
  generic <- d$value[d$component == "generic smoothing"]
  topo    <- d$value[d$component == "real topology"]
  data.frame(metric = m, generic = generic, topology = topo,
             total = generic + topo,
             topology_share = topo / (generic + topo))
}))
print(split_tbl, digits = 4, row.names = FALSE)

cat("\nFigures ->", out_dir, "\n")
