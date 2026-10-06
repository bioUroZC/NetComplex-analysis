library(ggplot2)

ROOT <- "/proj/c.zihao/work3/04parameter"
res <- read.csv(file.path(ROOT, "results/pca_alpha_sensitivity.csv"))
out_dir <- file.path(ROOT, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

by_cancer <- aggregate(cbind(mantel_r, procrustes_r) ~ cancer + alpha,
                       data = res, FUN = mean, na.rm = TRUE)
summary <- aggregate(cbind(mantel_r, procrustes_r) ~ alpha, data = by_cancer,
                     FUN = mean, na.rm = TRUE)
summary <- data.frame(
  alpha = summary$alpha,
  mantel_mean = summary$mantel_r,
  procrustes_mean = summary$procrustes_r
)

long <- rbind(
  data.frame(alpha = summary$alpha, metric = "Mantel r",
             mean = summary$mantel_mean),
  data.frame(alpha = summary$alpha, metric = "Symmetric Procrustes r",
             mean = summary$procrustes_mean)
)
long$metric <- factor(long$metric, levels = c("Mantel r", "Symmetric Procrustes r"))

p <- ggplot(long, aes(x = alpha, y = mean, colour = metric)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("Mantel r" = "#08519C",
                                 "Symmetric Procrustes r" = "#D55E00"),
                      name = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  labs(x = "RWR restart alpha", y = "Mean concordance (PC = 3)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "top"
  )
ggsave(file.path(out_dir, "Figure5A_alpha_pca_sensitivity.pdf"), p,
       width = 5, height = 5)
