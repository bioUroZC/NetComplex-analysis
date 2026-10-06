
library(data.table)
library(ggplot2)

base_dir <- "/proj/c.zihao/work3/03Survival"
font_size <- 9                                                                 
out_dir  <- file.path(base_dir, "01plots")
cancers  <- c("BLCA", "BRCA", "CHOL", "CRC", "GBM", "KIRC", "OV", "PAAD", "STAD")

methods <- c("raw_expr", "netcomplex", "ssgsea", "mean")
labels  <- c("Raw expression", "NetComplex", "ssGSEA", "Mean")
colours <- c("Raw expression" = "#8C8C8C", "NetComplex" = "#000000",
             "ssGSEA" = "#0072B2", "Mean" = "#D55E00")



res <- data.frame()
for (cancer in cancers) {
  wide <- as.data.frame(fread(file.path(base_dir, cancer, "02model", "lodo_cindex.csv")))
  for (m in methods) {
    res <- rbind(res, data.frame(cancer  = cancer,
                                 dataset = wide$dataset,
                                 method  = m,
                                 cindex  = wide[[m]]))
  }
}

res$cancer <- factor(res$cancer, levels = cancers)
res$method <- factor(res$method, levels = rev(methods), labels = rev(labels))



med <- aggregate(cindex ~ cancer + method, data = res, FUN = median)

base_med <- med[med$method == labels[1], c("cancer", "cindex")]



x_lo  <- min(0.50, min(res$cindex, na.rm = TRUE)) - 0.02
x_hi  <- max(res$cindex, na.rm = TRUE) + 0.02
lab_x <- x_hi + 0.04 * (x_hi - x_lo)                                               

med_head <- data.frame(cancer = factor(cancers, levels = cancers),
                       x = lab_x, y = length(methods) + 0.6)

p <- ggplot(res, aes(cindex, method)) +
  geom_vline(xintercept = 0.50, colour = "grey60", linewidth = 0.35, linetype = "dotted") +
  geom_vline(data = base_med, aes(xintercept = cindex),
             colour = "grey55", linewidth = 0.35, linetype = "dashed") +
  geom_point(aes(fill = method, shape = "One held-out dataset"),
             position = position_jitter(width = 0, height = 0.12, seed = 1),
             colour = "white", stroke = 0.3, size = 2, alpha = 0.55) +
  geom_point(data = med, aes(fill = method, shape = "Median across datasets"),
             colour = "white", stroke = 0.5, size = 3.4) +
  geom_text(data = med, aes(x = lab_x, label = sprintf("%.3f", cindex)),
            colour = "black", size = font_size, size.unit = "pt", hjust = 0) +
  geom_text(data = med_head, aes(x = x, y = y), label = "Median",
            colour = "black", size = font_size, size.unit = "pt", hjust = 0, inherit.aes = FALSE) +
  facet_wrap(~cancer, ncol = 3, axes = "all_x") +
  scale_fill_manual(values = colours, guide = "none") +
  scale_shape_manual(name = NULL, values = c("One held-out dataset" = 21, "Median across datasets" = 23)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey35", alpha = 1, size = c(2, 3.2)))) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = function(x) sprintf("%.1f", x),
                     expand = expansion(mult = 0.01)) +
  coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0.5, length(methods) + 0.5), clip = "off") +
  labs(title = "Leave-one-dataset-out survival prediction",
       x = "C-index values", y = NULL) +
  theme_classic(base_size = font_size, base_family = "sans") +
  theme(text               = element_text(size = font_size, colour = "black", face = "plain"),
        axis.line.x        = element_line(linewidth = 0.35, colour = "black"),
        axis.line.y        = element_blank(),
        axis.ticks.x       = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks.y       = element_blank(),
        axis.ticks.length  = unit(2.2, "pt"),
        axis.text          = element_text(size = font_size, colour = "black"),
        axis.text.y        = element_text(hjust = 1),
        panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey92"),
        strip.background   = element_blank(),
        strip.text         = element_text(size = font_size, colour = "black", hjust = 0,
                                          margin = margin(0, 0, 6, 0)),
        panel.spacing.x    = unit(40, "pt"),
        panel.spacing.y    = unit(14, "pt"),
        legend.position    = "top",
        legend.justification = "left",
        legend.location    = "plot",
        legend.margin      = margin(0, 0, 0, 0),
        legend.text        = element_text(size = font_size, colour = "black"),
        plot.title         = element_text(size = font_size, colour = "black"),
        plot.title.position = "plot",
        plot.margin        = margin(8, 40, 6, 6))

out_pdf <- file.path(out_dir, "supfigure7.pdf")
ggsave(out_pdf, p, width = 10, height = 8.4)
cat("saved", out_pdf, "\n")
