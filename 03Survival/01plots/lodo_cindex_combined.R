
library(data.table)
library(ggplot2)
library(patchwork)

base_dir <- "/proj/c.zihao/work3/03Survival"
font_size <- 9                                                                 
out_dir  <- file.path(base_dir, "01plots")
cancers  <- c("LUAD", "LIHC", "PRAD")
tags     <- c("D", "E", "F")                                                

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

plots <- list()
for (i in seq_along(cancers)) {
  cancer <- cancers[i]
  d   <- res[res$cancer == cancer, ]
  m   <- med[med$cancer == cancer, ]
  b   <- base_med$cindex[base_med$cancer == cancer]

  p <- ggplot(d, aes(cindex, method)) +
    geom_vline(xintercept = 0.50, colour = "grey60", linewidth = 0.35, linetype = "dotted") +
    geom_vline(xintercept = b, colour = "grey55", linewidth = 0.35, linetype = "dashed") +
    geom_point(aes(fill = method, shape = "One held-out dataset"),
               position = position_jitter(width = 0, height = 0.12, seed = 1),
               colour = "white", stroke = 0.3, size = 2, alpha = 0.55) +
    geom_point(data = m, aes(fill = method, shape = "Median across datasets"),
               colour = "white", stroke = 0.5, size = 3.4) +
    geom_text(data = m, aes(x = lab_x, label = sprintf("%.3f", cindex)),
              colour = "black", size = font_size, size.unit = "pt", hjust = 0) +
    annotate("text", x = lab_x, y = length(methods) + 0.6, label = "Median",
             colour = "black", size = font_size, size.unit = "pt", hjust = 0) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_shape_manual(name = NULL, values = c("One held-out dataset" = 21, "Median across datasets" = 23)) +
    guides(shape = guide_legend(override.aes = list(fill = "grey35", alpha = 1, size = c(2, 3.2)))) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = function(x) sprintf("%.1f", x),
                       expand = expansion(mult = 0.01)) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0.5, length(methods) + 0.5), clip = "off") +
    labs(title = paste0(tags[i], "    ", cancer, ": leave-one-dataset-out survival prediction"),
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
          legend.text        = element_text(size = font_size, colour = "black"),
          plot.title         = element_text(size = font_size, colour = "black",
                                            margin = margin(0, 0, 6, 0)),
          plot.title.position = "plot",
          plot.margin        = margin(6, 40, 6, 6))

  if (i > 1) p <- p + theme(axis.text.y = element_blank())
  plots[[i]] <- p
}

fig <- wrap_plots(plots, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

out_pdf <- file.path(out_dir, "figure3def.pdf")
ggsave(out_pdf, fig, width = 11, height = 3.4)
cat("saved", out_pdf, "\n")
