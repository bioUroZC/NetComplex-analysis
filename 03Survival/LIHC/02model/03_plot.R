
library(data.table)
library(ggplot2)

cancer    <- "LIHC"
model_dir <- "/proj/c.zihao/work3/03Survival/LIHC/02model"

methods <- c("raw_expr", "netcomplex", "ssgsea", "mean")
labels  <- c("raw expression  (baseline)", "NetComplex", "ssGSEA", "Mean")

res <- as.data.frame(fread(file.path(model_dir, "lodo_results.csv")))
res <- res[res$cancer == cancer & res$feature_set %in% methods, ]



wide <- data.frame(dataset = sort(unique(res$dataset)))
for (i in seq_len(nrow(wide))) {
  rows <- res[res$dataset == wide$dataset[i], ]
  wide$n_test[i]        <- max(rows$n_test)
  wide$n_events_test[i] <- max(rows$n_events_test)
  for (m in c("ssgsea", "netcomplex", "mean", "raw_expr")) {
    wide[i, m] <- rows$cindex[rows$feature_set == m]
  }
}
fwrite(wide, file.path(model_dir, "lodo_cindex.csv"))



res$method <- factor(res$feature_set, levels = rev(methods), labels = rev(labels))

med <- data.frame(method = levels(res$method))
for (i in seq_len(nrow(med))) {
  med$cindex[i] <- median(res$cindex[res$method == med$method[i]], na.rm = TRUE)
}
med$method <- factor(med$method, levels = levels(res$method))
med$colour <- ifelse(med$method == "NetComplex", "#C62F2F", "grey35")
base_med   <- med$cindex[med$method == labels[1]]

scored   <- wide$dataset[!is.na(wide$raw_expr)]
subtitle <- paste(length(scored), "datasets,",
                  sum(wide$n_test[wide$dataset %in% scored]), "samples,",
                  sum(wide$n_events_test[wide$dataset %in% scored]), "events")

x_lo  <- min(0.50, min(res$cindex, na.rm = TRUE) - 0.02)
x_hi  <- max(res$cindex, na.rm = TRUE) + 0.02
lab_x <- x_hi + 0.06 * (x_hi - x_lo)                                      

p <- ggplot(res, aes(cindex, method)) +
  geom_vline(xintercept = c(0.50, base_med), colour = "grey45", linewidth = 0.35) +
  geom_point(aes(shape = "one held-out dataset"),
             colour = "grey55", fill = NA, size = 2.1, stroke = 0.6) +
  geom_point(data = med, aes(colour = colour, shape = "median across datasets"),
             size = 3.4) +
  geom_text(data = med, aes(x = lab_x, label = sprintf("%.3f", cindex)),
            colour = "grey25", size = 3.1, hjust = 1) +
  annotate("text", x = lab_x, y = length(methods) + 0.75, label = "median",
           colour = "grey40", size = 3.1, hjust = 1, fontface = "bold") +
  scale_shape_manual(name = NULL,
                     values = c("one held-out dataset" = 21, "median across datasets" = 18),
                     breaks = c("one held-out dataset", "median across datasets")) +
  scale_colour_identity() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand = expansion(mult = 0.02)) +
  coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0.5, length(methods) + 0.5),
                  clip = "off") +
  labs(title = paste(cancer, "- leave-one-dataset-out survival prediction"),
       subtitle = subtitle, x = "C-index on the held-out dataset", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.line.x      = element_line(colour = "grey30", linewidth = 0.4),
        axis.ticks.x     = element_line(colour = "grey30", linewidth = 0.4),
        axis.text.y      = element_text(colour = "grey15", hjust = 0),
        legend.position  = "top",
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(colour = "grey35", hjust = 0.5),
        plot.margin      = margin(10, 58, 10, 10))

ggsave(file.path(model_dir, "lodo_cindex.pdf"), p, width = 6.4, height = 4)
cat("saved lodo_cindex.csv and lodo_cindex.pdf\n")
