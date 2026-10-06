
library(ggplot2)

BASE    <- "/proj/c.zihao/work3/01protein/01random"
IN_FILE <- file.path(BASE, "results", "concord.csv")

CV_FRACTION     <- 0.05                                                          
EXCLUDE_CANCERS <- character(0)                                         

OUT_PDF <- file.path(BASE, "results", "supfigure3.pdf")

METHODS <- c("netcomplex", "mean", "min", "zscore", "netmean",
             "perpas_adapted", "gsva", "ssgsea", "plage")
LABELS  <- c(netcomplex = "NetComplex", mean = "Mean", min = "Min", zscore = "Z-score",
             netmean = "NetMean", perpas_adapted = "PerPAS", gsva = "GSVA",
             ssgsea = "ssGSEA", plage = "PLAGE")
METRICS <- c(mantel_r = "Mantel r")


res <- read.csv(IN_FILE, stringsAsFactors = FALSE)
res <- res[res$n_pc == 3, ]
res <- res[!(res$cancer %in% EXCLUDE_CANCERS), ]
res <- res[round(res$cv_top_fraction, 4) %in% round(CV_FRACTION, 4), ]
if (nrow(res) == 0) stop("no rows left for CV_FRACTION = ", paste(CV_FRACTION, collapse = ", "))

n_cancers <- length(unique(res$cancer))
n_sets    <- length(unique(res$method[res$annotation == "random"]))
methods   <- METHODS[METHODS %in% res$target_method]
cv_text   <- paste0("top ", paste0(100 * CV_FRACTION, "%", collapse = " + "), " CV")

values <- aggregate(cbind(mantel_r, procrustes_r) ~ target_method + cancer + method + annotation,
                    data = res, FUN = mean, na.rm = TRUE)


rank_sum_distribution <- function(n_random) {
  prob <- 1                                                               
  for (n in n_random) {
    new_prob <- rep(0, length(prob) + n + 1)
    for (r in 1:(n + 1)) {                                             
      idx <- (1:length(prob)) + r
      new_prob[idx] <- new_prob[idx] + prob / (n + 1)
    }
    prob <- new_prob
  }
  prob
}

stratified_p <- function(d, metric) {
  ranks <- c()
  n_random <- c()
  for (cancer in unique(d$cancer)) {
    real <- d[[metric]][d$cancer == cancer & d$annotation == "real"]
    rand <- d[[metric]][d$cancer == cancer & d$annotation == "random"]
    if (length(real) != 1 || length(rand) < 2) next                                   
    if (any(rand == real)) stop("tie between real and random in ", cancer)
    ranks    <- c(ranks, 1 + sum(rand < real))
    n_random <- c(n_random, length(rand))
  }
  prob <- rank_sum_distribution(n_random)
  s <- sum(ranks)
  p_high <- sum(prob[(s + 1):length(prob)])                         
  p_low  <- sum(prob[1:(s + 1)])                                    
  min(1, 2 * min(p_high, p_low))
}


plot_data <- data.frame()
p_labels  <- data.frame()
for (i in seq_along(methods)) {
  m <- methods[i]
  d <- values[values$target_method == m, ]
  for (metric in names(METRICS)) {
    plot_data <- rbind(plot_data, data.frame(
      label  = LABELS[[m]],
      metric = METRICS[[metric]],
      group  = ifelse(d$annotation == "real", "Real", "Random"),
      x      = i + ifelse(d$annotation == "real", 0.2, -0.2),                      
      value  = d[[metric]]
    ))
    p <- stratified_p(d, metric)
    p_labels <- rbind(p_labels, data.frame(
      label  = LABELS[[m]],
      metric = METRICS[[metric]],
      x      = i,
      y      = max(d[[metric]]),
      text   = paste0("p=", trimws(formatC(p, format = "g", digits = 2)))
    ))
  }
}
plot_data$group  <- factor(plot_data$group, levels = c("Random", "Real"))
plot_data$metric <- factor(plot_data$metric, levels = METRICS)
p_labels$metric  <- factor(p_labels$metric, levels = METRICS)


idx <- ave(seq_len(nrow(p_labels)), p_labels$metric, FUN = seq_along)
panel_max <- tapply(plot_data$value, plot_data$metric, max)
panel_min <- tapply(plot_data$value, plot_data$metric, min)
panel_span <- panel_max - panel_min
p_labels$y <- panel_max[as.character(p_labels$metric)] +
  (0.10 + 0.16 * (idx %% 2)) * panel_span[as.character(p_labels$metric)]

fmt_p <- function(p) {
  if (p < 0.001) {
    e <- floor(log10(p))
    m <- signif(p / 10^e, 2)
    if (m >= 9.95) { m <- 1; e <- e + 1 }
    sprintf("italic(p) == %.1f %%*%% 10^{%d}", m, e)
  } else if (p < 0.01) {
    sprintf("italic(p) == %.4f", signif(p, 2))
  } else if (p < 0.1) {
    sprintf("italic(p) == %.3f", p)
  } else {
    sprintf("italic(p) == %.2f", p)
  }
}
p_labels$text <- vapply(p_labels$text, function(s) fmt_p(as.numeric(sub("^p=", "", s))), "")

fig <- ggplot(plot_data, aes(x = x, y = value)) +
  geom_jitter(data = plot_data[plot_data$group == "Random", ],
              width = 0.06, height = 0, size = 0.28, colour = "#B0B0B0", alpha = 0.28, stroke = 0) +
  geom_boxplot(aes(group = interaction(label, group), fill = group),
               colour = "#333333", width = 0.32, outlier.shape = NA, linewidth = 0.3) +
  geom_point(data = plot_data[plot_data$group == "Real", ],
             colour = "#D55E00", size = 1.15) +
  geom_text(data = p_labels, aes(x = x, y = y, label = text),
            vjust = 0, size = 2.15, parse = TRUE, colour = "black") +
  scale_x_continuous(breaks = seq_along(methods), labels = LABELS[methods]) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  scale_fill_manual(values = c(Random = "#E6E6E6", Real = "#F6CFA8"),
                    labels = c(Random = "Random", Real = "Real CORUM"), name = NULL) +
  labs(x = NULL, y = "Mantel r") +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks.length = unit(2.2, "pt"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.title = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", size = 9),
        panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey90"),
        panel.spacing.y = unit(12, "pt"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.width = unit(16, "pt"),
        legend.key.height = unit(10, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 7.5, colour = "grey20", hjust = 0, lineheight = 1.05),
        plot.caption.position = "plot",
        plot.margin = margin(4, 10, 4, 4))

ggsave(OUT_PDF, fig, width = 7.2, height = 4.5)
embedFonts(OUT_PDF, options = "-dPDFSETTINGS=/prepress")
print(p_labels[, c("metric", "label", "text")], row.names = FALSE)
cat("Figure ->", OUT_PDF, "\n")
