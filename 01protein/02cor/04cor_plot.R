
library(ggplot2)

IN_FILE  <- "/proj/c.zihao/work3/01protein/02cor/cor_complex_level.csv"
OUT_SUMMARY <- "/proj/c.zihao/work3/01protein/02cor/cor_diagonal_summary.csv"
OUT_TESTS   <- "/proj/c.zihao/work3/01protein/02cor/cor_diagonal_tests.csv"
OUT_FIG     <- "/proj/c.zihao/work3/01protein/02cor/supfigure1.pdf"

METHODS <- c("netcomplex", "mean", "min", "zscore", "netmean",
             "perpas_adapted", "gsva", "ssgsea", "plage")
METHOD_LABELS <- c(
  netcomplex = "NetComplex", mean = "Mean", min = "Min", zscore = "Z-score",
  netmean = "NetMean", perpas_adapted = "PerPAS", gsva = "GSVA",
  ssgsea = "ssGSEA", plage = "PLAGE"
)

if (!file.exists(IN_FILE)) stop("no correlation table; run 02cor_calc.R first: ", IN_FILE)
res <- read.csv(IN_FILE, stringsAsFactors = FALSE)

per_cohort <- aggregate(spearman_r ~ cancer + rna_method + protein_method,
                        data = res, FUN = mean, na.rm = TRUE)
overall <- aggregate(spearman_r ~ rna_method + protein_method,
                     data = per_cohort, FUN = mean, na.rm = TRUE)
write.csv(transform(overall, spearman_r = round(spearman_r, 5)), OUT_SUMMARY, row.names = FALSE)

cat("=== mean per-complex Spearman r (cohorts averaged), rows=RNA, cols=protein ref ===\n")
mat <- matrix(NA_real_, nrow = length(METHODS), ncol = length(METHODS),
              dimnames = list(METHOD_LABELS[METHODS], METHOD_LABELS[METHODS]))
for (i in seq_len(nrow(overall))) {
  mat[METHOD_LABELS[overall$rna_method[i]], METHOD_LABELS[overall$protein_method[i]]] <- overall$spearman_r[i]
}
print(round(mat, 4))

test_rows <- list()
for (pm in METHODS) {
  d <- res[res$protein_method == pm, c("cancer", "complex", "rna_method", "spearman_r")]
  wide <- reshape(d, idvar = c("cancer", "complex"), timevar = "rna_method", direction = "wide")
  colnames(wide) <- sub("^spearman_r\\.", "", colnames(wide))

  for (other in setdiff(METHODS, pm)) {
    paired <- wide[, c(pm, other)]
    paired <- paired[complete.cases(paired), ]
    if (nrow(paired) < 10) next
    wt <- suppressWarnings(wilcox.test(paired[[pm]], paired[[other]], paired = TRUE))
    test_rows[[length(test_rows) + 1]] <- data.frame(
      protein_method = pm,
      diagonal_rna = pm,
      compared_to_rna = other,
      n_paired_complexes = nrow(paired),
      median_diagonal_r = round(median(paired[[pm]]), 5),
      median_other_r = round(median(paired[[other]]), 5),
      median_diff = round(median(paired[[pm]] - paired[[other]]), 5),
      wilcoxon_p = signif(wt$p.value, 4),
      diagonal_wins = median(paired[[pm]]) > median(paired[[other]]),
      stringsAsFactors = FALSE
    )
  }
}
test_tab <- do.call(rbind, test_rows)
write.csv(test_tab, OUT_TESTS, row.names = FALSE)

cat("\n=== diagonal-preference paired tests (paired by cancer x complex) ===\n")
for (i in seq_len(nrow(test_tab))) {
  row <- test_tab[i, ]
  cat(sprintf("protein=%-11s  %s vs %-11s  n=%5d  median_r %.4f vs %.4f  diff=%+.4f  p=%.3g  -> diagonal %s\n",
              METHOD_LABELS[row$protein_method], METHOD_LABELS[row$diagonal_rna],
              METHOD_LABELS[row$compared_to_rna], row$n_paired_complexes,
              row$median_diagonal_r, row$median_other_r, row$median_diff, row$wilcoxon_p,
              if (row$diagonal_wins) "WINS" else "LOSES"))
}

n_diag_wins <- sum(test_tab$diagonal_wins)
n_diag_best <- sum(sapply(METHODS, function(pm) {
  col <- overall[overall$protein_method == pm, ]
  col$rna_method[which.max(col$spearman_r)] == pm
}))
cat(sprintf("\nDiagonal method wins its own column in %d of %d pairwise comparisons.\n",
            n_diag_wins, nrow(test_tab)))
cat(sprintf("Diagonal method is the column-wise maximum (cohort-averaged r) in %d of %d columns.\n",
            n_diag_best, length(METHODS)))
cat("(If the column-wise winner often changes with the protein reference,\n")
cat(" these protein-side scores are not method-neutral.)\n")

plot_df <- overall
plot_df$rna_label <- factor(METHOD_LABELS[plot_df$rna_method], levels = METHOD_LABELS[METHODS])
plot_df$protein_label <- factor(METHOD_LABELS[plot_df$protein_method], levels = METHOD_LABELS[METHODS])

p <- ggplot(plot_df, aes(protein_label, rna_label, fill = spearman_r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", spearman_r)), size = 2.6) +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027",
                       midpoint = mean(plot_df$spearman_r), name = "mean\nSpearman r") +
  labs(x = "protein-side reference", y = "RNA-side method",
       title = "Complex-level RNA-protein correlation depends on protein reference choice") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 11, face = "bold"))

ggsave(OUT_FIG, p, width = 8.5, height = 7.5)

cat("\nSummary       ->", OUT_SUMMARY, "\n")
cat("Diagonal tests ->", OUT_TESTS, "\n")
cat("Heatmap       ->", OUT_FIG, "\n")
