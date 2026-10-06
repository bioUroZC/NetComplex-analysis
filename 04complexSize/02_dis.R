
library(ggplot2)

OUT_DIR    <- "/proj/c.zihao/work3/04complexSize"
COMPLEX_DB <- "/proj/c.zihao/work3/00data/corum/complexData.csv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SIZE_BIN_LEVELS <- c("1", "2", "3", "4", "5-7", "8-10", "11+")

cdb <- read.csv(COMPLEX_DB, stringsAsFactors = FALSE)
complex_size <- vapply(strsplit(cdb$Genes, ";"),
                       function(x) sum(nzchar(trimws(x))), integer(1))
size_bin <- ifelse(complex_size <= 4L, as.character(complex_size),
            ifelse(complex_size <= 7L, "5-7",
            ifelse(complex_size <= 10L, "8-10", "11+")))

counts <- as.data.frame(table(factor(size_bin, levels = SIZE_BIN_LEVELS)),
                        stringsAsFactors = FALSE)
names(counts) <- c("size_bin", "n_complexes")
counts$size_bin <- factor(counts$size_bin, levels = SIZE_BIN_LEVELS)
counts$label <- format(counts$n_complexes, big.mark = ",", trim = TRUE)

p <- ggplot(counts, aes(x = size_bin, y = n_complexes)) +
  geom_col(fill = "#4e79a7", width = 0.72) +
  geom_text(aes(label = label), vjust = -0.35, size = 3.6, color = "#1f2933") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Number of Complexes in Each Size Bin",
    x = "Complex size bin",
    y = "Number of complexes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#243b53", fill = NA, linewidth = 0.7),
    axis.title = element_text(face = "bold", color = "#1f2933"),
    axis.text = element_text(color = "#243b53"),
    plot.title = element_text(face = "bold", color = "#102a43")
  )

out_file <- file.path(OUT_DIR, "Figure4E_distribution.pdf")
ggsave(out_file, p, width = 6, height = 5)

cat("Saved:", out_file, "\n")
for (i in seq_len(nrow(counts))) {
  cat(sprintf("  size %-5s %5d complexes\n",
              as.character(counts$size_bin[i]), counts$n_complexes[i]))
}
