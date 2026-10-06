library(data.table)
library(limma)

ROOT <- "/proj/c.zihao/work3"
STUDY_DIR <- file.path(ROOT, "09TCGAstudy")
RESULT_DIR <- file.path(STUDY_DIR, "results")
TCGA_SCORE_ROOT <- file.path(ROOT, "02TCGAnt", "02stringTest")
EXPR_DIR <- file.path(ROOT, "00data", "TCGAnt", "exprset")
CORUM_FILE <- file.path(ROOT, "00data", "corum", "complexData.csv")

limma_paired <- function(delta) {
  valid <- rowSums(is.finite(delta)) == ncol(delta) &
    apply(delta, 1L, sd) > 0

  valid_delta <- delta[valid, , drop = FALSE]
  fit <- eBayes(lmFit(valid_delta, design = matrix(1, ncol(valid_delta), 1L)))

  out <- data.table(
    complex = rownames(delta),
    mean_delta = NA_real_,
    dz = NA_real_,
    p_value = NA_real_
  )
  out[valid, mean_delta := rowMeans(valid_delta)]
  out[valid, dz := rowMeans(valid_delta) / apply(valid_delta, 1L, sd)]
  out[valid, p_value := fit$p.value[, 1L]]
  out[, fdr := p.adjust(p_value, "BH")]
}

CANCER <- "KIRC"
OUT_DIR <- RESULT_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NET_FDR <- 0.01
NET_DZ <- 0.8
MEAN_FDR <- 0.05
MEAN_ABS_DZ <- 0.3
MEMBER_MEDIAN_ABS_DZ <- 0.2
MEMBER_MAX_SIGNIFICANT_FRACTION <- 1 / 3
MIN_ANNOTATED_MEMBERS <- 3L
MIN_MEASURED_MEMBERS <- 2L
MIN_COVERAGE <- 0.5

net_raw <- fread(
  file.path(TCGA_SCORE_ROOT, "netcomplex", "complex_score", CANCER,
            sprintf("%s_netcomplex_complex_score.csv", CANCER)),
  check.names = FALSE
)
net_scores <- as.matrix(net_raw[, -1, with = FALSE])
storage.mode(net_scores) <- "numeric"
rownames(net_scores) <- net_raw[[1L]]

sample_ids <- colnames(net_scores)
patient <- sub("_(01A|11A)$", "", sample_ids)
tumour <- grepl("_01A$", sample_ids)
normal <- grepl("_11A$", sample_ids)
patient_id <- sort(intersect(patient[tumour], patient[normal]))

pairs <- data.table(
  patient_id,
  tumour_id = unname(setNames(sample_ids[tumour], patient[tumour])[patient_id]),
  normal_id = unname(setNames(sample_ids[normal], patient[normal])[patient_id])
)

net_delta <- net_scores[, pairs$tumour_id] - net_scores[, pairs$normal_id]
colnames(net_delta) <- pairs$patient_id

mean_raw <- fread(
  file.path(TCGA_SCORE_ROOT, "bench", "mean", CANCER, "mean.csv"),
  check.names = FALSE
)
mean_scores <- as.matrix(mean_raw[, -1, with = FALSE])
storage.mode(mean_scores) <- "numeric"
rownames(mean_scores) <- mean_raw[[1L]]

mean_delta <- mean_scores[, pairs$tumour_id] - mean_scores[, pairs$normal_id]
colnames(mean_delta) <- pairs$patient_id

expr_raw <- fread(
  file.path(EXPR_DIR, sprintf("%s_exprSet_filtered.csv", CANCER)),
  check.names = FALSE
)
expr <- as.matrix(expr_raw[, -1, with = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- expr_raw[[1L]]

net_stats <- limma_paired(net_delta)
setnames(net_stats, c("mean_delta", "dz", "p_value", "fdr"),
         c("net_delta", "net_dz", "net_p_value", "net_fdr"))
mean_stats <- limma_paired(mean_delta)
setnames(mean_stats, c("mean_delta", "dz", "p_value", "fdr"),
         c("mean_delta", "mean_dz", "mean_p_value", "mean_fdr"))
tab <- merge(net_stats, mean_stats, by = "complex")
corum <- fread(CORUM_FILE)
member_list <- setNames(lapply(corum$Genes, function(x) unique(strsplit(x, ";", fixed = TRUE)[[1L]])), corum$Complex)
coverage <- data.table(
  complex = names(member_list),
  n_annotated = lengths(member_list),
  n_measured = vapply(member_list, function(x) sum(x %in% rownames(expr)), integer(1))
)
coverage[, coverage := n_measured / n_annotated]
tab <- merge(tab, coverage, by = "complex", all.x = TRUE)
tab[, `:=`(cancer = CANCER, n_pairs = nrow(pairs),
           net_tumour_mean = rowMeans(net_scores[complex, pairs$tumour_id, drop = FALSE]),
           net_normal_mean = rowMeans(net_scores[complex, pairs$normal_id, drop = FALSE]))]

expr_delta <- expr[, pairs$tumour_id] - expr[, pairs$normal_id]
colnames(expr_delta) <- pairs$patient_id
gene_stats <- limma_paired(expr_delta)
setnames(gene_stats, "complex", "gene")
members <- rbindlist(lapply(names(member_list), function(cx) data.table(complex = cx, gene = member_list[[cx]])))
members <- merge(members, gene_stats, by = "gene", all.x = TRUE)
members[, significant := !is.na(fdr) & fdr < MEAN_FDR]
member_summary <- members[, .(
  n_member_genes = .N,
  n_member_gene_significant = sum(significant),
  member_gene_significant_fraction = mean(significant),
  member_gene_median_abs_dz = median(abs(dz), na.rm = TRUE)
), by = complex]
tab <- merge(tab, member_summary, by = "complex", all.x = TRUE)

candidates <- tab[
  net_fdr < NET_FDR & net_dz >= NET_DZ &
  mean_fdr > MEAN_FDR & abs(mean_dz) < MEAN_ABS_DZ &
  member_gene_median_abs_dz < MEMBER_MEDIAN_ABS_DZ &
  member_gene_significant_fraction <= MEMBER_MAX_SIGNIFICANT_FRACTION &
  n_annotated >= MIN_ANNOTATED_MEMBERS & n_measured >= MIN_MEASURED_MEMBERS & coverage >= MIN_COVERAGE
]
setorder(candidates, -net_dz, net_fdr, complex)

selected <- candidates[1L]

output_tables <- list(
  kirc_all_complexes.csv = tab,
  kirc_candidates.csv = candidates,
  kirc_selected_case.csv = selected,
  kirc_selected_member_genes.csv = members[complex == selected$complex]
)

for (file_name in names(output_tables)) {
  table_to_write <- copy(output_tables[[file_name]])
  numeric_columns <- names(table_to_write)[vapply(table_to_write, is.numeric, logical(1))]

  for (column in numeric_columns) {
    values <- table_to_write[[column]]
    tiny_nonzero <- is.finite(values) & values != 0 & abs(values) < 1e-5
    whole_number <- is.finite(values) & values == floor(values)

    table_to_write[[column]] <- ifelse(
      is.na(values),
      NA_character_,
      ifelse(tiny_nonzero,
             formatC(values, format = "e", digits = 5),
             ifelse(whole_number,
                    formatC(values, format = "f", digits = 0),
                    formatC(values, format = "f", digits = 5)))
    )
  }

  fwrite(table_to_write, file.path(OUT_DIR, file_name))
}

cat(sprintf("Selected KIRC case: %s\n", selected$complex))
cat(sprintf("NetComplex: dz=%.3f, FDR=%.3g; Mean: dz=%.3f, FDR=%.3g\n",
            selected$net_dz, selected$net_fdr, selected$mean_dz, selected$mean_fdr))
cat(sprintf("Member genes: median |dz|=%.3f; significant fraction=%.2f\n",
            selected$member_gene_median_abs_dz, selected$member_gene_significant_fraction))
