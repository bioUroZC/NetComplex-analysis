library(data.table)

RES_DIR <- "/proj/c.zihao/work3/10Cellstudy/results"
COMPARATORS <- c("ssgsea", "aucell", "ucell", "gsva", "netmean")

NC_MIN <- 0.80
COMP_MAX <- 0.60
MARGIN <- 0.15
PRIMARY_DATASET <- "GSE96583"
EXTERNAL_DATASET <- "GSE226572"
N_CASES <- 3L

wide <- fread(file.path(RES_DIR, "percomplex_auc_wide.csv"))
pos <- wide[set == "positive"]

comp_mat <- as.matrix(pos[, ..COMPARATORS])
comp_mat[is.na(comp_mat)] <- -Inf
pos[, best_comparator := apply(comp_mat, 1, max)]
pos[!is.finite(best_comparator), best_comparator := NA_real_]
pos[, best_comparator_name := COMPARATORS[max.col(comp_mat, ties.method = "first")]]
pos[is.na(best_comparator), best_comparator_name := NA_character_]
pos[, margin_over_best := round(netComplex - best_comparator, 3)]

sel <- pos[is.finite(netComplex) & netComplex > NC_MIN &
           is.finite(best_comparator) & best_comparator < COMP_MAX]
setorder(sel, -margin_over_best)
fwrite(sel, file.path(RES_DIR, "nc_detected_only.csv"))

both <- sel[, .N, by = complex][N == length(unique(pos$dataset)), complex]
sel_both <- sel[complex %in% both]
setorder(sel_both, complex, dataset)
fwrite(sel_both, file.path(RES_DIR, "nc_detected_both.csv"))

primary <- sel[dataset == PRIMARY_DATASET]
external <- pos[dataset == EXTERNAL_DATASET,
  .(complex, external_netComplex = netComplex,
    external_best_comparator = best_comparator,
    external_best_comparator_name = best_comparator_name,
    external_margin_over_best = margin_over_best)]
cases <- merge(primary, external, by = "complex", all.x = TRUE)
cases[, external_support_directional := is.finite(external_margin_over_best) &
        external_margin_over_best > 0]
setorder(cases, -external_margin_over_best, -margin_over_best, complex)
cases <- head(cases, N_CASES)
cases[, `:=`(primary_dataset = PRIMARY_DATASET,
             external_dataset = EXTERNAL_DATASET,
             strict_hit_replicated = FALSE)]
fwrite(cases, file.path(RES_DIR, "case_complexes.csv"))

counter <- pos[is.finite(best_comparator) & best_comparator - netComplex > MARGIN]
setorder(counter, margin_over_best)
fwrite(counter, file.path(RES_DIR, "counter_examples.csv"))

summ <- pos[, .(positive_complexes = .N,
                netComplex_mean = round(mean(netComplex), 3),
                best_comparator_mean = round(mean(best_comparator, na.rm = TRUE), 3),
                selected = sum(netComplex > NC_MIN & best_comparator < COMP_MAX, na.rm = TRUE),
                counter_examples = sum(best_comparator - netComplex > MARGIN, na.rm = TRUE)),
            by = dataset]
fwrite(summ, file.path(RES_DIR, "selection_summary.csv"))

cat(sprintf("Thresholds: netComplex > %.2f, best of {%s} < %.2f\n\n",
            NC_MIN, paste(COMPARATORS, collapse = ", "), COMP_MAX))
print(summ, row.names = FALSE)
cat("\n=== selected complexes ===\n")
print(sel[, .(dataset, complex, n_members, netComplex = round(netComplex, 3),
              best_comparator = round(best_comparator, 3), best_comparator_name,
              margin_over_best, genes)], row.names = FALSE)
cat(sprintf("\nSelected in both datasets: %d complexes\n", uniqueN(sel_both$complex)))
print(unique(sel_both$complex))
cat(sprintf("\nCounter-examples (a comparator ahead by > %.2f): %d\n", MARGIN, nrow(counter)))
print(counter[, .(dataset, complex, netComplex = round(netComplex, 3),
                  best_comparator = round(best_comparator, 3), best_comparator_name)],
      row.names = FALSE)
cat("\n=== mechanism case complexes (strict primary hits; external support is descriptive) ===\n")
print(cases[, .(complex, primary_dataset, netComplex,
                best_comparator, margin_over_best, external_dataset, external_netComplex,
                external_best_comparator, external_margin_over_best,
                external_support_directional)], row.names = FALSE)
