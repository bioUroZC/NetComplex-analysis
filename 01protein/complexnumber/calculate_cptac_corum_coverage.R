

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_arg)))
} else {
  normalizePath(getwd())
}
root_dir <- normalizePath(file.path(script_dir, "..", ".."))

corum_file <- file.path(root_dir, "00data", "corum", "corum_humanComplexes.txt")
expr_dir <- file.path(root_dir, "00data", "cptacT", "exprset")
summary_file <- file.path(script_dir, "cptac_corum_coverage_summary.csv")
detail_file <- file.path(script_dir, "cptac_corum_coverage_detail.csv")

stopifnot(file.exists(corum_file), dir.exists(expr_dir))
expr_files <- sort(list.files(expr_dir, pattern = "_exprSet_filtered\\.csv$",
                              full.names = TRUE))
if (length(expr_files) == 0L) stop("No CPTAC expression matrices found in: ", expr_dir)

raw_corum <- fread(corum_file, sep = "\t", quote = "", fill = Inf,
                   header = TRUE, data.table = FALSE)
required_cols <- c("complex_name", "subunits_gene_name")
if (!all(required_cols %in% colnames(raw_corum))) {
  stop("CORUM file does not contain required columns: ",
       paste(required_cols, collapse = ", "))
}

complexes <- raw_corum[
  !is.na(raw_corum$complex_name) & raw_corum$complex_name != "" &
    !is.na(raw_corum$subunits_gene_name) & raw_corum$subunits_gene_name != "",
  required_cols,
  drop = FALSE
]

gene_set_key <- vapply(strsplit(complexes$subunits_gene_name, ";", fixed = TRUE),
                       function(genes) {
                         paste(sort(unique(trimws(genes))), collapse = ";")
                       },
                       character(1))
complexes <- complexes[!duplicated(gene_set_key), , drop = FALSE]
gene_set_key <- gene_set_key[!duplicated(gene_set_key)]

has_overlap <- function(subunit_strings, panel_genes) {
  vapply(strsplit(subunit_strings, ";", fixed = TRUE),
         function(genes) any(trimws(genes) %in% panel_genes),
         logical(1))
}

panel_list <- lapply(expr_files, function(path) fread(path, select = 1)[[1]])
cohorts <- sub("_exprSet_filtered\\.csv$", "", basename(expr_files))
reference_panel <- panel_list[[1L]]

summary <- rbindlist(lapply(seq_along(expr_files), function(i) {
  retained <- has_overlap(complexes$subunits_gene_name, panel_list[[i]])
  data.table(
    cohort = cohorts[[i]],
    expression_genes = length(panel_list[[i]]),
    panel_identical_to_first = setequal(panel_list[[i]], reference_panel),
    corum_complexes_deduplicated = nrow(complexes),
    complexes_retained = sum(retained),
    complexes_excluded = sum(!retained),
    retained_percent = round(100 * mean(retained), 1),
    excluded_percent = round(100 * mean(!retained), 1)
  )
}))
setorder(summary, cohort)

union_panel <- unique(unlist(panel_list, use.names = FALSE))
detail <- data.table(
  complex_name = complexes$complex_name,
  subunit_genes = complexes$subunits_gene_name,
  gene_set_key = gene_set_key,
  retained_in_any_cptac_panel = has_overlap(complexes$subunits_gene_name, union_panel)
)

fwrite(summary, summary_file)
fwrite(detail, detail_file)

message("Deduplicated CORUM complexes: ", nrow(complexes))
message("Wrote: ", summary_file)
message("Wrote: ", detail_file)
