
rm(list = ls())

library(pROC)
library(PRROC)
library(dplyr)


BASE_DIR <- "/proj/c.zihao/work3/06singleCell/GSE96583"
OUT_DIR  <- "/proj/c.zihao/work3/06singleCell/GSE96583/04analysis/"

SAMPLES <- c("GSM2560248", "GSM2560249")

id_map <- c(
  GSM2560248 = "control", GSM2560249 = "stim"
)

METHOD_CFG <- list(
  list(name = "netComplex",        dir = "netcomplex/bulk",      file = "netComplex.csv",     fmt = "csv"),
  list(name = "netComplex_smooth", dir = "netcomplex/smooth",    file = "smooth.csv",         fmt = "csv"),
  list(name = "noRWR_bulk",        dir = "ablation/noRWR_bulk",        file = "noRWR_bulk.csv",        fmt = "csv"),
  list(name = "degRandPPI_bulk",   dir = "ablation/degRandPPI_bulk",   file = "degRandPPI_bulk.csv",   fmt = "csv"),
  list(name = "noRWR_smooth",      dir = "ablation/noRWR_smooth",      file = "noRWR_smooth.csv",      fmt = "csv"),
  list(name = "degRandPPI_smooth", dir = "ablation/degRandPPI_smooth", file = "degRandPPI_smooth.csv", fmt = "csv"),
  list(name = "netmean",           dir = "bench/netmean",        file = "netmean.csv",        fmt = "csv"),
  list(name = "perpas_adapted",    dir = "bench/perpas_adapted", file = "perpas_adapted.csv", fmt = "csv"),
  list(name = "gsva",              dir = "bench/gsva",           file = "gsva.csv",           fmt = "csv"),
  list(name = "ssgsea",            dir = "bench/ssgsea",         file = "ssgsea.csv",         fmt = "csv"),
  list(name = "aucell",            dir = "bench/aucell",         file = "aucell.csv",         fmt = "csv"),
  list(name = "ucell",             dir = "bench/ucell",          file = "ucell.csv",          fmt = "csv")
)


IFN_COMPLEXES <- c(
  "Comp_IFNAR1_IFNAR2_IFNW1_complex",
  "Comp_IFNAR1_IFNAR2_RACK1_STAT1_complex",
  "Comp_IFNAR1_IFNAR2_TYK2_complex",
  "Comp_IFNAR1_STAT1alpha_TYK2_complex",
  "Comp_IFNAR2_IFNA2_complex",

  "Comp_ISGF3_complex",
  "Comp_IRF1_STAT1_complex",
  "Comp_DTX3L_PARP9_STAT1_complex",

  "Comp_CBP_IRF3_complex",
  "Comp_DDX3_IKBKE_IRF3_complex",
  "Comp_HSP90_IRF3_MAVS_TBK1_TOMM70_complex",
  "Comp_HSP90_IRF3_TBK1_TOMM70_complex",
  "Comp_IKBKE_IRF3_TANK_TBK1_complex",
  "Comp_IRF7_MYD88_complex",
  "Comp_IRF7_MYD88_TRAF6_complex",
  "Comp_TANK_TBK1_IRF3_complex",

  "Comp_DDX58_MAVS_PPP6C_TRIM14_WHIP_complex",
  "Comp_DDX58_MAVS_TRIM31_complex",
  "Comp_DDX58_TRIM25_complex",
  "Comp_DDX58_TRIM25_YWHAE_complex",
  "Comp_PPP6C_TRIM14_WHIP_complex",
  "Comp_TRIM25_G3BP2_complex",

  "Comp_IFI16_BRCA1_complex",
  "Comp_IFI16_H2B_BRCA1_complex",
  "Comp_IFI16_inflammasome",
  "Comp_IFI1_IFI2_IFI3_complex",

  "Comp_DHX9_ADAR_vigilin_DNA_PK_Ku_antigen_complex"
)


read_sample_matrix <- function(path, fmt) {
  switch(fmt,
    tsv = as.matrix(read.table(path, sep = "\t", header = TRUE,
                               row.names = 1, check.names = FALSE)),
    csv = as.matrix(read.csv(path, row.names = 1, check.names = FALSE)),
    rds = as.matrix(readRDS(path))
  )
}

load_method_matrix <- function(cfg, samples, base_dir) {
  mat_list <- lapply(samples, function(sid) {
    path <- file.path(base_dir, cfg$dir, sid, cfg$file)
    read_sample_matrix(path, cfg$fmt)
  })
  common_rows <- Reduce(intersect, lapply(mat_list, rownames))
  mat_list <- lapply(mat_list, function(m) m[common_rows, , drop = FALSE])
  do.call(cbind, mat_list)
}

calc_auc <- function(group_bin, score_vec, positive = "stim") {
  tryCatch({
    roc_obj <- pROC::roc(
      response  = factor(group_bin, levels = c("control", positive)),
      predictor = score_vec,
      direction = "<",
      quiet     = TRUE
    )
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA_real_)
}

calc_prauc <- function(group_bin, score_vec, positive = "stim") {
  tryCatch({
    fg <- score_vec[group_bin == positive]
    bg <- score_vec[group_bin != positive]
    pr_obj <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE)
    as.numeric(pr_obj$auc.integral)
  }, error = function(e) NA_real_)
}


mat_list <- list()

for (cfg in METHOD_CFG) {
  cat("[", cfg$name, "] loading ...\n", sep = "")
  mat <- load_method_matrix(cfg, SAMPLES, BASE_DIR)
  cat("  matrix:", nrow(mat), "x", ncol(mat), "\n")
  mat_list[[cfg$name]] <- mat
}

keep <- intersect(IFN_COMPLEXES, Reduce(intersect, lapply(mat_list, rownames)))

is_all_na  <- function(m) vapply(keep, function(r) all(is.na(m[r, ])), logical(1))
drop_mask  <- Reduce(`|`, lapply(mat_list, is_all_na))

cat("\nIFN complexes listed:", length(IFN_COMPLEXES),
    "| present in all methods:", length(keep), "\n")
if (any(drop_mask)) {
  cat("Dropped", sum(drop_mask), "unevaluable (all-NA for >=1 method):\n")
  cat(paste0("  ", keep[drop_mask]), sep = "\n")
}
keep <- keep[!drop_mask]
cat("Evaluable complexes:", length(keep), "\n\n")
stopifnot(length(keep) > 0)

all_metric <- list()

for (cfg in METHOD_CFG) {
  sub_mat     <- mat_list[[cfg$name]][keep, , drop = FALSE]
  cell_ids    <- sub("_.*", "", colnames(sub_mat))
  cell_groups <- unname(id_map[cell_ids])
  data_t      <- t(sub_mat)

  auc_vec   <- apply(data_t, 2, function(s) calc_auc(cell_groups, s))
  prauc_vec <- apply(data_t, 2, function(s) calc_prauc(cell_groups, s))

  metric_df <- data.frame(
    method  = cfg$name,
    complex = rownames(sub_mat),
    auc     = as.numeric(auc_vec),
    pr_auc  = as.numeric(prauc_vec),
    stringsAsFactors = FALSE
  )

  all_metric[[cfg$name]] <- metric_df

  cat(sprintf("[%s] mean AUROC=%.4f  mean PR-AUC=%.4f\n", cfg$name,
              mean(metric_df$auc,    na.rm = TRUE),
              mean(metric_df$pr_auc, na.rm = TRUE)))
}


metric_data <- do.call(rbind, all_metric)
metric_data$auc    <- round(metric_data$auc,    3)
metric_data$pr_auc <- round(metric_data$pr_auc, 3)
write.csv(metric_data,
          file.path(OUT_DIR, "auc.csv"),
          row.names = FALSE)

summary_df <- metric_data %>%
  group_by(method) %>%
  summarise(
    n_complexes = n(),
    mean_auc    = round(mean(auc,    na.rm = TRUE), 3),
    mean_pr_auc = round(mean(pr_auc, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auc))

write.csv(summary_df,
          file.path(OUT_DIR, "summary.csv"),
          row.names = FALSE)

cat("\n-- IFN complexes: mean AUROC ranking --\n")
print(as.data.frame(summary_df))

cat("\n-- Per-complex mean AUC (across methods) --\n")
per_complex <- metric_data %>%
  group_by(complex) %>%
  summarise(
    mean_auc    = round(mean(auc,    na.rm = TRUE), 3),
    mean_pr_auc = round(mean(pr_auc, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_auc))
print(as.data.frame(per_complex))

cat("\nResults saved to:", OUT_DIR, "\n")
