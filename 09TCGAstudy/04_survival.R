library(data.table)
library(survival)
library(ggplot2)

ROOT <- "/proj/c.zihao/work3"
STUDY_DIR <- file.path(ROOT, "09TCGAstudy")
RESULT_DIR <- file.path(STUDY_DIR, "results")
PLOT_DIR <- file.path(STUDY_DIR, "plots")
SURVIVAL_DIR <- file.path(ROOT, "03Survival", "KIRC", "TCGAKIRC")
CLINICAL_FILE <- file.path(SURVIVAL_DIR, "data", "pd.csv")
NET_FILE <- file.path(SURVIVAL_DIR, "ablation", "netcomplex", "complex_score",
                      "TCGAKIRC_netcomplex_complex_score.csv")
MEAN_FILE <- file.path(SURVIVAL_DIR, "bench", "mean", "mean.csv")
EXPR_FILE <- file.path(SURVIVAL_DIR, "data", "exprSet_filtered.csv")
TN_EXPR_FILE <- file.path(ROOT, "00data", "TCGAnt", "exprset", "KIRC_exprSet_filtered.csv")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

cox_row <- function(data, method) {
  data <- data[is.finite(score) & !is.na(OS) & !is.na(OS_Time) & OS_Time > 0]

  fit <- coxph(Surv(OS_Time, OS) ~ scale(score), data = data, x = TRUE)
  s <- summary(fit)$coefficients[1L, ]
  ci <- summary(fit)$conf.int[1L, ]

  staged <- data[!is.na(Stage)]
  fit_staged <- coxph(Surv(OS_Time, OS) ~ scale(score) + strata(Stage), data = staged, x = TRUE)
  ss <- summary(fit_staged)$coefficients[1L, ]
  cis <- summary(fit_staged)$conf.int[1L, ]
  stage_stratified_score_ph_p <- cox.zph(fit_staged)$table[1L, "p"]

  list(
    row = data.table(
      method = method, n = nrow(data), n_events = sum(data$OS),
      cox_hr_per_sd = unname(ci["exp(coef)"]),
      cox_ci_low = unname(ci["lower .95"]), cox_ci_high = unname(ci["upper .95"]),
      cox_p_value = unname(s["Pr(>|z|)"]),
      stage_stratified_n = nrow(staged), stage_stratified_events = sum(staged$OS),
      stage_stratified_hr_per_sd = unname(cis["exp(coef)"]),
      stage_stratified_ci_low = unname(cis["lower .95"]),
      stage_stratified_ci_high = unname(cis["upper .95"]),
      stage_stratified_p_value = unname(ss["Pr(>|z|)"]),
      stage_stratified_score_ph_p = stage_stratified_score_ph_p
    ),
    data = data
  )
}

time_varying_curve <- function(data, method, years = seq(0.5, 10, by = 0.1)) {
  data <- copy(data)[!is.na(Stage)]
  data[, z_score := as.numeric(scale(score))]
  fit <- coxph(
    Surv(OS_Time, OS) ~ z_score + tt(z_score) + strata(Stage), data = data,
    tt = function(x, t, ...) x * log(pmax(t, 0.1))
  )
  beta <- coef(fit)
  variance <- vcov(fit)
  log_time <- log(years)
  linear_predictor <- beta[1L] + beta[2L] * log_time
  se <- sqrt(variance[1L, 1L] + log_time^2 * variance[2L, 2L] +
             2 * log_time * variance[1L, 2L])
  interaction_p <- summary(fit)$coefficients[2L, "Pr(>|z|)"]
  data.table(
    method = method, years = years,
    hr_per_sd = exp(linear_predictor),
    ci_low = exp(linear_predictor - 1.96 * se),
    ci_high = exp(linear_predictor + 1.96 * se),
    time_interaction_p = interaction_p
  )
}

selected <- fread(file.path(RESULT_DIR, "kirc_selected_case.csv"))
members <- fread(file.path(RESULT_DIR, "kirc_selected_member_genes.csv"))

clinical <- fread(CLINICAL_FILE)
clinical[, OS := as.numeric(OS)]
clinical[, OS_Time := as.numeric(OS_Time)]
clinical[, Stage := factor(Stage, levels = c("I", "II", "III", "IV"))]

files <- c(
  net = NET_FILE,
  mean_score = MEAN_FILE,
  expr = EXPR_FILE,
  tn_expr = TN_EXPR_FILE
)
matrices <- vector("list", length(files))
names(matrices) <- names(files)

for (name in names(files)) {
  x <- fread(files[[name]], check.names = FALSE)
  matrices[[name]] <- as.matrix(x[, -1, with = FALSE])
  storage.mode(matrices[[name]]) <- "numeric"
  rownames(matrices[[name]]) <- x[[1L]]
}
net <- matrices$net
mean_score <- matrices$mean_score
expr <- matrices$expr
tn_expr <- matrices$tn_expr

sample_ids <- colnames(tn_expr)
patient <- sub("_(01A|11A)$", "", sample_ids)
tumour <- grepl("_01A$", sample_ids)
normal <- grepl("_11A$", sample_ids)
patient_id <- sort(intersect(patient[tumour], patient[normal]))
tn_pairs <- list(
  patient_id = patient_id,
  tumour_id = unname(setNames(sample_ids[tumour], patient[tumour])[patient_id]),
  normal_id = unname(setNames(sample_ids[normal], patient[normal])[patient_id])
)

gene_long <- rbindlist(lapply(members$gene, function(gene) {
  data.table(
    patient = rep(tn_pairs$patient_id, each = 2L),
    condition = factor(rep(c("Normal", "Tumour"), length(tn_pairs$patient_id)),
                       levels = c("Normal", "Tumour")),
    value = as.numeric(c(rbind(tn_expr[gene, tn_pairs$normal_id],
                               tn_expr[gene, tn_pairs$tumour_id]))),
    panel = gene
  )
}))
gene_labels <- merge(
  members[, .(panel = gene, label = sprintf("d[z] = %.2f\nFDR = %.2g", dz, fdr))],
  gene_long[, .(value = max(value) + 0.035 * diff(range(value))), by = panel], by = "panel", sort = FALSE
)
gene_labels[, condition := factor("Normal", levels = c("Normal", "Tumour"))]
fig_member <- ggplot(gene_long, aes(condition, value, group = patient)) +
  geom_line(colour = "grey55", alpha = 0.28, linewidth = 0.35) +
  geom_point(aes(colour = condition), alpha = 0.62, size = 1.35) +
  geom_boxplot(aes(group = condition), width = 0.48, outlier.shape = NA,
               fill = NA, colour = "black", linewidth = 0.45) +
  geom_text(data = gene_labels, aes(x = condition, y = value, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1.15, size = 3.1) +
  facet_wrap(~panel, scales = "free_y") +
  scale_colour_manual(values = c(Normal = "#277DA1", Tumour = "#D1495B"), guide = "none") +
  labs(title = "KIRC case-member expression",
       x = NULL, y = "Expression") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

net_data <- merge(
  clinical,
  data.table(Sample = colnames(net), score = as.numeric(net[selected$complex, ])),
  by = "Sample"
)
net_result <- cox_row(net_data, "NetComplex")

mean_data <- merge(
  clinical,
  data.table(Sample = colnames(mean_score), score = as.numeric(mean_score[selected$complex, ])),
  by = "Sample"
)
mean_result <- cox_row(mean_data, "Member mean")

member_results <- setNames(lapply(members$gene, function(gene) {
  gene_data <- merge(
    clinical,
    data.table(Sample = colnames(expr), score = as.numeric(expr[gene, ])),
    by = "Sample"
  )
  cox_row(gene_data, gene)
}), members$gene)

survival_table <- rbindlist(c(
  list(net_result$row, mean_result$row),
  lapply(member_results, `[[`, "row")
))

time_varying <- rbindlist(list(
  time_varying_curve(net_result$data, "NetComplex"),
  time_varying_curve(mean_result$data, "Member mean"),
  time_varying_curve(member_results[["WRN"]]$data, "WRN"),
  time_varying_curve(member_results[["XRCC5"]]$data, "XRCC5"),
  time_varying_curve(member_results[["XRCC6"]]$data, "XRCC6")
))

output_tables <- list(
  kirc_survival_associations.csv = survival_table,
  kirc_time_varying_cox.csv = time_varying
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

  fwrite(table_to_write, file.path(RESULT_DIR, file_name))
}

cox_plot <- copy(survival_table)
cox_plot[, display_method := factor(method, levels = rev(c("NetComplex", "Member mean", "WRN", "XRCC5", "XRCC6")))]
cox_plot[, ph_note := ifelse(stage_stratified_score_ph_p < 0.05, "PH assumption violated", "")]
fig_e <- ggplot(cox_plot, aes(display_method, stage_stratified_hr_per_sd, colour = method)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(ymin = stage_stratified_ci_low, ymax = stage_stratified_ci_high), width = 0.12, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(y = stage_stratified_ci_high * 1.08, label = ph_note),
            colour = "#B2182B", hjust = 0, size = 3.1, show.legend = FALSE) +
  coord_flip(clip = "off") +
  scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2), limits = c(0.5, 2)) +
  scale_colour_manual(values = c("NetComplex" = "#1B9E77", "Member mean" = "#5F5F5F",
                                 "WRN" = "#7570B3", "XRCC5" = "#377EB8", "XRCC6" = "#E7298A"),
                      guide = "none") +
  labs(
    title = "KIRC survival association",
    x = NULL, y = "Hazard ratio for overall survival (log scale)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), plot.margin = margin(5.5, 110, 5.5, 5.5))
feature_colours <- c("NetComplex" = "#1B9E77", "Member mean" = "#5F5F5F",
                     "WRN" = "#7570B3", "XRCC5" = "#377EB8", "XRCC6" = "#E7298A")
feature_order <- c("NetComplex", "Member mean", "WRN", "XRCC5", "XRCC6")
ph_violators <- cox_plot[stage_stratified_score_ph_p < 0.05, method]
time_varying[, method := factor(method, levels = feature_order)]
time_ribbons <- time_varying[method %in% ph_violators]
fig_f <- ggplot(time_varying, aes(years, hr_per_sd, colour = method, fill = method)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_ribbon(data = time_ribbons, aes(ymin = ci_low, ymax = ci_high), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = feature_colours, name = NULL) +
  scale_fill_manual(values = feature_colours, guide = "none") +
  scale_y_log10(breaks = c(0.4, 0.5, 0.75, 1, 1.5, 2), limits = c(0.35, 2)) +
  labs(
    title = "Time-varying survival associations",
    x = "Overall survival time (years)", y = "Hazard ratio per 1 SD (log scale)",
    colour = NULL, fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = c(0.83, 0.84))

grDevices::pdf(file.path(PLOT_DIR, "SupplementaryFigure11.pdf"), width = 17, height = 5.5, useDingbats = FALSE)
grid::grid.newpage()
print(fig_member, vp = grid::viewport(x = 0.19, y = 0.47, width = 0.98 * 0.38, height = 0.91))
print(fig_e, vp = grid::viewport(x = 0.535, y = 0.47, width = 0.98 * 0.31, height = 0.91))
print(fig_f, vp = grid::viewport(x = 0.845, y = 0.47, width = 0.98 * 0.31, height = 0.91))
grid::grid.text("A", x = 0.01, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grid::grid.text("B", x = 0.39, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grid::grid.text("C", x = 0.70, y = 0.99, just = c("left", "top"),
                gp = grid::gpar(fontface = "bold", fontsize = 16))
grDevices::dev.off()

net_row <- net_result$row
cat(sprintf("NetComplex survival: n=%d, events=%d, HR/SD=%.3f, P=%.3g; stage-stratified HR=%.3f, P=%.3g\n",
            net_row$n, net_row$n_events, net_row$cox_hr_per_sd, net_row$cox_p_value,
            net_row$stage_stratified_hr_per_sd, net_row$stage_stratified_p_value))
