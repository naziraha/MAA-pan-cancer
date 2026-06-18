# ============================================================
# RQ2 — CLINICAL, GENOMIC & STAGE-TRANSITION ANALYSIS
# Questions:
#   (a) Which MAA genes predict survival independently?
#   (b) How does MAA change across Stage I→IV (cancer progression)?
#   (c) Which MAA genes are recurrently mutated / amplified in cancer?
# Outputs: Survival forest, KM plots, stage-trend plots, mutation table
# ============================================================

source("R startup/00_gene_panel.R")
library(survival); library(survminer); library(ggplot2); library(dplyr)
library(openxlsx); library(ggrepel); library(tidyr)

## ---------- 1. Load ----------
master     <- readRDS("xena_data/master.rds")
tub_tumor  <- readRDS("xena_data/tub_mat_tumor.rds")
gistic     <- readRDS("xena_data/gistic_thresh.rds")
mc3        <- readRDS("xena_data/mc3.rds")
found_genes <- rownames(tub_tumor)

## ================================================================
## PART A: SURVIVAL ANALYSIS
## ================================================================

## ---------- 2. Determine available covariates ----------
# stage_num may be all-NA if clinical file lacked stage; make Cox models adaptive
has_stage  <- "stage_num" %in% names(master) && sum(!is.na(master$stage_num)) > 100
has_age    <- "age_at_initial_pathologic_diagnosis" %in% names(master) &&
              sum(!is.na(master$age_at_initial_pathologic_diagnosis)) > 100
has_gender <- "gender" %in% names(master) && sum(!is.na(master$gender)) > 100
cat("Covariates available — stage:", has_stage,
    "| age:", has_age, "| gender:", has_gender, "\n")

## ---------- 3. Univariate Cox — Overall Survival (OS) ----------
cat("Running univariate Cox regressions...\n")
surv_results <- lapply(found_genes, function(g) {
  # Only require OS, OS.time, and gene expression (no covariate NAs blocking analysis)
  df <- master[, c("OS", "OS.time", g), drop = FALSE]
  df <- df[complete.cases(df) & df[["OS.time"]] > 0, ]
  if (nrow(df) < 50) return(NULL)
  form <- as.formula(paste0("Surv(OS.time, OS) ~ `", g, "`"))
  fit  <- tryCatch(coxph(form, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coef
  data.frame(gene = g, n = nrow(df),
             HR = exp(s[1,1]), HR_lo = exp(s[1,1] - 1.96*s[1,3]),
             HR_hi = exp(s[1,1] + 1.96*s[1,3]),
             p = s[1,5], stringsAsFactors = FALSE)
})
surv_os <- do.call(rbind, Filter(Negate(is.null), surv_results))
if (is.null(surv_os) || nrow(surv_os) == 0)
  stop("No univariate Cox results — check OS/OS.time columns in master table")
surv_os$q <- p.adjust(surv_os$p, method = "BH")
surv_os   <- left_join(surv_os, gene_layer_table, by = "gene") %>% arrange(p)
cat("Univariate Cox done:", nrow(surv_os), "genes | significant (q<0.05):",
    sum(surv_os$q < 0.05, na.rm=TRUE), "\n")

## ---------- 4. Multivariate Cox — adjusted for available covariates ----------
cat("Running multivariate Cox regressions...\n")
# Build covariate list dynamically based on what's available
cov_vars <- c(
  if (has_age)    "age_at_initial_pathologic_diagnosis",
  if (has_stage)  "stage_num",
  if (has_gender) "gender"
)
multiv_results <- lapply(found_genes, function(g) {
  keep_cols <- c("OS", "OS.time", cov_vars, g)
  df <- master[, keep_cols, drop = FALSE]
  names(df)[ncol(df)] <- "gene_expr"
  df <- df[complete.cases(df) & df[["OS.time"]] > 0, ]
  if (nrow(df) < 50) return(NULL)
  rhs <- paste(c("gene_expr", cov_vars), collapse = " + ")
  form <- as.formula(paste("Surv(OS.time, OS) ~", rhs))
  fit  <- tryCatch(coxph(form, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coef
  data.frame(gene = g, n = nrow(df),
             HR = exp(s[1,1]), HR_lo = exp(s[1,1] - 1.96*s[1,3]),
             HR_hi = exp(s[1,1] + 1.96*s[1,3]),
             p = s[1,5], stringsAsFactors = FALSE)
})
multiv_os <- do.call(rbind, Filter(Negate(is.null), multiv_results))
if (is.null(multiv_os) || nrow(multiv_os) == 0) {
  cat("Note: multivariate Cox returned no results — using univariate results instead\n")
  multiv_os <- surv_os
} else {
  multiv_os$q <- p.adjust(multiv_os$p, method = "BH")
  multiv_os   <- left_join(multiv_os, gene_layer_table, by = "gene") %>% arrange(p)
}

n_sig_mv <- sum(multiv_os$q < 0.05, na.rm = TRUE)
cat("Genes independently prognostic (multivariate q<0.05):", n_sig_mv, "/", nrow(multiv_os), "\n")

## ---------- 4. Forest plot — top 30 multivariate genes ----------
top30 <- head(multiv_os[multiv_os$q < 0.05, ], 30)
if (nrow(top30) > 0) {
  top30$gene <- factor(top30$gene, levels = rev(top30$gene))
  p_forest <- ggplot(top30, aes(x = HR, xmin = HR_lo, xmax = HR_hi,
                                 y = gene, colour = layer)) +
    geom_pointrange(size = 0.7) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_colour_manual(values = layer_colors, na.value = "grey60") +
    scale_x_log10() +
    labs(title = "RQ2: Multivariate Cox — Top Prognostic MAA Genes (OS)",
         x = "Hazard Ratio (95% CI)", y = NULL) +
    theme_bw(base_size = 12)
  ggsave("figures/Fig3_survival_forest_OS.pdf", p_forest, width = 10, height = 10)
}

## ---------- 5. KM plots — top 10 genes ----------
km_genes <- head(multiv_os$gene[multiv_os$q < 0.05], 10)
for (g in km_genes) {
  df <- master[complete.cases(master[, c("OS", "OS.time", g)]) &
                 master$OS.time > 0, ]
  if (nrow(df) < 50) next
  df$grp <- ifelse(df[[g]] >= median(df[[g]], na.rm = TRUE), "High", "Low")
  fit <- survfit(Surv(OS.time / 365.25, OS) ~ grp, data = df)
  p_km <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                     title = paste("OS by", g, "expression (median split)"),
                     xlab = "Time (years)", palette = c("#E41A1C","#377EB8"))
  pdf(paste0("figures/KM_OS_", g, ".pdf"), width = 8, height = 7)
  print(p_km)
  dev.off()
}
cat("KM plots saved for top genes.\n")

## ================================================================
## PART B: STAGE-TRANSITION ANALYSIS (Progression model)
## ================================================================

## ---------- 6. Stage I→IV expression trends ----------
cat("Analysing stage-expression trends...\n")
stage_df <- master[!is.na(master$stage_num) & master$stage_num %in% 1:4, ]

stage_trend <- lapply(found_genes, function(g) {
  df <- stage_df[!is.na(stage_df[[g]]), c("stage_num", g)]
  if (nrow(df) < 50) return(NULL)
  fit <- lm(as.formula(paste0("`", g, "` ~ stage_num")), data = df)
  s   <- summary(fit)$coef
  data.frame(gene = g, beta = s[2,1], se = s[2,2], t = s[2,3],
             p = s[2,4], stringsAsFactors = FALSE)
})
stage_trend <- do.call(rbind, Filter(Negate(is.null), stage_trend))
stage_trend$q <- p.adjust(stage_trend$p, method = "BH")
stage_trend <- left_join(stage_trend, gene_layer_table, by = "gene") %>%
  arrange(p)

# Genes with significant stage-linear trend (cancer progression markers)
stage_up   <- filter(stage_trend, q < 0.05, beta > 0) %>% arrange(desc(beta))
stage_down <- filter(stage_trend, q < 0.05, beta < 0) %>% arrange(beta)
cat("Genes increasing with stage:", nrow(stage_up), "\n")
cat("Genes decreasing with stage:", nrow(stage_down), "\n")

## ---------- 7. Stage boxplot — top 6 stage-increasing genes ----------
top_stage_genes <- head(stage_up$gene, 6)
if (length(top_stage_genes) > 0) {
  plot_df <- stage_df %>%
    select(stage_num, all_of(top_stage_genes)) %>%
    pivot_longer(-stage_num, names_to = "gene", values_to = "expr") %>%
    mutate(stage_num = factor(stage_num, levels = 1:4,
                               labels = c("I","II","III","IV")))
  p_stage <- ggplot(plot_df, aes(x = stage_num, y = expr, fill = stage_num)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("I"="#FFFFCC","II"="#FEB24C","III"="#FC4E2A","IV"="#800026")) +
    labs(title = "RQ2: MAA Expression Across Tumour Stages",
         x = "Pathological Stage", y = "log2 TPM") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave("figures/RQ2_stage_trends.pdf", p_stage, width = 10, height = 8)
}

## ================================================================
## PART C: GENOMIC INTEGRITY — MUTATION & COPY NUMBER
## ================================================================

## ---------- 8. Somatic mutation enrichment ----------
cat("Analysing somatic mutations...\n")
# mc3 expected columns: sample, gene (or Hugo_Symbol), variant_classification
if (is.data.frame(mc3)) {
  gene_col  <- intersect(c("Hugo_Symbol","gene","Gene"), names(mc3))[1]
  samp_col  <- intersect(c("sample","Sample","Tumor_Sample_Barcode"), names(mc3))[1]
  mc3$gene  <- mc3[[gene_col]]
  mc3$samp  <- substr(mc3[[samp_col]], 1, 12)  # 12-char TCGA barcode
  # Count mutation frequency for MAA genes
  maf_maa <- mc3[mc3$gene %in% found_genes, ]
  mut_freq <- maf_maa %>%
    group_by(gene) %>%
    summarise(n_mutated_samples = n_distinct(samp), .groups = "drop")
  total_samples <- n_distinct(substr(colnames(tub_tumor), 1, 12))
  mut_freq$mut_rate_pct <- round(100 * mut_freq$n_mutated_samples / total_samples, 2)
  mut_freq <- left_join(mut_freq, gene_layer_table, by = "gene") %>%
    arrange(desc(mut_rate_pct))
} else {
  mut_freq <- data.frame(gene = found_genes, n_mutated_samples = NA,
                          mut_rate_pct = NA, layer = gene_layer_table[found_genes,"layer"])
}

## ---------- 9. Copy number amplification/deletion ----------
gistic_maa <- gistic[intersect(found_genes, rownames(gistic)), , drop = FALSE]
cn_summary <- data.frame(
  gene     = rownames(gistic_maa),
  pct_amp  = round(100 * rowMeans(gistic_maa == 2,  na.rm = TRUE), 2),
  pct_gain = round(100 * rowMeans(gistic_maa == 1,  na.rm = TRUE), 2),
  pct_del  = round(100 * rowMeans(gistic_maa == -1, na.rm = TRUE), 2),
  pct_homdel = round(100 * rowMeans(gistic_maa == -2, na.rm = TRUE), 2)
) %>%
  left_join(gene_layer_table, by = "gene") %>%
  arrange(desc(pct_amp))

## ---------- 10. CIN correlation ----------
cin_scores <- master$cin_score
names(cin_scores) <- master$barcode
cin_cor <- sapply(found_genes, function(g) {
  x <- master[[g]]; y <- master$cin_score
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 50) return(c(rho = NA, p = NA))
  ct <- cor.test(x[ok], y[ok], method = "spearman")
  c(rho = ct$estimate, p = ct$p.value)
})
cin_df <- data.frame(
  gene = found_genes,
  CIN_rho = cin_cor["rho.rho", ],
  CIN_p   = cin_cor["p", ]
)
cin_df$CIN_q <- p.adjust(cin_df$CIN_p, method = "BH")
cin_df <- left_join(cin_df, gene_layer_table, by = "gene") %>%
  arrange(desc(abs(CIN_rho)))

## ---------- 11. Save RQ2 results ----------
wb <- createWorkbook()
addWorksheet(wb, "Univariate_Cox_OS");    writeData(wb, "Univariate_Cox_OS",  surv_os)
addWorksheet(wb, "Multivariate_Cox_OS");  writeData(wb, "Multivariate_Cox_OS", multiv_os)
addWorksheet(wb, "Stage_Trends");         writeData(wb, "Stage_Trends",         stage_trend)
addWorksheet(wb, "Stage_Up_Genes");       writeData(wb, "Stage_Up_Genes",       stage_up)
addWorksheet(wb, "Stage_Down_Genes");     writeData(wb, "Stage_Down_Genes",     stage_down)
addWorksheet(wb, "Mutation_Frequency");   writeData(wb, "Mutation_Frequency",   mut_freq)
addWorksheet(wb, "Copy_Number");          writeData(wb, "Copy_Number",           cn_summary)
addWorksheet(wb, "CIN_Correlation");      writeData(wb, "CIN_Correlation",       cin_df)
saveWorkbook(wb, "results/RQ2_clinical_genomic.xlsx", overwrite = TRUE)
saveRDS(multiv_os,   "xena_data/RQ2_survival.rds")
saveRDS(stage_trend, "xena_data/RQ2_stage_trend.rds")

cat("RQ2 complete. Outputs:\n")
cat("  results/RQ2_clinical_genomic.xlsx\n")
cat("  figures/Fig3_survival_forest_OS.pdf\n")
cat("  figures/RQ2_stage_trends.pdf\n")
cat("  figures/KM_OS_*.pdf\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
