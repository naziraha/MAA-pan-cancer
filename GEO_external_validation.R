# ============================================================
# MAA PAN-CANCER — TIER 3 EXTERNAL VALIDATION
# GEO Independent Cohort Validation
# Cohorts:
#   GSE68465 — LUAD (n=443, OS)
#   GSE14520 — HCC/LIHC (n=225, OS)
#   GSE2034  — BRCA (n=286, RFS)
# Purpose: Validate top MAA priority candidates (MEI top 10)
#          in datasets independent of TCGA
# ============================================================

source("R startup/00_gene_panel.R")

## ---------- 0. Install/load packages ----------
pkgs <- c("GEOquery","survival","survminer","ggplot2","dplyr",
          "openxlsx","ggrepel","tidyr","Biobase")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("GEOquery","Biobase")) {
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else install.packages(p)
  }
}))
library(GEOquery); library(survival); library(survminer)
library(ggplot2);  library(dplyr);    library(openxlsx)
library(ggrepel);  library(tidyr);    library(Biobase)

if (!dir.exists("Validation")) dir.create("Validation")
if (!dir.exists("figures"))    dir.create("figures")

## ---------- 1. Top MAA priority genes to validate ----------
# Load MEI rankings from RQ5 output
if (file.exists("xena_data/RQ5_priority_candidates.rds")) {
  priority <- readRDS("xena_data/RQ5_priority_candidates.rds")
  top_genes <- head(priority$gene, 15)
} else {
  # Fallback: use known top candidates
  top_genes <- c("GTSE1","KIFC1","KIF18B","CDK1","BIRC5",
                 "KIF23","KIF20A","KIF2C","KIF14","CENPE",
                 "SPAG5","KIF11","KIF15","TUBA1B","TUBA1C")
}
cat("Validating genes:", paste(top_genes, collapse=", "), "\n")

## ================================================================
## HELPER FUNCTIONS
## ================================================================

## Extract expression matrix and clinical data from a GEO dataset
load_geo_dataset <- function(gse_id, cache_dir = "Validation") {
  rds_file <- file.path(cache_dir, paste0(gse_id, ".rds"))
  if (file.exists(rds_file)) {
    cat("Loading cached:", gse_id, "\n")
    return(readRDS(rds_file))
  }
  cat("Downloading", gse_id, "...\n")
  options(timeout = 600)
  gse <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE,
                destdir = cache_dir)
  saveRDS(gse, rds_file)
  gse
}

## Map probe IDs to gene symbols and collapse to gene-level (mean)
probes_to_genes <- function(expr_mat, feature_df, symbol_col = "Gene.symbol") {
  if (!symbol_col %in% colnames(feature_df)) {
    alt <- grep("symbol|gene.?name|gene_symbol", colnames(feature_df),
                ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(alt)) symbol_col <- alt else stop("No gene symbol column found")
  }
  feature_df$symbol <- trimws(as.character(feature_df[[symbol_col]]))
  # Remove probes with empty or multi-mapping symbols
  feature_df$symbol[grepl("///|^$", feature_df$symbol)] <- NA
  keep <- !is.na(feature_df$symbol)
  expr_filt <- expr_mat[keep, , drop = FALSE]
  syms      <- feature_df$symbol[keep]
  # Collapse: keep probe with highest mean expression per gene
  gene_means <- rowMeans(expr_filt, na.rm = TRUE)
  idx        <- tapply(seq_along(syms), syms,
                       function(i) i[which.max(gene_means[i])])
  expr_gene  <- expr_filt[unlist(idx), , drop = FALSE]
  rownames(expr_gene) <- names(idx)
  expr_gene
}

## Run KM + univariate Cox for a single gene in a cohort
run_survival_gene <- function(gene, expr_mat, surv_time, surv_event,
                               cohort_name, time_unit = "days") {
  if (!gene %in% rownames(expr_mat)) return(NULL)
  gexpr <- as.numeric(expr_mat[gene, ])
  ok    <- !is.na(gexpr) & !is.na(surv_time) & !is.na(surv_event) & surv_time > 0
  if (sum(ok) < 30) return(NULL)
  df <- data.frame(expr  = gexpr[ok],
                   time  = as.numeric(surv_time[ok]),
                   event = as.numeric(surv_event[ok]))
  # Convert to years
  if (time_unit == "days")  df$time_yr <- df$time / 365.25
  if (time_unit == "months") df$time_yr <- df$time / 12
  if (time_unit == "years")  df$time_yr <- df$time
  # Cox regression
  fit <- tryCatch(
    coxph(Surv(time_yr, event) ~ expr, data = df),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coef
  # KM (median split)
  df$group <- ifelse(df$expr >= median(df$expr), "High", "Low")
  km_fit   <- survfit(Surv(time_yr, event) ~ group, data = df)
  lr       <- survdiff(Surv(time_yr, event) ~ group, data = df)
  lr_p     <- 1 - pchisq(lr$chisq, df = 1)
  data.frame(
    gene        = gene,
    cohort      = cohort_name,
    n           = nrow(df),
    HR          = exp(s[1,1]),
    HR_lo       = exp(s[1,1] - 1.96*s[1,3]),
    HR_hi       = exp(s[1,1] + 1.96*s[1,3]),
    cox_p       = s[1,5],
    km_logrank_p = lr_p,
    stringsAsFactors = FALSE
  )
}

## ================================================================
## COHORT 1: GSE68465 — LUAD (n=443, OS)
## ================================================================
cat("\n--- Processing GSE68465 (LUAD) ---\n")
tryCatch({
  gse68465 <- load_geo_dataset("GSE68465")
  eset68465 <- gse68465[[1]]

  expr68465 <- exprs(eset68465)
  feat68465 <- fData(eset68465)
  phno68465 <- pData(eset68465)

  # Gene-level expression
  gene_expr_68465 <- probes_to_genes(expr68465, feat68465)
  cat("GSE68465: ", nrow(gene_expr_68465), "genes x", ncol(gene_expr_68465), "samples\n")

  # Extract survival — column names vary; try common patterns
  time_col  <- grep("os.*time|survival.*time|time.*overall|days.*survival|months",
                     colnames(phno68465), ignore.case = TRUE, value = TRUE)[1]
  event_col <- grep("os.*event|overall.*surv|vital|status|censor|event",
                     colnames(phno68465), ignore.case = TRUE, value = TRUE)[1]
  cat("GSE68465 survival cols:", time_col, "|", event_col, "\n")

  surv_time_68465  <- as.numeric(as.character(phno68465[[time_col]]))
  surv_event_68465 <- as.numeric(as.character(phno68465[[event_col]]))
  # Recode event if needed (1=dead/event, 0=alive/censored)
  if (max(surv_event_68465, na.rm = TRUE) > 1)
    surv_event_68465 <- ifelse(surv_event_68465 == 2, 1, 0)

  res68465 <- do.call(rbind, lapply(top_genes, run_survival_gene,
                                     expr_mat = gene_expr_68465,
                                     surv_time = surv_time_68465,
                                     surv_event = surv_event_68465,
                                     cohort_name = "GSE68465_LUAD",
                                     time_unit = "days"))
  cat("GSE68465 results:", nrow(res68465), "genes\n")
}, error = function(e) {
  cat("GSE68465 error:", conditionMessage(e), "\n")
  res68465 <<- NULL
})

## ================================================================
## COHORT 2: GSE14520 — HCC/LIHC (n=225, OS)
## ================================================================
cat("\n--- Processing GSE14520 (HCC) ---\n")
tryCatch({
  gse14520 <- load_geo_dataset("GSE14520")
  eset14520 <- gse14520[[1]]

  expr14520 <- exprs(eset14520)
  feat14520 <- fData(eset14520)
  phno14520 <- pData(eset14520)

  gene_expr_14520 <- probes_to_genes(expr14520, feat14520)
  cat("GSE14520: ", nrow(gene_expr_14520), "genes x", ncol(gene_expr_14520), "samples\n")

  time_col  <- grep("os.*time|survival.*time|time.*overall|overall.*surv.*time",
                     colnames(phno14520), ignore.case = TRUE, value = TRUE)[1]
  event_col <- grep("os.*event|overall.*surv.*event|vital|status|censor|event",
                     colnames(phno14520), ignore.case = TRUE, value = TRUE)[1]
  cat("GSE14520 survival cols:", time_col, "|", event_col, "\n")

  surv_time_14520  <- as.numeric(as.character(phno14520[[time_col]]))
  surv_event_14520 <- as.numeric(as.character(phno14520[[event_col]]))
  if (!is.na(max(surv_event_14520, na.rm = TRUE)) &&
      max(surv_event_14520, na.rm = TRUE) > 1)
    surv_event_14520 <- ifelse(surv_event_14520 == 2, 1, 0)

  res14520 <- do.call(rbind, lapply(top_genes, run_survival_gene,
                                     expr_mat = gene_expr_14520,
                                     surv_time = surv_time_14520,
                                     surv_event = surv_event_14520,
                                     cohort_name = "GSE14520_HCC",
                                     time_unit = "months"))
  cat("GSE14520 results:", nrow(res14520), "genes\n")
}, error = function(e) {
  cat("GSE14520 error:", conditionMessage(e), "\n")
  res14520 <<- NULL
})

## ================================================================
## COHORT 3: GSE2034 — BRCA (n=286, RFS)
## ================================================================
cat("\n--- Processing GSE2034 (BRCA) ---\n")
tryCatch({
  gse2034 <- load_geo_dataset("GSE2034")
  eset2034 <- gse2034[[1]]

  expr2034 <- exprs(eset2034)
  feat2034 <- fData(eset2034)
  phno2034 <- pData(eset2034)

  gene_expr_2034 <- probes_to_genes(expr2034, feat2034)
  cat("GSE2034: ", nrow(gene_expr_2034), "genes x", ncol(gene_expr_2034), "samples\n")

  time_col  <- grep("rfs.*time|relapse.*time|time.*relapse|dmfs.*time",
                     colnames(phno2034), ignore.case = TRUE, value = TRUE)[1]
  event_col <- grep("rfs.*event|relapse.*event|event.*relapse|dmfs.*event",
                     colnames(phno2034), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(time_col))  time_col  <- colnames(phno2034)[grep("time",   colnames(phno2034), ignore.case=TRUE)[1]]
  if (is.na(event_col)) event_col <- colnames(phno2034)[grep("event|status", colnames(phno2034), ignore.case=TRUE)[1]]
  cat("GSE2034 survival cols:", time_col, "|", event_col, "\n")

  surv_time_2034  <- as.numeric(as.character(phno2034[[time_col]]))
  surv_event_2034 <- as.numeric(as.character(phno2034[[event_col]]))
  if (!is.na(max(surv_event_2034, na.rm = TRUE)) &&
      max(surv_event_2034, na.rm = TRUE) > 1)
    surv_event_2034 <- ifelse(surv_event_2034 == 2, 1, 0)

  res2034 <- do.call(rbind, lapply(top_genes, run_survival_gene,
                                    expr_mat = gene_expr_2034,
                                    surv_time = surv_time_2034,
                                    surv_event = surv_event_2034,
                                    cohort_name = "GSE2034_BRCA",
                                    time_unit = "months"))
  cat("GSE2034 results:", nrow(res2034), "genes\n")
}, error = function(e) {
  cat("GSE2034 error:", conditionMessage(e), "\n")
  res2034 <<- NULL
})

## ================================================================
## COMBINE & SUMMARISE VALIDATION RESULTS
## ================================================================
cat("\n--- Combining results ---\n")
all_results <- do.call(rbind, Filter(Negate(is.null),
                                      list(res68465, res14520, res2034)))

if (!is.null(all_results) && nrow(all_results) > 0) {
  all_results$q        <- p.adjust(all_results$cox_p, method = "BH")
  all_results$sig      <- all_results$q < 0.05
  all_results$direction <- ifelse(all_results$HR > 1, "Risk", "Protective")
  all_results <- left_join(all_results, gene_layer_table, by = "gene")

  ## ---------- Validation concordance with TCGA ----------
  tcga_surv <- readRDS("xena_data/RQ2_survival.rds")
  all_results$TCGA_HR  <- tcga_surv$HR[match(all_results$gene, tcga_surv$gene)]
  all_results$TCGA_q   <- tcga_surv$q[match(all_results$gene, tcga_surv$gene)]
  all_results$concordant <- with(all_results,
    !is.na(TCGA_HR) & !is.na(HR) &
    ((HR > 1 & TCGA_HR > 1) | (HR < 1 & TCGA_HR < 1)))

  ## Summary table: per gene, how many cohorts validate?
  validation_summary <- all_results %>%
    group_by(gene, layer) %>%
    summarise(
      n_cohorts_tested    = n(),
      n_cohorts_sig       = sum(sig, na.rm = TRUE),
      n_cohorts_concordant = sum(concordant, na.rm = TRUE),
      mean_HR_GEO         = round(mean(HR, na.rm = TRUE), 3),
      TCGA_HR             = round(mean(TCGA_HR, na.rm = TRUE), 3),
      TCGA_q              = round(mean(TCGA_q, na.rm = TRUE), 4),
      validation_rate_pct = round(100 * n_cohorts_sig / n_cohorts_tested, 1),
      .groups = "drop"
    ) %>%
    mutate(validated = n_cohorts_concordant >= 1 & n_cohorts_sig >= 1) %>%
    arrange(desc(n_cohorts_sig), desc(n_cohorts_concordant))

  cat("\n=== VALIDATION SUMMARY ===\n")
  print(validation_summary)

  n_validated <- sum(validation_summary$validated, na.rm = TRUE)
  cat("\nGenes validated in ≥1 external GEO cohort:", n_validated, "/",
      nrow(validation_summary), "\n")

  ## ---------- Forest plot: GEO validation ----------
  plot_df <- all_results %>%
    filter(!is.na(HR), !is.na(HR_lo), !is.na(HR_hi)) %>%
    mutate(gene_cohort = paste0(gene, "\n(", cohort, ")"))

  p_forest_val <- ggplot(plot_df,
    aes(x = HR, xmin = HR_lo, xmax = HR_hi,
        y = reorder(gene_cohort, HR), colour = cohort)) +
    geom_pointrange(size = 0.6) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_x_log10() +
    scale_colour_brewer(palette = "Set1") +
    labs(title = "External Validation: MAA Priority Genes — GEO Cohorts",
         subtitle = "GSE68465 (LUAD) | GSE14520 (HCC) | GSE2034 (BRCA)",
         x = "Hazard Ratio (95% CI)", y = NULL, colour = "GEO Cohort") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  ggsave("figures/Validation_GEO_forest.pdf", p_forest_val,
         width = 10, height = 12)

  ## ---------- Concordance heatmap (TCGA vs GEO) ----------
  wide_hr <- all_results %>%
    select(gene, cohort, HR) %>%
    pivot_wider(names_from = cohort, values_from = HR) %>%
    tibble::column_to_rownames("gene")

  if (ncol(wide_hr) > 0) {
    # Log2HR for display
    log2hr_mat <- log2(as.matrix(wide_hr))
    log2hr_mat[!is.finite(log2hr_mat)] <- 0
    # Add TCGA column
    tcga_col <- log2(tcga_surv$HR[match(rownames(log2hr_mat), tcga_surv$gene)])
    tcga_col[is.na(tcga_col) | !is.finite(tcga_col)] <- 0
    log2hr_mat <- cbind(log2hr_mat, TCGA = tcga_col)

    pheatmap::pheatmap(
      pmin(pmax(log2hr_mat, -2), 2),
      color = colorRampPalette(c("#053061","#FFFFFF","#67001F"))(100),
      breaks = seq(-2, 2, length.out = 101),
      cluster_rows = TRUE, cluster_cols = FALSE,
      main = "External Validation: log2(HR) — TCGA vs GEO Cohorts",
      filename = "figures/Validation_HR_concordance_heatmap.pdf",
      width = 8, height = 10
    )
  }

  ## ---------- Save results ----------
  wb <- createWorkbook()
  addWorksheet(wb, "All_GEO_Results");    writeData(wb,"All_GEO_Results",   all_results)
  addWorksheet(wb, "Validation_Summary"); writeData(wb,"Validation_Summary",validation_summary)
  saveWorkbook(wb, "results/Validation_GEO_external.xlsx", overwrite = TRUE)
  saveRDS(all_results,       "xena_data/Validation_GEO_results.rds")
  saveRDS(validation_summary,"xena_data/Validation_GEO_summary.rds")

  cat("\nValidation complete. Outputs:\n")
  cat("  results/Validation_GEO_external.xlsx\n")
  cat("  figures/Validation_GEO_forest.pdf\n")
  cat("  figures/Validation_HR_concordance_heatmap.pdf\n")

} else {
  cat("No validation results generated. Check GEO download and survival column names.\n")
}


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
