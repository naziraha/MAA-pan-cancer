# ═══════════════════════════════════════════════════════════════════════════
# MAA Pan-Cancer Study — GEO External Cohort Validation
# 12-Gene Validation Set
# Platforms: GSE14520 (LIHC), GSE96058 (BRCA), GSE72094 (LUAD),
#            GSE62254 (STAD), GSE39582 (COAD)
# ═══════════════════════════════════════════════════════════════════════════

# ── Step 1: Install packages (run once) ────────────────────────────────────
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GEOquery","limma","Biobase"))
# install.packages(c("survival","survminer","dplyr","ggplot2","writexl"))

# ── Step 2: Load libraries ─────────────────────────────────────────────────
library(GEOquery)
library(survival)
library(survminer)
library(limma)
library(dplyr)
library(ggplot2)
library(writexl)

# ── Step 3: Define 12 validation genes ────────────────────────────────────
val_genes <- c(
  "KIF20A", "BIRC5",  "KIF14",  "KIF23",  "KIF2C",  "CDK1",
  "CENPE",  "KIF18A", "TUBA1B", "FSCN1",  "TUBB3",  "TTL"
)

# ── Step 4: Define GEO datasets ───────────────────────────────────────────
gse_info <- data.frame(
  GSE       = c("GSE14520","GSE96058","GSE72094","GSE62254","GSE39582"),
  Cancer    = c("LIHC",    "BRCA",    "LUAD",    "STAD",    "COAD"),
  Platform  = c("Microarray","RNA-seq","Microarray","Microarray","Microarray"),
  n_samples = c(247,       3273,      442,       300,       585),
  stringsAsFactors = FALSE
)
print(gse_info)

# ── Step 5: Core validation function ──────────────────────────────────────
validate_gene_gse <- function(gse_id, gene_sym, verbose=TRUE) {

  if (verbose) message("  --> Processing: ", gse_id, " | Gene: ", gene_sym)

  # Download GSE
  gse_raw <- tryCatch(
    getGEO(gse_id, GSEMatrix=TRUE, AnnotGPL=TRUE),
    error = function(e) { message("Download failed: ", e$message); return(NULL) }
  )
  if (is.null(gse_raw)) return(NULL)
  gse <- gse_raw[[1]]

  em <- exprs(gse)          # expression matrix
  ph <- pData(gse)          # phenotype data
  fd <- fData(gse)          # feature/probe data

  # ── Find probes for gene ──
  gene_col <- grep("gene.symbol|gene_symbol|symbol|gene.name",
                   colnames(fd), ignore.case=TRUE, value=TRUE)[1]
  if (is.na(gene_col)) {
    message("    Cannot find gene symbol column in fData")
    return(data.frame(GSE=gse_id, Cancer=NA, Gene=gene_sym,
                      n=NA, p=NA, Significant=NA, Direction=NA, Note="No gene col"))
  }

  idx <- which(toupper(fd[[gene_col]]) == toupper(gene_sym))
  if (length(idx) == 0) {
    if (verbose) message("    Gene not found: ", gene_sym)
    return(data.frame(GSE=gse_id, Cancer=NA, Gene=gene_sym,
                      n=NA, p=NA, Significant=NA, Direction=NA, Note="Gene not found"))
  }

  # Average multiple probes
  gene_exp <- if (length(idx) > 1) colMeans(em[idx,,drop=FALSE], na.rm=TRUE) else em[idx,]

  # ── Assign high/low groups by median ──
  ph$gene_grp <- ifelse(gene_exp > median(gene_exp, na.rm=TRUE), "High", "Low")

  # ── Auto-detect survival columns ──
  os_time_col <- grep(
    "os[._]time|overall.*surv.*time|surv.*month|days.*surv|follow.*up.*day|months.*follow",
    colnames(ph), ignore.case=TRUE, value=TRUE)[1]

  os_event_col <- grep(
    "os[._]status|os[._]event|vital.*status|death|event.*os|censored",
    colnames(ph), ignore.case=TRUE, value=TRUE)[1]

  # Fallback: check characteristics columns
  if (is.na(os_time_col)) {
    char_cols <- grep("characteristics", colnames(ph), ignore.case=TRUE, value=TRUE)
    os_time_col  <- char_cols[grep("time|month|day", char_cols, ignore.case=TRUE)][1]
    os_event_col <- char_cols[grep("status|event|death|vital", char_cols, ignore.case=TRUE)][1]
  }

  if (is.na(os_time_col) || is.na(os_event_col)) {
    if (verbose) message("    No survival columns found. Available: ", paste(colnames(ph)[1:10], collapse=", "))
    return(data.frame(GSE=gse_id, Cancer=NA, Gene=gene_sym,
                      n=NA, p=NA, Significant=NA, Direction=NA, Note="No survival data"))
  }

  # Convert to numeric
  ph$os_time  <- suppressWarnings(as.numeric(as.character(ph[[os_time_col]])))
  ph$os_event <- suppressWarnings(as.numeric(as.character(ph[[os_event_col]])))

  # Remove NAs and keep only tumour samples (if os_time > 0)
  ph_clean <- ph[!is.na(ph$os_time) & !is.na(ph$os_event) & ph$os_time > 0, ]

  if (nrow(ph_clean) < 20) {
    return(data.frame(GSE=gse_id, Cancer=NA, Gene=gene_sym,
                      n=nrow(ph_clean), p=NA, Significant=NA, Direction=NA,
                      Note="Too few samples"))
  }

  # ── Kaplan-Meier + log-rank test ──
  fit   <- survfit(Surv(os_time, os_event) ~ gene_grp, data=ph_clean)
  pval  <- surv_pvalue(fit, ph_clean)$pval

  # Direction: does High expression = shorter survival?
  med_high <- median(ph_clean$os_time[ph_clean$gene_grp == "High"], na.rm=TRUE)
  med_low  <- median(ph_clean$os_time[ph_clean$gene_grp == "Low"],  na.rm=TRUE)
  direction <- ifelse(med_high < med_low, "Unfavorable (High=worse)", "Favorable (High=better)")

  # Get cancer type from gse_info
  cancer <- gse_info$Cancer[gse_info$GSE == gse_id]
  if (length(cancer)==0) cancer <- NA

  return(data.frame(
    GSE         = gse_id,
    Cancer      = cancer,
    Gene        = gene_sym,
    n           = nrow(ph_clean),
    p           = round(pval, 4),
    Significant = ifelse(!is.na(pval) & pval < 0.05, "YES", "NO"),
    Direction   = direction,
    Note        = "OK",
    stringsAsFactors = FALSE
  ))
}

# ── Step 6: Run all GSE × gene combinations ───────────────────────────────
message("\n========================================")
message("Starting GEO validation: ",
        length(gse_info$GSE), " datasets x ", length(val_genes), " genes")
message("========================================\n")

all_results <- list()

for (gse_id in gse_info$GSE) {
  message("\n--- Processing dataset: ", gse_id, " ---")
  for (gene in val_genes) {
    res <- tryCatch(
      validate_gene_gse(gse_id, gene),
      error = function(e) {
        data.frame(GSE=gse_id, Cancer=NA, Gene=gene,
                   n=NA, p=NA, Significant="ERROR", Direction=NA,
                   Note=as.character(e$message), stringsAsFactors=FALSE)
      }
    )
    if (!is.null(res)) all_results[[paste(gse_id, gene, sep="_")]] <- res
  }
}

# ── Step 7: Compile results ────────────────────────────────────────────────
results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

message("\n========================================")
message("RESULTS SUMMARY")
message("========================================")
print(results_df)

# ── Step 8: Summary pivot table ───────────────────────────────────────────
summary_wide <- results_df %>%
  select(GSE, Cancer, Gene, Significant, p) %>%
  tidyr::pivot_wider(names_from=Gene, values_from=Significant)

message("\nSignificance Summary (YES = p<0.05):")
print(summary_wide)

# Count validations per gene
gene_summary <- results_df %>%
  filter(!is.na(p)) %>%
  group_by(Gene) %>%
  summarise(
    n_cohorts_tested = n(),
    n_significant    = sum(Significant=="YES", na.rm=TRUE),
    pct_validated    = round(100*n_significant/n_cohorts_tested, 1),
    direction_unfav  = sum(grepl("Unfavorable", Direction), na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(n_significant))

message("\nGene-level validation summary:")
print(gene_summary)

# ── Step 9: Save outputs ───────────────────────────────────────────────────
output_list <- list(
  "All_Results"    = results_df,
  "Gene_Summary"   = gene_summary,
  "Wide_Summary"   = as.data.frame(summary_wide)
)

write_xlsx(output_list, "GEO_validation_results.xlsx")
write.csv(results_df, "GEO_validation_results.csv", row.names=FALSE)
message("\nResults saved: GEO_validation_results.xlsx")

# ── Step 10: Generate KM plots for significant results ────────────────────
message("\nGenerating KM plots for significant results...")

pdf("GEO_validation_KM_plots.pdf", width=8, height=6)

for (gse_id in gse_info$GSE) {
  gse_raw <- tryCatch(getGEO(gse_id, GSEMatrix=TRUE), error=function(e) NULL)
  if (is.null(gse_raw)) next
  gse <- gse_raw[[1]]
  em <- exprs(gse); ph <- pData(gse); fd <- fData(gse)
  cancer <- gse_info$Cancer[gse_info$GSE == gse_id]

  gene_col <- grep("gene.symbol|gene_symbol|symbol",
                   colnames(fd), ignore.case=TRUE, value=TRUE)[1]
  os_time_col  <- grep("os[._]time|overall.*surv|surv.*month|days.*surv",
                        colnames(ph), ignore.case=TRUE, value=TRUE)[1]
  os_event_col <- grep("os[._]status|vital.*status|death|event",
                        colnames(ph), ignore.case=TRUE, value=TRUE)[1]

  if (is.na(gene_col) || is.na(os_time_col) || is.na(os_event_col)) next

  ph$os_time  <- suppressWarnings(as.numeric(ph[[os_time_col]]))
  ph$os_event <- suppressWarnings(as.numeric(ph[[os_event_col]]))
  ph <- ph[!is.na(ph$os_time) & !is.na(ph$os_event) & ph$os_time > 0, ]

  for (gene in val_genes) {
    # Check if significant
    sig_row <- results_df[results_df$GSE==gse_id & results_df$Gene==gene,]
    if (nrow(sig_row)==0 || is.na(sig_row$p) || sig_row$p >= 0.05) next

    idx <- which(toupper(fd[[gene_col]]) == toupper(gene))
    if (length(idx)==0) next
    gene_exp <- if (length(idx)>1) colMeans(em[idx,,drop=FALSE]) else em[idx,]
    gene_exp_subset <- gene_exp[rownames(ph)]
    ph_plot <- ph
    ph_plot$gene_grp <- ifelse(gene_exp_subset > median(gene_exp_subset,na.rm=TRUE),
                               "High","Low")

    fit <- survfit(Surv(os_time, os_event) ~ gene_grp, data=ph_plot)
    p <- ggsurvplot(fit, data=ph_plot,
                    pval=TRUE, risk.table=TRUE,
                    palette=c("#C00000","#2E75B6"),
                    title=paste0(gene, " — ", cancer, " (", gse_id, ")"),
                    xlab="Time", ylab="Overall Survival",
                    legend.labs=c("High","Low"),
                    ggtheme=theme_minimal())
    print(p)
  }
}
dev.off()
message("KM plots saved: GEO_validation_KM_plots.pdf")

message("\n========================================")
message("GEO VALIDATION COMPLETE")
message("========================================")
