# ============================================================
# MAA PAN-CANCER ANALYSIS — MASTER RUNNER v2.0
# Run this file from the R project root to execute the full pipeline.
# Each section is self-contained and can be run independently.
# ============================================================
# Estimated runtimes (8-core workstation, 32 GB RAM):
#   Setup     : <1 min
#   Download  : 20–40 min (1.5 GB download + matrix parse)
#   Master    : 5–10 min
#   RQ0       : 10–15 min
#   RQ1       : 5 min
#   RQ2       : 15–20 min
#   RQ3       : 10–15 min
#   RQ4       : 10 min
#   RQ5       : 10 min
#   Total     : ~1.5–2 h (excl. download)

## ---- STEP 0: Install required packages (run once) ----
required_packages <- c(
  # Core data handling
  "data.table","readr","openxlsx","R.utils","tidyr","dplyr","tibble",
  # Matrix math
  "matrixStats",
  # Visualisation
  "ggplot2","ggrepel","pheatmap","RColorBrewer","ggraph","corrplot",
  "survminer","patchwork",
  # Statistics / survival
  "survival","limma",
  # Network
  "igraph",
  # Enrichment
  "clusterProfiler","org.Hs.eg.db","ReactomePA","msigdbr","enrichplot"
)
new_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly=TRUE)]
if (length(new_pkgs) > 0) {
  message("Installing: ", paste(new_pkgs, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  bioc_pkgs <- c("limma","clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot","msigdbr")
  cran_new  <- setdiff(new_pkgs, bioc_pkgs)
  bioc_new  <- intersect(new_pkgs, bioc_pkgs)
  if (length(cran_new)  > 0) install.packages(cran_new)
  if (length(bioc_new)  > 0) BiocManager::install(bioc_new, update = FALSE)
}
message("All packages available.")

## ---- STEP 1: Gene panel definition ----
source("R startup/00_gene_panel.R")

## ---- STEP 2: Download data (ONCE only — skip if already done) ----
if (!file.exists("xena_data/rna_tpm_pancan.rds")) {
  source("R startup/01_download_data.R")
} else {
  message("TPM data already downloaded: xena_data/rna_tpm_pancan.rds")
}

## ---- STEP 3: Build master table ----
if (!file.exists("xena_data/master.rds")) {
  source("R startup/02_build_master.R")
} else {
  message("Master table already exists. Loading...")
  source("R startup/00_gene_panel.R")
  master     <- readRDS("xena_data/master.rds")
  tub_tumor  <- readRDS("xena_data/tub_mat_tumor.rds")
  tub_normal <- readRDS("xena_data/tub_mat_normal.rds")
  gene_layer_table <- readRDS("xena_data/gene_layer_table.rds")
  found_genes <- rownames(tub_tumor)
}

## ---- STEP 4: RQ0 — Normal vs Tumor (MAA Initiation) ----
message("\n=== Running RQ0: Normal vs Tumor ===")
source("RQ_analysis/RQ0_normal_vs_tumor.R")

## ---- STEP 5: RQ1 — Pan-cancer Expression Atlas ----
message("\n=== Running RQ1: Expression Atlas ===")
source("RQ_analysis/RQ1_expression_atlas.R")

## ---- STEP 6: RQ2 — Clinical, Genomic & Stage Analysis ----
message("\n=== Running RQ2: Clinical & Genomic ===")
source("RQ_analysis/RQ2_clinical_genomic.R")

## ---- STEP 7: RQ3 — Fitness Traits + Drug Resistance ----
message("\n=== Running RQ3: Fitness Traits & Drug Resistance ===")
source("RQ_analysis/RQ3_fitness_drug_resistance.R")

## ---- STEP 8: RQ4 — Network & Pathway Enrichment ----
message("\n=== Running RQ4: Network & Pathways ===")
source("RQ_analysis/RQ4_network_pathways.R")

## ---- STEP 9: RQ5 — Adaptive States & Evolution Trajectory ----
message("\n=== Running RQ5: Adaptive States & Evolution ===")
source("RQ_analysis/RQ5_adaptive_states_evolution.R")

## ---- STEP 10: Final summary ----
priority <- readRDS("xena_data/RQ5_priority_candidates.rds")
message("\n============================================================")
message("ANALYSIS COMPLETE — MAA Pan-Cancer Pipeline v2.0")
message("============================================================")
message("Top 10 MAA Priority Candidates (MEI score):")
print(head(priority[, c("gene","layer","MEI_score","C6_drug_resistance",
                          "n_resist_traits")], 10))
message("\nAll results saved in: results/")
message("All figures saved in: figures/")
message("R objects saved in:   xena_data/*.rds")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
