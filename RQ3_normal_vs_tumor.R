# ============================================================
# RQ0 — MAA INITIATION: NORMAL → TUMOR TRANSITION
# Question: Which MAA genes are dysregulated at cancer initiation?
# Data: TCGA matched normal samples (solid normal -11A/-11B)
# Output: Differential expression table, volcano plots, heatmap
# ============================================================
# Addresses: "How MAA evolves from normal cells to cancer cells"
# Prereq: 02_build_master.R completed

source("R startup/00_gene_panel.R")
library(ggplot2); library(ggrepel); library(pheatmap); library(dplyr)
library(openxlsx); library(limma)

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

## ---------- 1. Load data ----------
tub_tumor  <- readRDS("xena_data/tub_mat_tumor.rds")
tub_normal <- readRDS("xena_data/tub_mat_normal.rds")
master     <- readRDS("xena_data/master.rds")
pheno      <- readRDS("xena_data/tcga_sample_types.rds")

found_genes <- rownames(tub_tumor)
names(pheno) <- tolower(names(pheno))
cancer_col <- intersect(c("_primary_disease","cancer_type","_cohort"), names(pheno))[1]
pheno$cancer_type <- toupper(gsub("TCGA-|tcga-", "", as.character(pheno[[cancer_col]])))

## ---------- 2. Annotate normal samples with cancer type ----------
norm_barcodes <- colnames(tub_normal)

# Match normal samples to cancer type via phenotype table (NOT barcode position 2,
# which gives numeric TSS codes not matching master$cancer_type abbreviations)
pheno_sample_col <- names(pheno)[1]   # first column = sample barcode
norm_pheno  <- pheno[match(substr(norm_barcodes, 1, 15),
                            substr(pheno[[pheno_sample_col]], 1, 15)), ]
norm_cancer <- norm_pheno$cancer_type

# Fallback for unmatched: check master's cancer_type via patient barcode
unmatched <- is.na(norm_cancer)
if (any(unmatched)) {
  # Try matching to tumor samples from the same patient (first 12 chars)
  patient_12 <- substr(norm_barcodes[unmatched], 1, 12)
  master_12  <- substr(master$barcode, 1, 12)
  m_idx      <- match(patient_12, master_12)
  norm_cancer[unmatched] <- master$cancer_type[m_idx]
}
norm_cancer[is.na(norm_cancer)] <- "UNKNOWN"

cat("Normal sample cancer type distribution:\n")
print(sort(table(norm_cancer), decreasing = TRUE))

# Cancer types with ≥5 matched normals AND ≥5 tumor samples
norm_count     <- table(norm_cancer)
valid_ct_tumor <- unique(master$cancer_type)
ct_with_normals <- names(norm_count)[norm_count >= 5 &
                                       names(norm_count) %in% valid_ct_tumor]
cat("\nCancer types usable for DE (≥5 normal + tumor):", length(ct_with_normals), "\n")
cat(paste(ct_with_normals, collapse = ", "), "\n")

## ---------- 3. Pan-cancer differential expression (tumor vs normal) ----------
# Use limma-voom framework
de_results <- list()

for (ct in ct_with_normals) {
  t_idx <- which(master$cancer_type == ct)
  n_idx <- which(norm_cancer == ct)
  if (length(t_idx) < 5 || length(n_idx) < 5) next

  t_mat <- tub_tumor[, t_idx, drop = FALSE]
  n_mat <- tub_normal[, n_idx, drop = FALSE]
  combined <- cbind(t_mat, n_mat)
  group    <- factor(c(rep("tumor", ncol(t_mat)), rep("normal", ncol(n_mat))),
                      levels = c("normal", "tumor"))

  design <- model.matrix(~group)
  fit    <- lmFit(combined, design)
  fit    <- eBayes(fit)
  tt     <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$gene        <- rownames(tt)
  tt$cancer_type <- ct
  de_results[[ct]] <- tt
}

de_all <- do.call(rbind, de_results)

## ---------- 4. Summary: pan-cancer log2FC per gene ----------
de_summary <- de_all %>%
  group_by(gene) %>%
  summarise(
    mean_logFC     = mean(logFC,       na.rm = TRUE),
    median_logFC   = median(logFC,     na.rm = TRUE),
    n_up           = sum(logFC > 1  & adj.P.Val < 0.05, na.rm = TRUE),
    n_down         = sum(logFC < -1 & adj.P.Val < 0.05, na.rm = TRUE),
    n_sig          = sum(adj.P.Val < 0.05, na.rm = TRUE),
    n_cancers_tested = n(),
    .groups = "drop"
  ) %>%
  mutate(pct_sig = round(100 * n_sig / n_cancers_tested, 1)) %>%
  left_join(gene_layer_table, by = "gene") %>%
  arrange(desc(abs(mean_logFC)))

## ---------- 5. Figure: volcano plot (pan-cancer mean logFC) ----------
de_summary$label <- ifelse(abs(de_summary$mean_logFC) > 1 | de_summary$n_sig >= 5,
                            de_summary$gene, "")

p_volcano <- ggplot(de_summary, aes(x = mean_logFC, y = pct_sig,
                                    colour = layer, label = label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = layer_colors, na.value = "grey60") +
  labs(title = "RQ0: MAA Normal vs Tumor — Pan-cancer Fold Change",
       subtitle = paste0("n = ", length(ct_with_normals), " cancer types with matched normals"),
       x = "Mean log2 Fold Change (Tumor / Normal)",
       y = "% Cancer Types Significantly Dysregulated") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")
ggsave("figures/RQ0_volcano_tumor_vs_normal.pdf", p_volcano, width = 10, height = 7)

## ---------- 6. Figure: heatmap of log2FC across cancer types ----------
# Pivot: genes × cancer types, fill = logFC
fc_mat <- de_all %>%
  select(gene, cancer_type, logFC) %>%
  tidyr::pivot_wider(names_from = cancer_type, values_from = logFC) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Cap values for display
fc_mat_capped <- pmin(pmax(fc_mat, -4), 4)

ann_row <- gene_layer_table[rownames(fc_mat_capped), "layer", drop = FALSE]
colnames(ann_row) <- "Layer"
ann_colors <- list(Layer = layer_colors)

pheatmap::pheatmap(
  fc_mat_capped,
  annotation_row     = ann_row,
  annotation_colors  = ann_colors,
  clustering_method  = "ward.D2",
  color              = colorRampPalette(c("#053061","#2166AC","#F7F7F7","#D6604D","#67001F"))(100),
  breaks             = seq(-4, 4, length.out = 101),
  show_colnames      = TRUE,
  show_rownames      = TRUE,
  fontsize_row       = 7,
  fontsize_col       = 8,
  main               = "RQ0: MAA log2FC (Tumor vs Normal) by Cancer Type",
  filename           = "figures/RQ0_heatmap_tumor_vs_normal.pdf",
  width = 14, height = 16
)

## ---------- 7. Top consistently upregulated MAA genes (initiation signature) ----------
initiation_sig <- de_summary %>%
  filter(mean_logFC > 0.5, pct_sig >= 50) %>%
  arrange(desc(mean_logFC)) %>%
  select(gene, layer, mean_logFC, median_logFC, n_up, n_sig,
         n_cancers_tested, pct_sig)

cat("\nTop MAA initiation genes (upregulated in ≥50% of cancer types):\n")
print(head(initiation_sig, 20))

## ---------- 8. Save results ----------
wb <- createWorkbook()
addWorksheet(wb, "DE_Summary")
writeData(wb, "DE_Summary", de_summary)
addWorksheet(wb, "DE_All_CancerTypes")
writeData(wb, "DE_All_CancerTypes", de_all)
addWorksheet(wb, "Initiation_Signature")
writeData(wb, "Initiation_Signature", initiation_sig)
saveWorkbook(wb, "results/RQ0_normal_vs_tumor.xlsx", overwrite = TRUE)

saveRDS(de_summary,    "xena_data/RQ0_de_summary.rds")
saveRDS(initiation_sig,"xena_data/RQ0_initiation_signature.rds")

cat("\nRQ0 complete. Outputs:\n")
cat("  results/RQ0_normal_vs_tumor.xlsx\n")
cat("  figures/RQ0_volcano_tumor_vs_normal.pdf\n")
cat("  figures/RQ0_heatmap_tumor_vs_normal.pdf\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
