# ============================================================
# RQ1 — PAN-CANCER MAA EXPRESSION ATLAS
# Question: Do MAA genes show coordinated dysregulation across
#           33 cancer types, and what modules co-regulate?
# Outputs: Heatmap (genes × cancer types), co-expression matrix,
#          module assignments, expression summary table
# ============================================================

source("R startup/00_gene_panel.R")
library(ggplot2); library(pheatmap); library(dplyr); library(openxlsx)
library(corrplot); library(RColorBrewer)

## ---------- 1. Load data ----------
master     <- readRDS("xena_data/master.rds")
tub_tumor  <- readRDS("xena_data/tub_mat_tumor.rds")
found_genes <- rownames(tub_tumor)

## ---------- 2. Per-cancer median expression ----------
cancer_types <- sort(unique(master$cancer_type))
cancer_types <- cancer_types[cancer_types != "UNKNOWN"]

med_expr <- sapply(cancer_types, function(ct) {
  idx <- which(master$cancer_type == ct)
  if (length(idx) < 10) return(rep(NA, nrow(tub_tumor)))
  rowMedians(tub_tumor[, idx, drop = FALSE])
})
rownames(med_expr) <- found_genes

# Remove cancer types with too few samples
valid_ct <- colSums(!is.na(med_expr)) == nrow(med_expr)
med_expr <- med_expr[, valid_ct, drop = FALSE]

# Remove genes where SD = 0 across cancer types (z-score would produce NaN)
gene_sd  <- apply(med_expr, 1, sd, na.rm = TRUE)
med_expr <- med_expr[gene_sd > 0 & !is.na(gene_sd), , drop = FALSE]

# Z-score each gene across cancer types (per-row)
z_expr <- t(scale(t(med_expr)))
# Replace any residual NaN/Inf with 0 (shouldn't occur after SD filter, but safety net)
z_expr[!is.finite(z_expr)] <- 0

## ---------- 3. Pan-cancer heatmap (Fig1) ----------
ann_row <- gene_layer_table[rownames(z_expr), "layer", drop = FALSE]
colnames(ann_row) <- "Layer"
ann_colors <- list(Layer = layer_colors)

pheatmap::pheatmap(
  z_expr,
  annotation_row    = ann_row,
  annotation_colors = ann_colors,
  clustering_method = "ward.D2",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-3, 3, length.out = 101),
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 6,
  fontsize_col      = 9,
  main              = "RQ1: MAA Pan-Cancer Expression Atlas (z-score)",
  filename          = "figures/Fig1_pancancer_heatmap.pdf",
  width = 16, height = 20
)
cat("Fig1 saved.\n")

## ---------- 4. Co-expression Spearman matrix ----------
# Across all tumor samples (only genes present in matrix)
cat("Computing pan-cancer co-expression matrix...\n")
genes_for_cor <- intersect(found_genes, rownames(tub_tumor))
cor_mat <- cor(t(tub_tumor[genes_for_cor, , drop = FALSE]),
               method = "spearman", use = "pairwise.complete.obs")

# Replace any NA/NaN in correlation matrix (can occur for zero-variance genes)
cor_mat[!is.finite(cor_mat)] <- 0
diag(cor_mat) <- 1   # ensure diagonal is exactly 1

# Cluster genes by co-expression
hc <- hclust(as.dist(1 - cor_mat), method = "ward.D2")
gene_clusters <- cutree(hc, k = 5)  # 5 co-expression modules
module_df <- data.frame(
  gene   = names(gene_clusters),
  module = paste0("Module_", gene_clusters),
  stringsAsFactors = FALSE
)
module_df <- left_join(module_df, gene_layer_table, by = "gene")

# Co-expression heatmap (Fig2)
ann_row2 <- data.frame(
  Layer  = gene_layer_table[names(gene_clusters), "layer"],
  Module = paste0("M", gene_clusters),
  row.names = names(gene_clusters)
)
module_colors <- setNames(
  RColorBrewer::brewer.pal(5, "Set2"),
  paste0("M", 1:5)
)
ann_colors2 <- list(Layer = layer_colors, Module = module_colors)

pheatmap::pheatmap(
  cor_mat[hc$order, hc$order],
  annotation_row    = ann_row2,
  annotation_col    = ann_row2,
  annotation_colors = ann_colors2,
  color = colorRampPalette(c("#053061","#FFFFFF","#67001F"))(100),
  breaks            = seq(-1, 1, length.out = 101),
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 6,
  main              = "RQ1: MAA Co-expression Matrix (Spearman rho)",
  filename          = "figures/Fig2_coexpression_matrix.pdf",
  width = 14, height = 13
)
cat("Fig2 saved.\n")

## ---------- 5. Expression summary statistics ----------
expr_summary <- data.frame(
  gene        = found_genes,
  layer       = gene_layer_table[found_genes, "layer"],
  module      = module_df$module[match(found_genes, module_df$gene)],
  mean_expr   = rowMeans(tub_tumor, na.rm = TRUE),
  median_expr = rowMedians(tub_tumor),
  sd_expr     = apply(tub_tumor, 1, sd, na.rm = TRUE),
  cv_expr     = apply(tub_tumor, 1, sd, na.rm = TRUE) / rowMeans(tub_tumor, na.rm = TRUE),
  stringsAsFactors = FALSE
) %>%
  arrange(layer, desc(mean_expr))

## ---------- 6. Save ----------
wb <- createWorkbook()
addWorksheet(wb, "Expression_Summary");  writeData(wb, "Expression_Summary", expr_summary)
addWorksheet(wb, "Per_Cancer_Median");   writeData(wb, "Per_Cancer_Median",
                                                     cbind(gene=rownames(med_expr), med_expr))
addWorksheet(wb, "Coexpression_Modules"); writeData(wb, "Coexpression_Modules", module_df)
saveWorkbook(wb, "results/RQ1_expression_atlas.xlsx", overwrite = TRUE)
saveRDS(cor_mat,    "xena_data/RQ1_cor_mat.rds")
saveRDS(module_df,  "xena_data/RQ1_modules.rds")

cat("RQ1 complete. Outputs:\n")
cat("  results/RQ1_expression_atlas.xlsx\n")
cat("  figures/Fig1_pancancer_heatmap.pdf\n")
cat("  figures/Fig2_coexpression_matrix.pdf\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
