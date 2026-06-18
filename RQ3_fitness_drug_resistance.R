# ============================================================
# RQ3 — FITNESS TRAIT COUPLING + DRUG RESISTANCE ANALYSIS
# Questions:
#   (a) Is MAA comprehensively coupled to oncogenic fitness traits?
#   (b) Which MAA genes correlate with drug-resistance signatures?
#   (c) Which tumours use MAA as a multitrait fitness strategy?
# Fitness traits: EMT, Proliferation, Stemness, Immune infiltration,
#                 Immune evasion, Mechanical/invasion, Focal adhesion,
#                 Drug resistance (taxane + MDR + anti-apoptosis)
# ============================================================

source("R startup/00_gene_panel.R")
library(ggplot2); library(dplyr); library(openxlsx); library(tidyr)
library(pheatmap); library(ggrepel); library(RColorBrewer)

## ---------- 1. Load ----------
master    <- readRDS("xena_data/master.rds")
tub_tumor <- readRDS("xena_data/tub_mat_tumor.rds")
found_genes <- rownames(tub_tumor)

## ---------- 2. Fitness trait gene sets ----------
# Each score = mean expression of the gene set members found in data

trait_genesets <- list(
  # A. Epithelial-Mesenchymal Transition
  EMT = list(
    mesenchymal = c("VIM","FN1","CDH2","TWIST1","TWIST2","SNAI1","SNAI2",
                    "ZEB1","ZEB2","MMP2","MMP9","ACTA2","S100A4","SPARC"),
    epithelial  = c("CDH1","EPCAM","KRT19","KRT18","KRT8","OCLN","CLDN4")
  ),
  # B. Proliferation
  Proliferation = c("MKI67","CDK1","PCNA","TOP2A","CCNB1","CCNE1","AURKA",
                     "AURKB","BUB1","BIRC5","CDC20","PLK1","MCM2","MCM7"),
  # C. Stemness — supplemented by mRNAsi Xena scores
  Stemness      = c("SOX2","OCT4","NANOG","KLF4","MYC","ALDH1A1","CD44",
                     "PROM1","ABCG2","NOTCH1","WNT5A","LGR5"),
  # D. Immune Infiltration
  Immune_Infiltration = c("CD3E","CD4","CD8A","CD68","FOXP3","NCAM1",
                            "CD19","MS4A1","GZMB","PRF1","IFNG"),
  # E. Immune Evasion (checkpoint)
  Immune_Evasion = c("CD274","PDCD1LG2","CTLA4","PDCD1","LAG3","HAVCR2",
                      "TIGIT","IDO1","TGFB1","IL10","VEGFA"),
  # F. Mechanical / Invasion
  Invasion = c("MMP1","MMP2","MMP3","MMP9","MMP14","ITGB1","ITGAV",
                "RHOA","ROCK1","ROCK2","CTTN","WAVE2","CDC42","RAC1",
                "PAK1","LAMC1","COL1A1","FBN1"),
  # G. Focal Adhesion / ECM
  Focal_Adhesion = c("FAK1","PTK2","SRC","PAXILLIN","VINCULIN","ILK",
                      "ITGA5","ITGB5","LAMA4","FN1","TNC","THBS1","VCL"),
  # H. Drug Resistance — EXPANDED (key innovation vs old analysis)
  # H1: Taxane resistance (microtubule-targeting drug resistance)
  Taxane_Resistance = c("TUBB3","TUBB4A","MAP4","STMN1","STMN2","STMN3","STMN4",
                          "PTPN11","ABCB1","ABCC1","ABCG2","MAP2","BCL2"),
  # H2: Multi-drug resistance (MDR / efflux pumps)
  MDR_Efflux = c("ABCB1","ABCC1","ABCC2","ABCC3","ABCG2","MVP","LRP",
                   "GSTP1","YBX1","ABCA1"),
  # H3: Anti-apoptosis / survival
  Anti_Apoptosis = c("BCL2","BCL2L1","MCL1","BIRC5","BIRC2","BIRC3",
                      "XIAP","BCL2A1","CFLAR","TRAF2","MAP3K14"),
  # H4: Epithelial-to-mesenchymal plasticity drug resistance
  EMT_Drug_Resistance = c("SNAI1","ZEB1","VIM","TWIST1","NOTCH1",
                            "AXL","GAS6","EGFR","ERBB2","ALDH1A1")
)

## ---------- 3. Calculate trait scores for each sample ----------
calc_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) < 2) return(rep(NA, ncol(mat)))
  colMeans(mat[g, , drop = FALSE], na.rm = TRUE)
}

trait_scores <- data.frame(
  barcode    = master$barcode,
  cancer_type = master$cancer_type,
  stringsAsFactors = FALSE
)

# EMT: mesenchymal − epithelial
emt_mes <- calc_score(tub_tumor, trait_genesets$EMT$mesenchymal)
emt_epi <- calc_score(tub_tumor, trait_genesets$EMT$epithelial)
trait_scores$EMT_score   <- emt_mes - emt_epi

# Simple mean for all other traits
simple_traits <- c("Proliferation","Stemness","Immune_Infiltration","Immune_Evasion",
                    "Invasion","Focal_Adhesion","Taxane_Resistance","MDR_Efflux",
                    "Anti_Apoptosis","EMT_Drug_Resistance")
for (tr in simple_traits) {
  trait_scores[[paste0(tr, "_score")]] <- calc_score(tub_tumor, trait_genesets[[tr]])
}

# Add stemness from Xena (mRNAsi)
if ("stemness" %in% names(master)) {
  trait_scores$mRNAsi <- master$stemness[match(trait_scores$barcode, master$barcode)]
}

## ---------- 4. Spearman correlations: each MAA gene vs each trait ----------
score_cols <- grep("_score$|mRNAsi", names(trait_scores), value = TRUE)
cat("Computing Spearman correlations for", length(found_genes), "genes x",
    length(score_cols), "traits...\n")

cor_results <- lapply(found_genes, function(g) {
  gexpr <- master[[g]]
  lapply(score_cols, function(tr) {
    ok <- !is.na(gexpr) & !is.na(trait_scores[[tr]])
    if (sum(ok) < 100) return(NULL)
    ct <- cor.test(gexpr[ok], trait_scores[[tr]][ok], method = "spearman")
    data.frame(gene = g, trait = tr, rho = ct$estimate, p = ct$p.value,
               stringsAsFactors = FALSE)
  })
})
cor_df <- do.call(rbind, do.call(c, cor_results))
cor_df  <- cor_df[!is.null(cor_df), ]
cor_df$q <- p.adjust(cor_df$p, method = "BH")
cor_df <- left_join(cor_df, gene_layer_table, by = "gene")

pct_sig <- mean(cor_df$q < 0.05, na.rm = TRUE) * 100
cat(sprintf("%.1f%% of MAA gene-trait correlations significant (BH q<0.05)\n", pct_sig))

## ---------- 5. Multi-trait gene ranking ----------
multitrait <- cor_df %>%
  filter(q < 0.05) %>%
  group_by(gene) %>%
  summarise(
    n_traits_sig    = n_distinct(trait),
    mean_abs_rho    = mean(abs(rho), na.rm = TRUE),
    # Drug resistance specifically
    n_resist_traits = sum(grepl("Taxane|MDR|Anti_Apop|EMT_Drug", trait)),
    max_resist_rho  = max(abs(rho[grepl("Taxane|MDR|Anti_Apop|EMT_Drug", trait)]),
                           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(gene_layer_table, by = "gene") %>%
  arrange(desc(n_traits_sig), desc(mean_abs_rho))

cat("Top multi-trait MAA genes:\n")
print(head(select(multitrait, gene, layer, n_traits_sig, mean_abs_rho,
                   n_resist_traits, max_resist_rho), 15))

## ---------- 6. Drug resistance score per sample ----------
resist_genes   <- c("Taxane_Resistance_score","MDR_Efflux_score",
                     "Anti_Apoptosis_score","EMT_Drug_Resistance_score")
trait_scores$drug_resistance_score <- rowMeans(
  trait_scores[, resist_genes], na.rm = TRUE)

## ---------- 7. Figure: heatmap of MAA–trait correlations ----------
rho_mat <- cor_df %>%
  select(gene, trait, rho) %>%
  pivot_wider(names_from = trait, values_from = rho) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Order genes by mean absolute rho
gene_order <- names(sort(rowMeans(abs(rho_mat), na.rm = TRUE), decreasing = TRUE))
rho_mat    <- rho_mat[gene_order, , drop = FALSE]

ann_row <- gene_layer_table[rownames(rho_mat), "layer", drop = FALSE]
colnames(ann_row) <- "Layer"

pheatmap::pheatmap(
  rho_mat,
  annotation_row    = ann_row,
  annotation_colors = list(Layer = layer_colors),
  color             = colorRampPalette(c("#053061","#FFFFFF","#67001F"))(100),
  breaks            = seq(-0.7, 0.7, length.out = 101),
  clustering_method = "ward.D2",
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  fontsize_row      = 6,
  fontsize_col      = 9,
  main              = "RQ3: MAA–Fitness Trait Correlation Matrix (Spearman rho)",
  filename          = "figures/Fig5_fitness_trait_correlations.pdf",
  width = 14, height = 18
)

## ---------- 8. Drug resistance bubble plot ----------
drug_summary <- cor_df %>%
  filter(grepl("Taxane|MDR|Anti_Apop|EMT_Drug", trait)) %>%
  group_by(gene, layer) %>%
  summarise(mean_rho = mean(rho, na.rm = TRUE),
            n_sig    = sum(q < 0.05, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(mean_rho)) %>%
  arrange(desc(mean_rho))

top_resist <- head(drug_summary$gene, 20)
drug_plot_df <- cor_df %>%
  filter(gene %in% top_resist, grepl("Taxane|MDR|Anti_Apop|EMT_Drug", trait)) %>%
  mutate(trait = gsub("_score","",trait))

p_drug <- ggplot(drug_plot_df, aes(x = trait, y = reorder(gene, rho),
                                    size = abs(rho), colour = rho)) +
  geom_point(alpha = 0.85) +
  scale_colour_gradient2(low = "#053061", mid = "white", high = "#67001F",
                          midpoint = 0, limits = c(-0.7, 0.7)) +
  scale_size_continuous(range = c(1, 8)) +
  labs(title = "RQ3: MAA Genes — Drug Resistance Trait Correlations",
       x = NULL, y = NULL, colour = "Spearman rho", size = "|rho|") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("figures/RQ3_drug_resistance_correlations.pdf", p_drug, width = 9, height = 10)

## ---------- 9. Save ----------
wb <- createWorkbook()
addWorksheet(wb, "Trait_Correlations_All"); writeData(wb,"Trait_Correlations_All",cor_df)
addWorksheet(wb, "Multitrait_Ranking");     writeData(wb,"Multitrait_Ranking",multitrait)
addWorksheet(wb, "Drug_Resistance_Summary");writeData(wb,"Drug_Resistance_Summary",drug_summary)
addWorksheet(wb, "Trait_Scores_Samples");   writeData(wb,"Trait_Scores_Samples",trait_scores)
saveWorkbook(wb, "results/RQ3_fitness_drug_resistance.xlsx", overwrite = TRUE)
saveRDS(trait_scores, "xena_data/RQ3_trait_scores.rds")
saveRDS(multitrait,   "xena_data/RQ3_multitrait.rds")
saveRDS(cor_df,       "xena_data/RQ3_cor_df.rds")

cat("RQ3 complete. Outputs:\n")
cat("  results/RQ3_fitness_drug_resistance.xlsx\n")
cat("  figures/Fig5_fitness_trait_correlations.pdf\n")
cat("  figures/RQ3_drug_resistance_correlations.pdf\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
