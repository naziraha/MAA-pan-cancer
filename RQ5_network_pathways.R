# ============================================================
# RQ4 — NETWORK ARCHITECTURE, PATHWAY ENRICHMENT & DRUG TARGETS
# Questions:
#   (a) What is the network topology of the MAA system?
#   (b) What oncogenic pathways do MAA hub genes connect?
#   (c) Which MAA genes have druggable targets (approved/clinical)?
# Outputs: Network graph, enrichment tables, drug-target map
# ============================================================

source("R startup/00_gene_panel.R")
library(igraph); library(ggraph); library(ggplot2); library(dplyr)
library(clusterProfiler); library(org.Hs.eg.db); library(openxlsx)
library(enrichplot); library(msigdbr); library(tidyr)

## ---------- 1. Load ----------
master     <- readRDS("xena_data/master.rds")
tub_tumor  <- readRDS("xena_data/tub_mat_tumor.rds")
cor_mat    <- readRDS("xena_data/RQ1_cor_mat.rds")
found_genes <- rownames(tub_tumor)

## ================================================================
## PART A: CO-EXPRESSION NETWORK
## ================================================================

## ---------- 2. Build network (|rho| > 0.3) ----------
adj     <- abs(cor_mat) > 0.3 & upper.tri(cor_mat)
diag(adj) <- FALSE
edges   <- which(adj, arr.ind = TRUE)
edge_df <- data.frame(
  from  = rownames(cor_mat)[edges[,1]],
  to    = colnames(cor_mat)[edges[,2]],
  weight= cor_mat[edges],
  sign  = ifelse(cor_mat[edges] > 0, "positive", "negative"),
  stringsAsFactors = FALSE
)

g <- igraph::graph_from_data_frame(edge_df, directed = FALSE,
                                     vertices = found_genes)
V(g)$layer  <- gene_layer_table[V(g)$name, "layer"]
V(g)$degree <- igraph::degree(g)
# betweenness requires positive weights — use absolute correlation values
E(g)$abs_weight <- abs(E(g)$weight)
V(g)$betweenness <- igraph::betweenness(g, weights = E(g)$abs_weight,
                                         normalized = TRUE)

cat("Network: ", vcount(g), "nodes,", ecount(g), "edges\n")

## ---------- 3. Hub genes ----------
hub_genes <- data.frame(
  gene        = V(g)$name,
  layer       = V(g)$layer,
  degree      = V(g)$degree,
  betweenness = V(g)$betweenness,
  stringsAsFactors = FALSE
) %>%
  mutate(hub_score = scale(degree)[,1] + scale(betweenness)[,1]) %>%
  arrange(desc(hub_score))

cat("Top 15 network hub genes:\n")
print(head(hub_genes[, c("gene","layer","degree","betweenness","hub_score")], 15))

## ---------- 4. Network visualisation ----------
# Focus on largest connected component
components <- igraph::components(g)
main_comp  <- induced_subgraph(g, which(components$membership ==
                                          which.max(components$csize)))

# FR layout requires positive weights — set abs_weight attribute on subgraph
E(main_comp)$abs_weight <- abs(E(main_comp)$weight)
set.seed(42)
p_network <- ggraph(main_comp, layout = "fr", weights = abs_weight) +
  geom_edge_link(aes(alpha = abs(weight), colour = sign), width = 0.4) +
  geom_node_point(aes(size = degree, colour = layer), alpha = 0.9) +
  geom_node_text(aes(label = ifelse(degree >= quantile(degree, 0.8), name, "")),
                  size = 2.5, repel = TRUE) +
  scale_colour_manual(values = c(layer_colors,
                                   positive = "#999999", negative = "#CCCCCC")) +
  scale_edge_colour_manual(values = c(positive = "#AAAAAA", negative = "#CC6600")) +
  scale_size_continuous(range = c(1, 9)) +
  labs(title = "RQ4: MAA Co-expression Network (|rho| > 0.3)") +
  theme_void() + theme(legend.position = "right")
ggsave("figures/RQ4_coexpression_network.pdf",  p_network, width = 14, height = 12)
ggsave("figures/RQ4_coexpression_network.png",  p_network, width = 14, height = 12, dpi = 150)

## ================================================================
## PART B: PATHWAY ENRICHMENT
## ================================================================

## ---------- 5. Convert to Entrez IDs ----------
entrez <- bitr(found_genes, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
entrez_ids <- entrez$ENTREZID

## ---------- 6. KEGG enrichment ----------
cat("KEGG enrichment...\n")
kegg_res <- enrichKEGG(gene = entrez_ids, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
kegg_df  <- as.data.frame(kegg_res@result)

## ---------- 7. Reactome enrichment ----------
if (requireNamespace("ReactomePA", quietly = TRUE)) {
  library(ReactomePA)
  react_res <- enrichPathway(gene = entrez_ids, organism = "human",
                               pAdjustMethod = "BH", pvalueCutoff = 0.05)
  react_df  <- as.data.frame(react_res@result)
} else {
  react_df <- data.frame(note = "Install ReactomePA for Reactome enrichment")
}

## ---------- 8. GO enrichment (BP + MF) ----------
cat("GO enrichment...\n")
go_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
go_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
go_bp_df <- as.data.frame(go_bp@result)
go_mf_df <- as.data.frame(go_mf@result)

## ---------- 9. Cancer hallmark enrichment (MSigDB) ----------
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")
hallmarks_list <- split(hallmarks$gene_symbol, hallmarks$gs_name)
hall_res <- clusterProfiler::enricher(
  gene   = found_genes,
  TERM2GENE = hallmarks[, c("gs_name","gene_symbol")],
  pAdjustMethod = "BH", pvalueCutoff = 0.05
)
hall_df <- as.data.frame(hall_res@result)

## ---------- 10. Publication-quality bar charts ----------
# Shared bar chart function: horizontal bars, sorted by -log10(padj), coloured by gene ratio
enrich_barchart <- function(df, title, top_n = 20, padj_col = "p.adjust",
                              ratio_col = "GeneRatio", desc_col = "Description") {
  df <- df[df[[padj_col]] < 0.05 & !is.na(df[[padj_col]]), ]
  if (nrow(df) == 0) { message("No significant terms for: ", title); return(NULL) }
  df <- head(df[order(df[[padj_col]]), ], top_n)
  # Parse GeneRatio (e.g. "12/132") to numeric
  if (is.character(df[[ratio_col]])) {
    df$ratio_num <- sapply(strsplit(df[[ratio_col]], "/"),
                            function(x) as.numeric(x[1]) / as.numeric(x[2]))
  } else {
    df$ratio_num <- as.numeric(df[[ratio_col]])
  }
  df$log10padj <- -log10(df[[padj_col]])
  # Shorten long pathway names
  df[[desc_col]] <- gsub("HALLMARK_", "", df[[desc_col]])
  df[[desc_col]] <- gsub("_", " ", df[[desc_col]])
  df[[desc_col]] <- stringr::str_wrap(df[[desc_col]], width = 45)
  df[[desc_col]] <- factor(df[[desc_col]], levels = rev(df[[desc_col]]))

  ggplot(df, aes(x = log10padj, y = .data[[desc_col]], fill = ratio_num)) +
    geom_bar(stat = "identity", width = 0.75) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    scale_fill_gradient(low = "#FEE0D2", high = "#CB181D",
                         name = "Gene\nratio") +
    labs(title = title,
         x = expression(-log[10](adjusted~italic(p))),
         y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.y  = element_text(size = 9),
          plot.title   = element_text(face = "bold", size = 12),
          panel.grid.major.y = element_blank())
}

# KEGG bar chart
p_kegg_bar <- enrich_barchart(kegg_df,
                               title = "KEGG Pathway Enrichment — MAA Genes",
                               top_n = 20)
if (!is.null(p_kegg_bar))
  ggsave("figures/RQ4_KEGG_barchart.pdf", p_kegg_bar, width = 9, height = 8)

# Hallmark bar chart
p_hall_bar <- enrich_barchart(hall_df,
                               title = "MSigDB Hallmark Enrichment — MAA Genes",
                               top_n = 20)
if (!is.null(p_hall_bar))
  ggsave("figures/RQ4_Hallmark_barchart.pdf", p_hall_bar, width = 9, height = 8)

# GO-BP bar chart (bonus)
p_go_bar <- enrich_barchart(go_bp_df,
                              title = "GO Biological Process Enrichment — MAA Genes",
                              top_n = 20)
if (!is.null(p_go_bar))
  ggsave("figures/RQ4_GOBP_barchart.pdf", p_go_bar, width = 9, height = 9)

cat("Bar charts saved: KEGG, Hallmark, GO-BP\n")

## ================================================================
## PART C: DRUG TARGET ANALYSIS
## ================================================================

## ---------- 11. Known drug targets & clinical drugs ----------
# Curated list of MAA genes with existing drugs (approved / Phase I-III)
drug_targets <- data.frame(
  gene = c("CDK1","HDAC6","HDAC1","HDAC2","HDAC8","SIRT2",
            "AURKA","AURKB","KIF11","CENPE","KIF18A","KIF15",
            "KIFC1","KIF20A","KIF2C","BIRC5","BRCA1","PLK1",
            "KIF5B","TUBB3","TUBB","TUBB4A","VASH1","VASH2",
            "SETD2","KAT2A","PRKAA1","NF1","TTL","ATAT1"),
  drug_class = c("CDK inhibitor","HDAC inhibitor","HDAC inhibitor","HDAC inhibitor",
                  "HDAC inhibitor","SIRT inhibitor",
                  "Aurora kinase inhibitor","Aurora kinase inhibitor",
                  "Kinesin Eg5 inhibitor","CENP-E inhibitor",
                  "Kinesin KIF18A inhibitor","Kinesin KIF15 inhibitor",
                  "HSET/KIFC1 inhibitor","KIF20A inhibitor","MCAK/KIF2C modulator",
                  "IAP inhibitor / Survivin inhibitor","PARP inhibitor target",
                  "PLK1 inhibitor (indirect)","Microtubule stabilizer context",
                  "Taxane resistance marker","Taxane target","Taxane resistance marker",
                  "Vasohibin inhibitor","Vasohibin inhibitor",
                  "EZH2/methylation context","KAT2A/GCN5 inhibitor",
                  "AMPK activator","RAS pathway inhibitor","TTL restoration therapy",
                  "Alpha-tubulin acetyltransferase inhibitor"),
  example_drugs = c("Dinaciclib,AT7519","Vorinostat,Tubastatin A","Entinostat,Romidepsin",
                     "Entinostat","PCI-34051","Sirtinol",
                     "Alisertib,MK-5108","AZD1152,Barasertib",
                     "Ispinesib,SB-743921","GSK923295",
                     "AMG-900 context","Kinesin inhibitors in development",
                     "CW069","CW069 analogs","In development",
                     "YM155,LY2181308","Olaparib,Rucaparib (BRCA1 mut)",
                     "Volasertib (via pathway)","Paclitaxel,Docetaxel resistance",
                     "Paclitaxel,Docetaxel","Paclitaxel,Docetaxel","Cabazitaxel resistance",
                     "In development","In development",
                     "Tazemetostat context","MB-3,CPTH2",
                     "Metformin,AICAR","Trametinib context","Investigational",
                     "Tubastatin A context"),
  development_stage = c("Approved/Clinical","Approved","Approved","Clinical","Clinical","Preclinical",
                          "Clinical","Approved","Clinical","Clinical","Preclinical","Preclinical",
                          "Preclinical","Preclinical","Preclinical","Clinical","Approved",
                          "Clinical","Target context","Biomarker","Target","Biomarker",
                          "Preclinical","Preclinical","Context","Preclinical",
                          "Approved(Metformin)","Context","Investigational","Preclinical"),
  stringsAsFactors = FALSE
)
drug_targets <- left_join(drug_targets, gene_layer_table, by = "gene")

## ---------- 12. Convergent genes (top RQ2 ∩ top RQ3) ----------
surv_obj   <- readRDS("xena_data/RQ2_survival.rds")
multitrait <- readRDS("xena_data/RQ3_multitrait.rds")

top_surv   <- head(surv_obj$gene[surv_obj$q < 0.05], 20)
top_multi  <- head(multitrait$gene, 20)
convergent <- intersect(top_surv, top_multi)
cat("Convergent genes (top RQ2 ∩ top RQ3):", paste(convergent, collapse = ", "), "\n")

## ---------- 13. Save ----------
wb <- createWorkbook()
addWorksheet(wb, "Hub_Genes");     writeData(wb, "Hub_Genes",     hub_genes)
addWorksheet(wb, "KEGG");          writeData(wb, "KEGG",          kegg_df)
addWorksheet(wb, "Reactome");      writeData(wb, "Reactome",      react_df)
addWorksheet(wb, "GO_BP");         writeData(wb, "GO_BP",         go_bp_df)
addWorksheet(wb, "GO_MF");         writeData(wb, "GO_MF",         go_mf_df)
addWorksheet(wb, "Hallmarks");     writeData(wb, "Hallmarks",     hall_df)
addWorksheet(wb, "Drug_Targets");  writeData(wb, "Drug_Targets",  drug_targets)
addWorksheet(wb, "Convergent_Genes"); writeData(wb, "Convergent_Genes",
                                                 data.frame(gene = convergent))
saveWorkbook(wb, "results/RQ4_network_pathways.xlsx", overwrite = TRUE)
saveRDS(hub_genes,    "xena_data/RQ4_hub_genes.rds")
saveRDS(drug_targets, "xena_data/RQ4_drug_targets.rds")
write.csv(edge_df, "results/Cytoscape_edge_table.csv", row.names = FALSE, quote = FALSE)
write.csv(hub_genes,"results/Cytoscape_node_table.csv",row.names = FALSE, quote = FALSE)

cat("RQ4 complete. Outputs:\n")
cat("  results/RQ4_network_pathways.xlsx\n")
cat("  figures/RQ4_coexpression_network.pdf/.png\n")
cat("  figures/RQ4_KEGG_dotplot.pdf\n")
cat("  figures/RQ4_Hallmark_dotplot.pdf\n")
cat("  results/Cytoscape_edge_table.csv (updated)\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
