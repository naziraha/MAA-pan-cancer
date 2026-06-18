# ============================================================
# RQ5 — ADAPTIVE STATES & MAA EVOLUTIONARY TRAJECTORY
# Questions:
#   (a) What tumour archetypes emerge from MAA expression?
#   (b) Does MAA state predict resistance to therapy?
#   (c) Can we model the normal→cancer→metastatic→resistant trajectory?
# Outputs: Adaptive state map, survival by state, evolution model,
#          drug resistance state analysis, priority candidates
# ============================================================

source("R startup/00_gene_panel.R")
library(ggplot2); library(dplyr); library(openxlsx); library(tidyr)
library(survival); library(survminer); library(pheatmap)
library(RColorBrewer); library(ggrepel); library(cluster)

## ---------- 1. Load ----------
master       <- as.data.frame(readRDS("xena_data/master.rds"))  # ensure plain df
tub_tumor    <- readRDS("xena_data/tub_mat_tumor.rds")
trait_scores <- readRDS("xena_data/RQ3_trait_scores.rds")
multitrait   <- readRDS("xena_data/RQ3_multitrait.rds")
hub_genes    <- readRDS("xena_data/RQ4_hub_genes.rds")
surv_obj     <- readRDS("xena_data/RQ2_survival.rds")
stage_trend  <- readRDS("xena_data/RQ2_stage_trend.rds")
de_summary   <- readRDS("xena_data/RQ0_de_summary.rds")
found_genes  <- rownames(tub_tumor)

## ================================================================
## PART A: K-MEANS ADAPTIVE STATES
## ================================================================

## ---------- 2. MAA signature score (mean of top MAA genes) ----------
top_maa_genes <- head(multitrait$gene, 50)
maa_sig       <- colMeans(tub_tumor[intersect(top_maa_genes, found_genes), ,
                                     drop = FALSE], na.rm = TRUE)
master$maa_score <- maa_sig[master$barcode]

## ---------- 3. Feature matrix for clustering ----------
# Primary: MAA expression matrix (always available)
expr_part <- t(tub_tumor[found_genes, master$barcode, drop = FALSE])

# Secondary: trait scores — only include columns with sufficient non-NA data (>50%)
trait_match <- trait_scores[match(master$barcode, trait_scores$barcode), ]
trait_cols  <- grep("_score$|mRNAsi", names(trait_match), value = TRUE)
usable_trait_cols <- trait_cols[
  sapply(trait_cols, function(col) {
    vals <- trait_match[[col]]
    sum(!is.na(vals)) > 0.5 * nrow(trait_match)
  })
]
cat("Usable trait score columns for clustering:", length(usable_trait_cols), "\n")

if (length(usable_trait_cols) > 0) {
  trait_part <- as.matrix(trait_match[, usable_trait_cols, drop = FALSE])
  feat_mat   <- cbind(expr_part, trait_part)
} else {
  cat("No trait scores available — clustering on MAA expression only\n")
  feat_mat <- as.matrix(expr_part)
}

# Add CIN and stemness if available
if (sum(!is.na(master$cin_score)) > nrow(master) * 0.5)
  feat_mat <- cbind(feat_mat, cin_score = master$cin_score[match(rownames(feat_mat), master$barcode)])
if (sum(!is.na(master$stemness)) > nrow(master) * 0.5)
  feat_mat <- cbind(feat_mat, stemness = master$stemness[match(rownames(feat_mat), master$barcode)])

feat_mat   <- feat_mat[complete.cases(feat_mat), ]
cat("Clustering on", nrow(feat_mat), "samples x", ncol(feat_mat), "features\n")
feat_scaled <- scale(feat_mat)
# Replace any NaN from zero-variance columns
feat_scaled[!is.finite(feat_scaled)] <- 0

## ---------- 4. Determine optimal k ----------
set.seed(42)
wss <- sapply(2:8, function(k) {
  kmeans(feat_scaled, centers = k, nstart = 10, iter.max = 100)$tot.withinss
})
# Use k=5 for richer resolution (5 states: normal-like, early, proliferative,
#                                         mesenchymal/invasive, drug-resistant)
k_opt <- 5
km    <- kmeans(feat_scaled, centers = k_opt, nstart = 25, iter.max = 300)
master$adaptive_state <- NA
master[rownames(feat_mat), "adaptive_state"] <- paste0("State_", km$cluster)

cat("Adaptive state distribution:\n")
print(table(master$adaptive_state))

## ---------- 5. Characterise each state ----------
state_profiles <- master %>%
  filter(!is.na(adaptive_state)) %>%
  group_by(adaptive_state) %>%
  summarise(
    n = n(),
    mean_maa_score  = mean(maa_score, na.rm = TRUE),
    mean_cin        = mean(cin_score, na.rm = TRUE),
    mean_stemness   = mean(stemness,  na.rm = TRUE),
    pct_stage3_4    = mean(stage_num %in% 3:4, na.rm = TRUE) * 100,
    pct_OS_event    = mean(OS, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Add trait score means per state
trait_by_state <- master %>%
  filter(!is.na(adaptive_state)) %>%
  left_join(trait_scores[, c("barcode","EMT_score","Proliferation_score",
                               "Taxane_Resistance_score","Anti_Apoptosis_score",
                               "drug_resistance_score")],
            by = "barcode") %>%
  group_by(adaptive_state) %>%
  summarise(across(ends_with("_score"), ~mean(.x, na.rm=TRUE)), .groups="drop")

state_profiles <- left_join(state_profiles, trait_by_state, by = "adaptive_state")

# Label states by phenotype
state_labels <- c(
  State_1 = "Normal-like / Quiescent",
  State_2 = "Epithelial / Early",
  State_3 = "Proliferative / Stem-like",
  State_4 = "Mesenchymal / Invasive",
  State_5 = "Drug-Resistant / Aggressive"
)
# Assign labels based on mean MAA + EMT + resistance scores
score_rank <- state_profiles %>%
  mutate(
    emt_rank    = rank(EMT_score),
    resist_rank = rank(drug_resistance_score),
    prolif_rank = rank(Proliferation_score),
    maa_rank    = rank(mean_maa_score)
  ) %>%
  arrange(desc(resist_rank))

cat("\nState profiles:\n")
print(state_profiles)

## ---------- 6. Adaptive state heatmap ----------
# Sample 200 per state for display
set.seed(42)
samp_idx <- master %>%
  filter(!is.na(adaptive_state)) %>%
  group_by(adaptive_state) %>%
  slice_sample(n = 200, replace = FALSE) %>%
  pull(barcode)

heat_mat  <- t(tub_tumor[found_genes, samp_idx, drop = FALSE])
heat_mat  <- scale(heat_mat)
heat_mat  <- t(heat_mat)

state_pal <- setNames(RColorBrewer::brewer.pal(5,"Set1"), paste0("State_",1:5))

# Annotation: use only State + Cancer (drop stage_num — numeric with NAs breaks pheatmap)
ann_col <- data.frame(
  State  = master$adaptive_state[match(samp_idx, master$barcode)],
  Cancer = master$cancer_type[match(samp_idx, master$barcode)],
  row.names = samp_idx,
  stringsAsFactors = FALSE
)
ann_col <- ann_col[!is.na(ann_col$State), , drop = FALSE]
samp_idx <- rownames(ann_col)   # keep only matched rows

ann_row <- gene_layer_table[found_genes, "layer", drop = FALSE]
colnames(ann_row) <- "Layer"

pheatmap::pheatmap(
  pmin(pmax(heat_mat[, samp_idx, drop=FALSE], -3), 3),
  annotation_col    = ann_col,
  annotation_row    = ann_row,
  annotation_colors = list(State = state_pal, Layer = layer_colors),
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 5,
  clustering_method = "ward.D2",
  color             = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
  main              = "RQ5: MAA Adaptive State Heatmap",
  filename          = "figures/RQ5_adaptive_state_heatmap.pdf",
  width = 18, height = 20
)

## ---------- 7. Survival by adaptive state ----------
df_surv <- master[!is.na(master$adaptive_state) &
                    !is.na(master$OS) &
                    !is.na(master$OS.time) &
                    master$OS.time > 0, ]
df_surv$state <- df_surv$adaptive_state   # clean column, no formula parsing issues

# Build KM curves manually with broom::tidy for ggplot compatibility
km_list <- lapply(sort(unique(df_surv$state)), function(s) {
  d   <- df_surv[df_surv$state == s, ]
  fit <- survfit(Surv(OS.time / 365.25, OS) ~ 1, data = d)
  data.frame(
    time  = fit$time,
    surv  = fit$surv,
    lower = fit$lower,
    upper = fit$upper,
    state = s,
    stringsAsFactors = FALSE
  )
})
km_df <- do.call(rbind, km_list)
state_pal6 <- setNames(RColorBrewer::brewer.pal(max(5, length(unique(km_df$state))), "Set1"),
                        sort(unique(km_df$state)))

p_surv <- ggplot(km_df, aes(x = time, y = surv, colour = state, fill = state)) +
  geom_step(size = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
  scale_colour_manual(values = state_pal6) +
  scale_fill_manual(values   = state_pal6) +
  coord_cartesian(xlim = c(0, 15)) +
  labs(title = "RQ5: Overall Survival by MAA Adaptive State",
       x = "Time (years)", y = "Survival probability",
       colour = "State", fill = "State") +
  theme_bw(base_size = 12)
ggsave("figures/RQ5_survival_by_state.pdf", p_surv, width = 10, height = 7)
cat("Survival plot saved.\n")

## ================================================================
## PART B: EVOLUTIONARY TRAJECTORY MODEL
## ================================================================

## ---------- 8. MAA trajectory score (pseudo-progression index) ----------
# Integrate: normal→tumor FC (RQ0), stage trend (RQ2), fitness breadth (RQ3),
#            network centrality (RQ4), and survival HR (RQ2)
# Each component normalised 0–1, then summed

norm_01 <- function(x) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

# Helper: safely extract a numeric column from a data frame by gene matching
pull_num <- function(df, gene_col, val_col, genes) {
  as.numeric(df[[val_col]][match(genes, df[[gene_col]])])
}

# Component 1: upregulation in tumor vs normal (initiation)
comp1 <- pull_num(as.data.frame(de_summary), "gene", "mean_logFC", found_genes)
comp1[is.na(comp1)] <- 0
comp1 <- norm_01(pmax(comp1, 0))

# Component 2: stage-linear increase (progression)
comp2 <- pull_num(stage_trend, "gene", "beta", found_genes)
comp2[is.na(comp2)] <- 0
comp2 <- norm_01(pmax(comp2, 0))

# Component 3: multi-trait fitness breadth (adaptability)
comp3 <- pull_num(as.data.frame(multitrait), "gene", "n_traits_sig", found_genes)
comp3[is.na(comp3)] <- 0
comp3 <- norm_01(comp3)

# Component 4: network hub centrality
comp4 <- pull_num(hub_genes, "gene", "hub_score", found_genes)
comp4[is.na(comp4)] <- 0
comp4 <- norm_01(comp4)

# Component 5: survival hazard ratio
comp5 <- pull_num(as.data.frame(surv_obj), "gene", "HR", found_genes)
comp5[is.na(comp5)] <- 1
comp5 <- norm_01(pmax(log(pmax(comp5, 1e-6)), 0))

# Component 6: drug resistance coupling
mt_df    <- as.data.frame(readRDS("xena_data/RQ3_multitrait.rds"))
drug_rho <- pull_num(mt_df, "gene", "max_resist_rho", found_genes)
drug_rho[is.na(drug_rho) | !is.finite(drug_rho)] <- 0
comp6 <- norm_01(drug_rho)

# Layer vector — safe plain-vector extraction
layer_vec <- as.character(gene_layer_table$layer[match(found_genes, gene_layer_table$gene)])

# Composite MAA Evolutionary Index (MEI)
mei_df <- data.frame(
  gene        = found_genes,
  layer       = layer_vec,
  C1_initiation   = comp1,
  C2_progression  = comp2,
  C3_fitness_breadth = comp3,
  C4_network_hub  = comp4,
  C5_survival_HR  = comp5,
  C6_drug_resistance = comp6,
  MEI_score = (comp1 + comp2 + comp3 + comp4 + comp5 + comp6) / 6,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(MEI_score))

cat("\nTop 20 MAA genes by Evolutionary Index (MEI):\n")
print(head(mei_df[, c("gene","layer","MEI_score","C1_initiation","C2_progression",
                        "C3_fitness_breadth","C6_drug_resistance")], 20))

## ---------- 9. MEI bubble plot (trajectory visualisation) ----------
mei_df$label <- ifelse(mei_df$MEI_score >= quantile(mei_df$MEI_score, 0.85),
                         mei_df$gene, "")
p_mei <- ggplot(mei_df, aes(x = C2_progression, y = C6_drug_resistance,
                              size = MEI_score, colour = layer, label = label)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_colour_manual(values = layer_colors, na.value = "grey60") +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "RQ5: MAA Evolutionary Trajectory Map",
       subtitle = "x = Stage-progression driver | y = Drug resistance coupling",
       x = "Stage-Progression Score (C2)", y = "Drug Resistance Score (C6)",
       size = "MEI Score") +
  theme_bw(base_size = 12)
ggsave("figures/RQ5_MAA_evolutionary_trajectory.pdf", p_mei, width = 11, height = 9)

## ---------- 10. PCA fitness landscape ----------
master_df  <- as.data.frame(master)   # ensure plain data.frame (not tibble)
state_mask <- !is.na(master_df$adaptive_state)
pca_feats  <- master_df[state_mask, intersect(found_genes, names(master_df)),
                          drop = FALSE]
pca_feats  <- pca_feats[complete.cases(pca_feats), ]
pca_res    <- prcomp(pca_feats, scale. = TRUE, center = TRUE)

# Use match() to extract metadata — avoids tibble matrix-subscript bug
row_idx <- match(rownames(pca_feats), master_df$barcode)
pca_df  <- data.frame(
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  state  = master_df$adaptive_state[row_idx],
  cancer = master_df$cancer_type[row_idx],
  stage  = master_df$stage_num[row_idx],
  stringsAsFactors = FALSE
)
var_exp   <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = state)) +
  geom_point(size = 0.8, alpha = 0.5) +
  stat_ellipse(level = 0.90, size = 1) +
  scale_colour_manual(values = state_pal) +
  labs(title = "RQ5: MAA Fitness Landscape (PCA)",
       x = paste0("PC1 (", var_exp[1], "% var)"),
       y = paste0("PC2 (", var_exp[2], "% var)")) +
  theme_bw(base_size = 12)
ggsave("figures/RQ5_PCA_fitness_landscape.pdf", p_pca, width = 10, height = 8)

# PCA coloured by cancer type
p_pca_ct <- p_pca + aes(colour = cancer) +
  scale_colour_discrete() +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 7))
ggsave("figures/RQ5_PCA_by_cancer_type.pdf", p_pca_ct, width = 14, height = 9)

## ================================================================
## PART C: PRIORITY CANDIDATES
## ================================================================

## ---------- 11. Priority MAA candidates (integrated ranking) ----------
priority_candidates <- mei_df %>%
  left_join(hub_genes[, c("gene","degree","betweenness")], by = "gene") %>%
  left_join(multitrait[, c("gene","n_traits_sig","n_resist_traits")], by = "gene") %>%
  left_join(surv_obj[, c("gene","HR","q")], by = "gene") %>%
  arrange(desc(MEI_score)) %>%
  select(gene, layer, MEI_score, C1_initiation, C2_progression,
         C3_fitness_breadth, C4_network_hub, C5_survival_HR,
         C6_drug_resistance, degree, betweenness, n_traits_sig,
         n_resist_traits, HR, q)

cat("\nTop 15 MAA Priority Candidates:\n")
print(head(priority_candidates[, c("gene","layer","MEI_score",
                                     "C6_drug_resistance","n_resist_traits")], 15))

## ---------- 12. Save ----------
wb <- createWorkbook()
addWorksheet(wb, "MEI_Rankings");        writeData(wb,"MEI_Rankings",       mei_df)
addWorksheet(wb, "Priority_Candidates"); writeData(wb,"Priority_Candidates",priority_candidates)
addWorksheet(wb, "State_Profiles");      writeData(wb,"State_Profiles",     state_profiles)
addWorksheet(wb, "Adaptive_States");
  state_assign <- master[!is.na(master$adaptive_state),
                          c("barcode","cancer_type","adaptive_state","maa_score",
                            "cin_score","stemness","stage_num")]
  writeData(wb,"Adaptive_States", state_assign)
saveWorkbook(wb, "results/RQ5_adaptive_states_evolution.xlsx", overwrite = TRUE)
saveRDS(mei_df,            "xena_data/RQ5_mei.rds")
saveRDS(priority_candidates,"xena_data/RQ5_priority_candidates.rds")
saveRDS(master,            "xena_data/master_final.rds")  # master with states

cat("RQ5 complete. Outputs:\n")
cat("  results/RQ5_adaptive_states_evolution.xlsx\n")
cat("  figures/RQ5_adaptive_state_heatmap.pdf\n")
cat("  figures/RQ5_survival_by_state.pdf\n")
cat("  figures/RQ5_MAA_evolutionary_trajectory.pdf\n")
cat("  figures/RQ5_PCA_fitness_landscape.pdf\n")
cat("  figures/RQ5_PCA_by_cancer_type.pdf\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
