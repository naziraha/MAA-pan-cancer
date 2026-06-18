# ============================================================
# MAA PAN-CANCER ANALYSIS — MASTER TABLE CONSTRUCTION v2.0
# Merges: TPM expression + clinical + GISTIC2 + MC3 + stemness + subtype
# Separates: tumor vs normal samples
# Output: master.rds, tub_mat_tumor.rds, tub_mat_normal.rds
# ============================================================
# Prereq: Run 00_gene_panel.R and 01_download_data.R first.

source("R startup/00_gene_panel.R")

## ---------- 1. Load full TPM matrix ----------
cat("Loading TPM matrix...\n")
rna_tpm <- readRDS("xena_data/rna_tpm_pancan.rds")
cat("Matrix:", nrow(rna_tpm), "genes x", ncol(rna_tpm), "samples\n")

## ---------- 2. Identify tumor vs normal samples ----------
# TCGA barcode: TCGA-XX-XXXX-[sample_type][vial]
# Sample types: 01=Primary tumor, 06=Metastatic, 10=Blood normal, 11=Solid normal
sample_type_code <- substr(colnames(rna_tpm), 14, 15)
is_tumor  <- sample_type_code %in% c("01", "02", "03", "05", "06", "07", "08", "09")
is_normal <- sample_type_code %in% c("10", "11", "12", "13", "14")

cat("Tumor samples:", sum(is_tumor), "\n")
cat("Normal samples:", sum(is_normal), "\n")

tumor_samples  <- colnames(rna_tpm)[is_tumor]
normal_samples <- colnames(rna_tpm)[is_normal]

## ---------- 3. Load sample phenotype / cancer type ----------
pheno <- readRDS("xena_data/tcga_sample_types.rds")

# Standardise column names
names(pheno) <- tolower(names(pheno))
# Key columns expected: sample, _primary_disease or cancer_type
if ("sample" %in% names(pheno)) {
  pheno$barcode <- pheno$sample
} else {
  pheno$barcode <- rownames(pheno)
}

# Cancer type column (try multiple names)
cancer_col <- intersect(c("_primary_disease", "primary_disease", "cancer_type",
                           "_cohort", "cohort", "project_id"), names(pheno))[1]
pheno$cancer_type <- as.character(pheno[[cancer_col]])
pheno$cancer_type <- toupper(gsub("TCGA-|tcga-", "", pheno$cancer_type))

## ---------- 4. Subset to MAA genes found in data ----------
found_genes <- intersect(all_maa_genes, rownames(rna_tpm))
missing_genes <- setdiff(all_maa_genes, rownames(rna_tpm))
cat("\nMAA genes found in data:", length(found_genes), "/", length(all_maa_genes), "\n")
if (length(missing_genes) > 0) {
  cat("Missing from data:", paste(missing_genes, collapse = ", "), "\n")
}

# Update gene layer table to found genes only
gene_layer_table <- gene_layer_table[gene_layer_table$gene %in% found_genes, ]

## ---------- 5. Extract MAA expression matrices ----------
tub_mat_tumor  <- rna_tpm[found_genes, tumor_samples,  drop = FALSE]
tub_mat_normal <- rna_tpm[found_genes, normal_samples, drop = FALSE]

# Remove full matrix from memory
rm(rna_tpm); gc()

## ---------- 6. Load ancillary clinical data ----------
clinical_raw <- readRDS("xena_data/clinical.rds")
if (!is.data.frame(clinical_raw)) clinical_raw <- as.data.frame(clinical_raw)

# Barcode is in the 'sample' column (rownames are just sequential integers)
barcode_col <- intersect(c("sample", "bcr_patient_barcode", "barcode",
                            "submitter_id"), names(clinical_raw))[1]
if (is.na(barcode_col)) stop("Cannot find barcode column in clinical data")
clinical_raw$barcode <- toupper(substr(clinical_raw[[barcode_col]], 1, 15))
cat("Clinical barcode column used:", barcode_col,
    "| Example:", clinical_raw$barcode[1], "\n")

## ---------- 7. Build cancer type vector for tumor samples ----------
pheno_match <- pheno[match(substr(tumor_samples, 1, 16), substr(pheno$barcode, 1, 16)), ]
cancer_vec  <- pheno_match$cancer_type
names(cancer_vec) <- tumor_samples

# Fallback: parse from barcode (TCGA-[TYPE]-...)
missing_ct <- is.na(cancer_vec)
if (any(missing_ct)) {
  cancer_vec[missing_ct] <- toupper(sapply(
    strsplit(tumor_samples[missing_ct], "-"), `[`, 2))
}
cancer_vec[is.na(cancer_vec)] <- "UNKNOWN"

## ---------- 8. CIN score from GISTIC2 ----------
gistic <- readRDS("xena_data/gistic_thresh.rds")
# Remove any non-numeric columns (e.g. "Sample" identifier column)
if (is.data.frame(gistic)) {
  num_cols <- sapply(gistic, is.numeric)
  if (!all(num_cols)) {
    cat("Removing non-numeric GISTIC columns:", paste(names(gistic)[!num_cols], collapse=", "), "\n")
    # If a Sample/rowname column exists, set it as rownames first
    id_col <- names(gistic)[!num_cols][1]
    rownames(gistic) <- gistic[[id_col]]
    gistic <- gistic[, num_cols, drop = FALSE]
  }
  gistic <- as.matrix(gistic)
}
# CIN = mean absolute copy number deviation across all genes (per sample = column)
cin_scores <- colMeans(abs(gistic), na.rm = TRUE)
names(cin_scores) <- colnames(gistic)
cin_scores <- cin_scores[!is.na(cin_scores)]

## ---------- 9. Stemness scores ----------
# Xena stemness format: wide table, 2 rows × N samples
#   Column "sample" = score type label (mRNAsi, mDNAsi)
#   Remaining columns = TCGA barcodes with score values
stem_raw <- readRDS("xena_data/stemness.rds")
stem_df  <- as.data.frame(stem_raw)

if ("sample" %in% names(stem_df) && nrow(stem_df) <= 10 && ncol(stem_df) > 100) {
  # Wide format: pick the mRNAsi row (or first row as fallback)
  use_row <- which(stem_df$sample == "mRNAsi")
  if (length(use_row) == 0) use_row <- 1L
  barcode_cols <- setdiff(names(stem_df), "sample")
  stem_scores  <- setNames(
    as.numeric(unlist(stem_df[use_row[1], barcode_cols])),
    barcode_cols
  )
  cat("Stemness (wide format): mRNAsi scores for",
      length(stem_scores), "samples\n")
} else if (is.numeric(stem_raw) && !is.null(names(stem_raw))) {
  # Named numeric vector
  stem_scores <- stem_raw
  cat("Stemness (vector format):", length(stem_scores), "samples\n")
} else {
  # Tall format: samples as rows, find score column
  id_col    <- intersect(c("sample","barcode","SampleID"), names(stem_df))[1]
  score_col <- intersect(c("mRNAsi","stemness","score","value"), names(stem_df))[1]
  if (is.na(score_col)) {
    num_c <- names(stem_df)[sapply(stem_df, is.numeric)]
    score_col <- num_c[1]
  }
  barcodes    <- if (!is.na(id_col)) stem_df[[id_col]] else rownames(stem_df)
  stem_scores <- setNames(as.numeric(stem_df[[score_col]]), barcodes)
  cat("Stemness (tall format):", length(stem_scores), "samples\n")
}

## ---------- 10. Molecular subtypes ----------
subtype_df <- readRDS("xena_data/subtype.rds")
if (!is.data.frame(subtype_df)) subtype_df <- as.data.frame(subtype_df)
sub_col <- intersect(c("Subtype_Selected", "subtype", "Subtype", "SUBTYPE"),
                     names(subtype_df))[1]
if (!is.na(sub_col)) subtype_df$subtype <- as.character(subtype_df[[sub_col]])

## ---------- 11. Assemble master table (1 row per tumor sample) ----------
cat("Assembling master table...\n")
master <- data.frame(
  barcode     = tumor_samples,
  cancer_type = cancer_vec[tumor_samples],
  stringsAsFactors = FALSE
)
rownames(master) <- tumor_samples

# Add MAA expression (transposed: samples x genes)
tub_expr_df <- as.data.frame(t(tub_mat_tumor))
colnames(tub_expr_df) <- found_genes
master <- cbind(master, tub_expr_df[tumor_samples, , drop = FALSE])

# Add CIN
master$cin_score <- cin_scores[match(substr(master$barcode, 1, 15),
                                      substr(names(cin_scores), 1, 15))]

# Add stemness
master$stemness <- stem_scores[match(substr(master$barcode, 1, 15),
                                      substr(names(stem_scores), 1, 15))]

# Add clinical variables
clin_15 <- substr(clinical_raw$barcode, 1, 15)
m15      <- substr(master$barcode, 1, 15)
idx      <- match(m15, clin_15)

# Copy all available clinical variables (take everything except the barcode column)
clin_to_copy <- setdiff(names(clinical_raw), c(barcode_col, "barcode"))
for (v in clin_to_copy) {
  if (!v %in% names(master)) {   # don't overwrite gene expression columns
    master[[v]] <- clinical_raw[[v]][idx]
  }
}
# Rename common alternate column names to standard names
rename_map <- c(
  "OS"        = "vital_status",
  "OS.time"   = "days_to_death",
  "gender"    = "gender"
)
for (std in names(rename_map)) {
  alt <- rename_map[[std]]
  if (!std %in% names(master) && alt %in% names(master)) {
    master[[std]] <- master[[alt]]
  }
}
# Report what survival data is available
cat("OS non-NA:", sum(!is.na(master$OS)),
    "| OS.time non-NA:", sum(!is.na(master$OS.time)),
    "| age non-NA:", sum(!is.na(master$age_at_initial_pathologic_diagnosis)), "\n")

# Standardise stage column — prioritise AJCC pathologic stage
stage_priority <- c("ajcc_pathologic_tumor_stage", "pathologic_stage",
                     "clinical_stage", "tumor_stage")
stage_candidates <- c(
  intersect(stage_priority, names(master)),
  grep("stage", names(master), ignore.case=TRUE, value=TRUE)
)
stage_candidates <- unique(stage_candidates)
stage_candidates <- stage_candidates[!stage_candidates %in% c("stage","stage_num")]
stage_col <- if (length(stage_candidates) > 0) stage_candidates[1] else NA_character_
if (!is.na(stage_col)) {
  master$stage <- as.character(master[[stage_col]])
  master$stage_num <- as.integer(
    ifelse(grepl("IV|4",      master$stage, ignore.case = TRUE), 4,
    ifelse(grepl("III|3",     master$stage, ignore.case = TRUE), 3,
    ifelse(grepl("II[^I]|2",  master$stage, ignore.case = TRUE), 2,
    ifelse(grepl("I[^VX]|1",  master$stage, ignore.case = TRUE), 1, NA)))))
  cat("Stage column used:", stage_col, "| non-NA:", sum(!is.na(master$stage_num)), "\n")
} else {
  master$stage     <- NA_character_
  master$stage_num <- NA_integer_
  cat("Note: no stage column found in clinical data\n")
}

# Add subtype
if (!is.na(sub_col)) {
  sub15 <- substr(rownames(subtype_df), 1, 15)
  master$subtype <- subtype_df$subtype[match(m15, sub15)]
}

cat("Master table: ", nrow(master), "samples x", ncol(master), "variables\n")
cat("Cancer types:", length(unique(master$cancer_type)), "\n")

## ---------- 12. Save outputs ----------
saveRDS(master,          "xena_data/master.rds")
saveRDS(tub_mat_tumor,   "xena_data/tub_mat_tumor.rds")
saveRDS(tub_mat_normal,  "xena_data/tub_mat_normal.rds")
saveRDS(gene_layer_table,"xena_data/gene_layer_table.rds")
write.csv(master, "results/master_table.csv", row.names = TRUE, quote = FALSE)

cat("\nOutputs saved:\n")
cat("  xena_data/master.rds\n")
cat("  xena_data/tub_mat_tumor.rds\n")
cat("  xena_data/tub_mat_normal.rds\n")
cat("  results/master_table.csv\n")
cat("\nRun RQ scripts next.\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
