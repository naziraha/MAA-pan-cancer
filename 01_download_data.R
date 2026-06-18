# ============================================================
# MAA PAN-CANCER ANALYSIS — DATA DOWNLOAD SCRIPT v2.1
# Downloads: tcga_RSEM_gene_tpm (TCGA PANCAN Toil hub, UCSC Xena)
# NOTE: Toil TPM file uses Ensembl IDs (ENSG...) — script maps to HUGO
# ============================================================

## ---------- 0. Setup ----------
if (!dir.exists("xena_data")) dir.create("xena_data")
if (!dir.exists("results"))   dir.create("results")
if (!dir.exists("figures"))   dir.create("figures")

required_pkgs <- c("data.table", "R.utils", "AnnotationDbi", "org.Hs.eg.db")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("AnnotationDbi", "org.Hs.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else {
      install.packages(p)
    }
  }
}))
library(data.table); library(R.utils)
library(AnnotationDbi); library(org.Hs.eg.db)

## ---------- 1. Load gene panel ----------
source("R startup/00_gene_panel.R")   # defines all_maa_genes

## ---------- 2. Build Ensembl ID → HUGO symbol map ----------
cat("Building Ensembl → HUGO symbol map for MAA genes...\n")
sym2ens <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_maa_genes,
  columns = c("SYMBOL", "ENSEMBL"),
  keytype = "SYMBOL"
)
sym2ens <- sym2ens[!is.na(sym2ens$ENSEMBL), ]
# Keep only the first Ensembl ID per symbol (some have multiple)
sym2ens <- sym2ens[!duplicated(sym2ens$SYMBOL), ]

# Reverse map: ENSEMBL → SYMBOL (for fast lookup during file scan)
ens2sym <- setNames(sym2ens$SYMBOL, sym2ens$ENSEMBL)
target_ensembl <- sym2ens$ENSEMBL
cat("MAA genes with Ensembl IDs:", length(target_ensembl), "/", length(all_maa_genes), "\n")

## ---------- 3. Download TPM file (if needed) ----------
tpm_gz  <- "xena_data/tcga_RSEM_gene_tpm.gz"
tpm_rds <- "xena_data/rna_tpm_pancan.rds"
xena_url <- "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"

if (!file.exists(tpm_rds)) {
  if (!file.exists(tpm_gz)) {
    cat("Downloading tcga_RSEM_gene_tpm.gz (~1.5 GB)...\n")
    options(timeout = 3600)
    download.file(xena_url, destfile = tpm_gz, method = "libcurl", mode = "wb")
    cat("Download complete.\n")
  }

  ## ---------- 4. Decompress (Windows-safe: R.utils::gunzip) ----------
  tpm_tsv <- sub("\\.gz$", "", tpm_gz)
  if (!file.exists(tpm_tsv)) {
    cat("Decompressing (~8 GB uncompressed, 5-10 min)...\n")
    R.utils::gunzip(tpm_gz, destname = tpm_tsv, remove = FALSE, overwrite = TRUE)
    cat("Decompression complete.\n")
  } else {
    cat("Decompressed TSV already exists, skipping gunzip.\n")
  }

  ## ---------- 5. Scan file in chunks, match on Ensembl IDs ----------
  # The Toil file uses versioned Ensembl IDs (e.g. ENSG00000123456.3)
  # Strip version suffix before matching
  cat("Scanning file for MAA gene rows (Ensembl ID matching)...\n")
  cat("(Reading 2000 lines at a time to stay within R memory limits)\n")

  con <- file(tpm_tsv, open = "r", encoding = "UTF-8")
  header_line <- readLines(con, n = 1)   # sample barcode header
  kept_lines  <- character(0)
  found_symbols <- character(0)
  chunk_size  <- 2000L
  total_read  <- 0L

  repeat {
    chunk <- readLines(con, n = chunk_size)
    if (length(chunk) == 0L) break
    total_read <- total_read + length(chunk)

    # Extract Ensembl ID (first field), strip version number (.XX)
    raw_ids      <- sub("\t.*", "", chunk)
    stripped_ids <- sub("\\.[0-9]+$", "", raw_ids)   # ENSG00000123456.3 → ENSG00000123456

    keep_idx <- stripped_ids %in% target_ensembl
    if (any(keep_idx)) {
      kept_lines    <- c(kept_lines,    chunk[keep_idx])
      found_symbols <- c(found_symbols, ens2sym[stripped_ids[keep_idx]])
    }
    if (total_read %% 10000L == 0L)
      cat("  Lines scanned:", total_read, "| MAA genes found:", length(kept_lines), "\n")
  }
  close(con)
  cat("Scan complete. Lines scanned:", total_read,
      "| MAA gene rows found:", length(kept_lines), "\n")

  if (length(kept_lines) == 0)
    stop("No MAA genes found in file. Check Ensembl ID mapping.")

  ## ---------- 6. Parse extracted rows ----------
  cat("Parsing MAA gene expression matrix...\n")
  rna_tpm <- data.table::fread(
    text       = paste(c(header_line, kept_lines), collapse = "\n"),
    sep        = "\t",
    header     = TRUE,
    data.table = FALSE
  )

  # Replace Ensembl IDs with HUGO gene symbols as row names
  row_ids_stripped   <- sub("\\.[0-9]+$", "", rna_tpm[[1]])
  rna_tpm[[1]]       <- ens2sym[row_ids_stripped]
  # Drop rows where mapping failed (NA)
  rna_tpm <- rna_tpm[!is.na(rna_tpm[[1]]), ]
  # Remove duplicate gene symbols (keep first)
  rna_tpm <- rna_tpm[!duplicated(rna_tpm[[1]]), ]
  rownames(rna_tpm)  <- rna_tpm[[1]]
  rna_tpm            <- rna_tpm[, -1, drop = FALSE]
  rna_tpm            <- as.matrix(rna_tpm)

  cat("Final matrix:", nrow(rna_tpm), "MAA genes x", ncol(rna_tpm), "samples\n")
  cat("Gene symbols in matrix:\n"); print(sort(rownames(rna_tpm)))

  cat("Saving RDS...\n")
  saveRDS(rna_tpm, tpm_rds, compress = "gzip")
  cat("Saved:", tpm_rds, "\n")

} else {
  cat("TPM RDS already exists:", tpm_rds, "\n")
}

## ---------- 7. Verify ancillary files ----------
ancillary <- c("xena_data/clinical.rds","xena_data/gistic_thresh.rds",
               "xena_data/mc3.rds","xena_data/stemness.rds","xena_data/subtype.rds")
cat("\nVerifying ancillary data files:\n")
for (f in ancillary) {
  if (file.exists(f)) {
    cat(" ✓", f, "(", round(file.size(f)/1e6,1), "MB )\n")
  } else {
    cat(" ✗ MISSING:", f, "\n")
  }
}

## ---------- 8. Download phenotype table (Windows-safe) ----------
phenotype_gz  <- "xena_data/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"
phenotype_tsv <- "xena_data/TCGA_phenotype_denseDataOnlyDownload.tsv"
phenotype_url <- "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz"

if (!file.exists("xena_data/tcga_sample_types.rds")) {
  cat("\nDownloading TCGA phenotype table...\n")
  options(timeout = 600)
  if (!file.exists(phenotype_gz))
    download.file(phenotype_url, destfile = phenotype_gz, method = "libcurl", mode = "wb")
  # Decompress with R.utils (no system gzip needed)
  if (!file.exists(phenotype_tsv))
    R.utils::gunzip(phenotype_gz, destname = phenotype_tsv,
                    remove = FALSE, overwrite = TRUE)
  pheno <- data.table::fread(phenotype_tsv, sep = "\t",
                               header = TRUE, data.table = FALSE)
  saveRDS(pheno, "xena_data/tcga_sample_types.rds", compress = "gzip")
  cat("Saved: xena_data/tcga_sample_types.rds\n")
} else {
  cat("Sample type annotation already present.\n")
}

cat("\nData download complete. Run 02_build_master.R next.\n")


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
