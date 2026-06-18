# ============================================================
# MAA PAN-CANCER ANALYSIS — GENE PANEL DEFINITION v2.0
# 131 unique canonical HGNC genes | 4 functional layers
# ============================================================
# Run this first; all downstream scripts source it.

## LAYER 1: Tubulin Structural Isoforms (21 genes)
# Alpha-tubulins
tubulin_alpha <- c("TUBA1A", "TUBA1B", "TUBA1C",
                   "TUBA3C", "TUBA3D", "TUBA3E",
                   "TUBA4A", "TUBA8", "TUBAL3")
# Beta-tubulins
tubulin_beta  <- c("TUBB",   "TUBB1",  "TUBB2A", "TUBB2B",
                   "TUBB3",  "TUBB4A", "TUBB4B", "TUBB6",
                   "TUBB8",  "TUBB8B")
# Non-alpha/beta
tubulin_other <- c("TUBD1", "TUBE1")

tubulin_isoforms <- unique(c(tubulin_alpha, tubulin_beta, tubulin_other))
# n = 21

## LAYER 2: MAPs, Kinesins, Motors & Scaffolding Proteins (85 genes)
# Classical structural MAPs
maps_structural <- c("MAP1A",  "MAP1B",  "MAP1S",  "MAP2",   "MAP4",
                     "MAP6",   "MAP7",   "MAP9",   "MAP10",  "MAP11",
                     "MAPT",   "MAP7D1", "MAP7D3")
# Doublecortin / MT nucleation MAPs
maps_dcx        <- c("DCX", "DCLK1", "DCDC1", "DCDC5")
# Autophagy-linked MAPs (LC3 family)
maps_autophagy  <- c("MAP1LC3A", "MAP1LC3B", "GABARAPL1")
# CAMSAP / minus-end stabilizers
maps_camsap     <- c("CAMSAP1", "CAMSAP2", "CAMSAP3")
# Plus-end tracking / anchoring
maps_plusend    <- c("CLIP1",  "GTSE1",  "SPAG5",  "EML4",  "EML6",
                     "SNPH",   "MTCL1",  "MTUS1",  "PLEKHA5")
# Kinesin motors (mitotic & interphase)
kinesins_mitotic     <- c("KIF11",  "KIF14",  "KIF15",  "KIF18A", "KIF18B",
                           "KIF20A", "KIF20B", "KIF23",  "KIFC1",  "KIFC3",
                           "CENPE",  "KIF2C",  "SPAG5")
kinesins_interphase  <- c("KIF1A",  "KIF2A",  "KIF2B",  "KIF3A",  "KIF3C",
                           "KIF5B",  "KIF16B", "KIF17",  "KIF21B", "KIF26A",
                           "KIF26B", "KIF27")
# Centrosomal / scaffolding
maps_centrosome <- c("CEP72", "CEP295", "CDK5RAP3")
# Signalling / cytoskeletal regulators co-localizing with MTs
maps_signalling <- c("BIRC5",  "RACGAP1", "EZR",    "MYO10",
                     "ARL8A",  "ARL8B",   "ARL4C",  "FSCN1",
                     "BRCA1",  "MASTL",   "PRKAA1", "RAD51D",
                     "NF1",    "SS18",    "SOGA1",  "LZTS1",
                     "CACYBP", "DIXDC1",  "FES",    "HDGFL3",
                     "FNTA",   "TBCC",    "TCP1",   "KIAA1493")

maps_genes <- unique(c(maps_structural, maps_dcx, maps_autophagy, maps_camsap,
                       maps_plusend,    kinesins_mitotic, kinesins_interphase,
                       maps_centrosome, maps_signalling))
# n = ~85

## LAYER 3: PTM Writers — Tubulin Modification Enzymes (19 genes)
# Polyglutamylases (TTLL family, writers)
ttll_writers   <- c("TTLL1",  "TTLL3",  "TTLL4",  "TTLL5",  "TTLL6",
                    "TTLL7",  "TTLL8",  "TTLL9",  "TTLL10", "TTLL11",
                    "TTLL12", "TTLL13")
# Acetylation writers
ace_writers    <- c("ATAT1", "KAT2A")
# Detyrosinase writers (vasohibin family)
vash_writers   <- c("VASH1", "VASH2", "SVBP")
# Methylation writer
meth_writers   <- c("SETD2")
# Kinase (phosphorylation)
kinase_writers <- c("CDK1")

ptm_writers <- unique(c(ttll_writers, ace_writers, vash_writers,
                         meth_writers, kinase_writers))
# n = 19

## LAYER 4: PTM Erasers — Tubulin Demodification Enzymes (14 genes)
# Deacetylases
deacetylases <- c("HDAC1", "HDAC2", "HDAC5", "HDAC6", "HDAC8", "SIRT2")
# Tubulin tyrosine ligase (tyrosination — reverse of VASH)
ttl_fam      <- c("TTL")
# Carboxypeptidase CCPs (deglutamylation / detyrosination)
ccp_fam      <- c("AGTPBP1",  # CCP1
                  "AGBL2",    # CCP2
                  "AGBL3",    # CCP3
                  "AGBL4",    # CCP4
                  "AGBL5",    # CCP5
                  "AGBL1")    # CCP6

ptm_erasers <- unique(c(deacetylases, ttl_fam, ccp_fam))
# n = 13

## FULL PANEL — 4 layers merged & deduplicated
all_maa_genes <- unique(c(tubulin_isoforms, maps_genes, ptm_writers, ptm_erasers))
cat("Total MAA genes:", length(all_maa_genes), "\n")
# Expected: ~131

## Gene-layer lookup table (for colour-coding, stratified analyses)
gene_layer_table <- rbind(
  data.frame(gene = tubulin_isoforms, layer = "tubulin_isoform",  stringsAsFactors = FALSE),
  data.frame(gene = setdiff(maps_genes,    c(tubulin_isoforms, ptm_writers, ptm_erasers)),
             layer = "maps_kinesin",       stringsAsFactors = FALSE),
  data.frame(gene = setdiff(ptm_writers,   c(tubulin_isoforms)),
             layer = "ptm_writer",         stringsAsFactors = FALSE),
  data.frame(gene = setdiff(ptm_erasers,   c(tubulin_isoforms, maps_genes, ptm_writers)),
             layer = "ptm_eraser",         stringsAsFactors = FALSE)
)
# remove any inadvertent duplicates
gene_layer_table <- gene_layer_table[!duplicated(gene_layer_table$gene), ]
rownames(gene_layer_table) <- gene_layer_table$gene

## Layer colour palette (used in all figures)
layer_colors <- c(
  tubulin_isoform = "#E41A1C",   # red
  maps_kinesin    = "#377EB8",   # blue
  ptm_writer      = "#4DAF4A",   # green
  ptm_eraser      = "#FF7F00"    # orange
)

cat("Gene panel loaded. Summary:\n")
print(table(gene_layer_table$layer))
cat("\nAll genes:\n")
print(sort(all_maa_genes))


# =============================================================================
# Session Information (for reproducibility)
# =============================================================================
cat("\n--- Session Information ---\n")
print(sessionInfo())
