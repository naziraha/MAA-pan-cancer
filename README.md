# Microtubule Adaptive Architecture (MAA) in Pan-Cancer

**The Microtubule Adaptive Architecture in Cancer: A Pan-Cancer Systems Analysis Reveals Coordinated Remodelling of Tubulin Isoforms, PTM Enzymes and MAP Regulators as a Multi-Trait Evolutionary Fitness Strategy**

---

## Overview

This repository contains all R analysis code supporting the above manuscript. The study characterises a 136-gene Microtubule Adaptive Architecture (MAA) — comprising tubulin isoforms, post-translational modification (PTM) enzymes, and microtubule-associated protein (MAP)/kinesin regulators — across 33 cancer types in the TCGA Pan-Cancer (PANCAN) cohort (n = 10,535 patients).

**Data source:** UCSC Xena TCGA Pan-Cancer dataset (Toil uniform re-alignment pipeline; Goldman et al., 2020; Vivian et al., 2017)

---

## Repository Structure

```
MAA-pan-cancer/
├── README.md
├── R startup/
│   ├── 00_gene_panel.R          # MAA gene panel definition (136 genes, 4 layers)
│   ├── 01_download_data.R       # Download TCGA TPM data from UCSC Xena
│   ├── 02_build_master.R        # Build master expression + clinical table
│   └── MAA_master_pipeline.R    # Master runner: sources all scripts in order
├── RQ_analysis/
│   ├── RQ0_normal_vs_tumor.R    # Differential expression: tumour vs matched normal
│   ├── RQ1_expression_atlas.R   # Pan-cancer expression heatmap + co-expression modules
│   ├── RQ2_clinical_genomic.R   # Survival, stage trends, CIN, mutation analysis
│   ├── RQ3_fitness_drug_resistance.R  # Multi-trait fitness + drug resistance correlations
│   ├── RQ4_network_pathways.R   # Co-expression network + pathway enrichment
│   └── RQ5_adaptive_states_evolution.R  # Adaptive state clustering + MEI score
└── Validation/
    ├── GEO_external_validation.R  # External GEO cohort validation (5 cohorts)
    └── GEO_validation.R           # Supporting GEO validation utilities
```

---

## Analysis Pipeline

The full analysis is divided into six research questions (RQ0–RQ5):

| Script | Research Question | Key Output |
|--------|------------------|------------|
| `RQ0_normal_vs_tumor.R` | MAA initiation: which genes are dysregulated at cancer onset? | Volcano plot, limma DE results |
| `RQ1_expression_atlas.R` | Pan-cancer expression coordination and co-expression architecture | Heatmap, co-expression modules (k=5) |
| `RQ2_clinical_genomic.R` | Clinical prognosis, stage progression, CIN correlation | Survival forest plots, Cox regression |
| `RQ3_fitness_drug_resistance.R` | Coupling to 7 oncogenic fitness traits and drug resistance | Trait correlation matrix |
| `RQ4_network_pathways.R` | Network topology, hub genes, KEGG/Reactome/Hallmark enrichment | Co-expression network, pathway dotplots |
| `RQ5_adaptive_states_evolution.R` | Adaptive state classification and MAA Evolutionary Index (MEI) | 5-state k-means, MEI ranking |

---

## Requirements

### R Version
R ≥ 4.2.0 (tested on R 4.3.x)

### Core Packages

**Data handling:**
```r
install.packages(c("tidyverse", "data.table", "readxl", "openxlsx"))
```

**Bioconductor packages:**
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("limma", "org.Hs.eg.db", "clusterProfiler",
                       "enrichplot", "ReactomePA", "msigdbr"))
```

**Survival analysis:**
```r
install.packages(c("survival", "survminer"))
```

**Visualisation:**
```r
install.packages(c("ggplot2", "pheatmap", "ComplexHeatmap", "ggrepel",
                   "RColorBrewer", "viridis", "patchwork", "cowplot"))
```

**Network analysis:**
```r
install.packages(c("igraph", "ggraph", "tidygraph"))
```

---

## How to Run

### Option 1: Run the full pipeline (recommended)

```r
# In RStudio, open R project.Rproj, then:
source("R startup/MAA_master_pipeline.R")
```

This sources all scripts in the correct order (00 → 01 → 02 → RQ0 → RQ1 → RQ2 → RQ3 → RQ4 → RQ5 → Validation).

### Option 2: Run individual scripts

```r
# Step 1: Define gene panel
source("R startup/00_gene_panel.R")

# Step 2: Download data (~1.5 GB, ~30 min depending on connection)
source("R startup/01_download_data.R")

# Step 3: Build master table
source("R startup/02_build_master.R")

# Step 4–9: Run each RQ script
source("RQ_analysis/RQ0_normal_vs_tumor.R")
source("RQ_analysis/RQ1_expression_atlas.R")
source("RQ_analysis/RQ2_clinical_genomic.R")
source("RQ_analysis/RQ3_fitness_drug_resistance.R")
source("RQ_analysis/RQ4_network_pathways.R")
source("RQ_analysis/RQ5_adaptive_states_evolution.R")

# Step 10: External validation
source("Validation/GEO_external_validation.R")
```

### Data Download Note

The TCGA Pan-Cancer TPM dataset (~1.5 GB compressed) is downloaded automatically by `01_download_data.R` from the UCSC Xena Toil hub:

```
https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz
```

Due to file size, raw data files are **not** included in this repository. All processed `.rds` intermediate files are saved locally in `xena_data/` after the first run.

---

## Output Files

All results are saved to:
- `results/` — Excel files (`.xlsx`) for all analysis outputs
- `figures/` — Publication-quality PDF/PNG figures

Key output files per RQ:

| RQ | Results file | Key figures |
|----|-------------|-------------|
| RQ0 | `RQ0_normal_vs_tumor.xlsx` | `RQ0_volcano_tumor_vs_normal.pdf` |
| RQ1 | `RQ1_expression_atlas.xlsx` | `Fig1_pancancer_heatmap.pdf`, `Fig2_coexpression_matrix.pdf` |
| RQ2 | `survival_results_all_genes.xlsx`, `multivariate_cox_all_genes.xlsx` | `Fig3_survival_forest_OS.pdf` |
| RQ3 | `RQ3_fitness_drug_resistance.xlsx` | `Fig5_fitness_trait_correlations.pdf` |
| RQ4 | `RQ4_KEGG_enrichment.xlsx`, `RQ4_hub_genes.xlsx` | `RQ4_coexpression_network.pdf` |
| RQ5 | `RQ5_adaptive_states.xlsx`, `RQ5_priority_candidates.xlsx` | `RQ5_adaptive_state_heatmap.pdf`, `RQ5_PCA_fitness_landscape.pdf` |

---

## MAA Gene Panel

The 136-gene MAA panel is defined in `R startup/00_gene_panel.R` and organised into four functional layers:

| Layer | n genes | Examples |
|-------|---------|---------|
| Layer 1 — Tubulin isoforms | 21 | TUBA1A, TUBA1B, TUBB3, TUBB4A, TUBD1, TUBE1 |
| Layer 2 — MAPs / Kinesins / Motors | 83 | GTSE1, KIFC1, KIF18B, KIF23, BIRC5, EML4, CAMSAP1-3 |
| Layer 3 — PTM writer enzymes | 19 | CDK1, ATAT1, TTLL1–TTLL13, VASH1, VASH2, SETD2 |
| Layer 4 — PTM eraser enzymes | 13 | HDAC6, SIRT2, TTL, AGTPBP1, AGBL1–AGBL5 |

Of the 136 panel genes, 132 were matched in the TCGA RNA-Seq matrix. Four genes (MAP11, DCDC5, SOGA1, KIAA1493) were absent from the expression dataset and excluded from downstream analyses.

---

## External Validation Cohorts (GEO)

| GEO Accession | Cancer Type | n |
|--------------|------------|---|
| GSE14520 | Liver hepatocellular carcinoma (LIHC) | 247 |
| GSE96058 | Breast cancer (BRCA) | 3,273 |
| GSE72094 | Lung adenocarcinoma (LUAD) | 442 |
| GSE62254 | Stomach adenocarcinoma (STAD) | 300 |
| GSE39582 | Colon adenocarcinoma (COAD) | 585 |

---

## Citation

If you use this code, please cite:

> [Author list]. The Microtubule Adaptive Architecture in Cancer: A Pan-Cancer Systems Analysis Reveals Coordinated Remodelling of Tubulin Isoforms, PTM Enzymes and MAP Regulators as a Multi-Trait Evolutionary Fitness Strategy. *[Journal]*, [Year]. DOI: [to be updated upon publication]

---

## Key References

- Goldman MJ, et al. (2020). Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology*, 38, 675–678.
- Vivian J, et al. (2017). Toil enables reproducible, open source, big biomedical data analyses. *Nature Biotechnology*, 35, 314–316.
- Liu J, et al. (2018). An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell*, 173, 400–416.

---

## License

This code is released under the MIT License. See `LICENSE` for details.

---

## Contact

For questions regarding the analysis pipeline, please open a GitHub Issue or contact the corresponding author.

