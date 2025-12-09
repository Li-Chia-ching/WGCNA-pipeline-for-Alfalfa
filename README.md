[![DOI](https://zenodo.org/badge/1109997346.svg)](https://doi.org/10.5281/zenodo.17831030)
# WGCNA Analysis Pipeline for *Medicago sativa*

**Current Version:** v5.1 (Targeted/Custom) | v7.0-Final (Unbiased/Automated)  
**Author:** Jiaqing Li, Huanwei Lei  
**Date:** December 8<sup>th</sup>, 2025  

## Overview

This repository hosts a robust, publication-ready R pipeline for **Weighted Gene Co-expression Network Analysis (WGCNA)**. Tailored for *Medicago sativa* (Alfalfa) and supports two distinct analytical strategies:

1.  **Targeted Analysis (v5.1):** An end-to-end automated pipeline designed for **manually pre-selected gene sets** (e.g., ~5,000 specific genes). Features modern visualization, automatic phenotype association (Stem Color/Photoperiod), and external enrichment integration.
2.  **Unbiased Discovery (v7.0-Final):** A complete, single-script pipeline for **raw transcriptome data** with robust filtering, dynamic module merging, and comprehensive visualization optimized for 15GB RAM systems.

### Core Philosophy

* **System Optimized:** v5.1 for 16GB RAM; v7.0-Final specifically tuned for 15GB AMD Ryzen systems with memory-aware handling.
* **Scientific Integrity:** Enforces strict data handling. All visualizations derived from real experimental data.
* **Species Specific:** Addresses *M. sativa* annotation gaps by facilitating rigorous external validation (KEGG/GO).

## 🚀 Key Features

### New in v7.0-Final (Complete Automation & Memory Optimization)
* **Full Pipeline Integration:** Merges previously separate steps into one seamless workflow.
* **Robust Data Handling:** Smart detection/cleaning of `:fpkm` or `:tpm` suffixes; automatic sample alignment.
* **Advanced QC:** Z-score based outlier detection with dynamic labeling to avoid overlap.
* **Dynamic Module Merging:** Auto-merges when >15 modules detected (threshold: 0.30) for manageable results.
* **Memory-Optimized Filtering:** Keeps top 25,000 variable genes (TPM>1 threshold) for 15GB RAM systems.
* **Comprehensive Visualization:** Modern color schemes (viridis), combined QC dashboards, and publication-ready PDFs.

### Features in v5.1 (Modern Visualization for Targeted Analysis)
* **Automatic Phenotype Mapping:** Parses sample names (e.g., `LLGS`, `SLRS`) to create binary trait matrices for **Stem Color** (Red/Green) and **Photoperiod** (Long/Short).
* **Next-Gen Visualization:** Uses **ComplexHeatmap** for dynamic row/column annotation and significance starring.
* **Hub Gene Visualization:** Generates high-resolution expression heatmaps (`pheatmap`) and interaction networks (`igraph`) for key modules.
* **Third-Party Integration:** Automatically formats and exports gene lists for **KEGG Mapper**, **KOBAS**, and **ShinyGO**, with template scripts for result visualization in R.

## 🛠 Dependencies & Installation

Requires **R 4.5.2+**.

### 1. Standard Dependencies
```r
# Set Mirrors (Tsinghua)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioc/")

# Core Packages
install.packages(c("WGCNA", "dplyr", "stringr", "ggplot2", "RColorBrewer", 
                   "viridis", "pheatmap", "igraph", "circlize", "patchwork", 
                   "ggrepel", "gridExtra", "corrplot", "scales", "data.table", 
                   "reshape2", "devtools"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GO.db", "impute", "preprocessCore"))
```

### 2. ComplexHeatmap (Critical for v5.1 Heatmaps)
```r
library(devtools)
install_github("jokergoo/ComplexHeatmap")
```

## 📖 Workflow Usage Guide

### Mode A: Targeted Analysis (v5.1)
*Best for: Analyzing specific gene subsets (manually filtered lists) with known phenotype groupings.*

1. **Input Data:** Load as `Summary_of_manually_filtered_genes` (Expression) and `Transcriptome_grouping_summary` (Phenotype).
2. **Run Script:** Execute `WGCNA_v5.1_Modern.R`.
3. **Key Outputs:**
   * `04_Module_Trait_Correlations/`: Master heatmap showing module-trait relationships.
   * `05_Hub_Genes/`: Detailed analysis of specific modules (e.g., Magenta, Darkgreen).
   * `09_Third_Party_Analysis/`: Gene lists ready for KEGG Mapper.

### Mode B: Unbiased Discovery (v7.0-Final)
*Best for: Raw RNA-seq data requiring noise filtering and novel pattern discovery.*

1. **Input Data:** Requires `Expression_with_annotation` (raw expression with FPKM/TPM) and `Transcriptome_grouping_summary`.
2. **Run Script:** Execute `WGCNA_v7.0_Final.R` – performs everything in one run.
3. **Key Outputs:**
   * `01_Data_Clean/`: QC plots, outlier detection, filtered data.
   * `03_Modules/`: Module assignments with kME values.
   * `04_Module_Trait_Correlations/`: Association heatmaps.
   * `05_Hub_Genes/`: Module-specific hub gene analyses.
   * `06_Enrichment_Materials/`: Gene lists for external KEGG/GO analysis.

## 📂 Output Directory Structure

### v7.0-Final Structure
```text
Project_Root/
├── 01_Data_Clean/                  # QC plots, filtered data, outlier detection
├── 02_Network_Topology/            # Soft-thresholding power analysis
├── 03_Modules/                     # Module dendrograms & statistics
├── 04_Module_Trait_Correlations/   # Module-trait association heatmaps
├── 05_Hub_Genes/                   # Hub gene networks & expression heatmaps
│   ├── magenta/
│   └── darkgreen/
├── 06_Enrichment_Materials/        # Gene lists for external KEGG/GO tools
├── 08_Summary_Reports/             # Executive Summary & Session Info
└── Input_Data.RData                # Checkpoint for processed data
```

### v5.1 Structure
```text
Project_Root/
├── 01_Data_Clean/                  # Preprocessed expression & trait matrices
├── 02_Network_Topology/            # Soft-thresholding power analysis plots
├── 03_Modules/                     # Module dendrograms & statistics
├── 04_Module_Trait_Correlations/   # COMPLEX HEATMAPS (Key Result)
├── 05_Hub_Genes/                   # Hub gene networks & expression heatmaps
├── 06_Enrichment_Analysis/         # Guides for external tools
├── 07_Cytoscape_Files/             # Edges/Nodes files for Cytoscape
├── 08_Summary_Reports/             # Executive Summary & Session Info
└── 09_Third_Party_Analysis/        # Ready-to-use inputs for KEGG/GO tools
```

## 📝 Version History

### v7.0-Final - Complete Pipeline Integration
*Release Date: December 2025*
* **Full Workflow:** Consolidated previously separate steps into single automated script.
* **Memory Optimization:** Specifically tuned for 15GB RAM (AMD Ryzen 7 PRO).
* **Robust Data Handling:** Automatic suffix cleaning (`:fpkm`/`:tpm`), sample alignment, and outlier detection.
* **Dynamic Module Management:** Auto-merges modules when count >15 for manageable results.
* **Enhanced Visualization:** Modern color schemes, combined QC dashboards, and comprehensive PDF outputs.
* **KEGG Integration:** Exports formatted gene lists with species guidance (*M. truncatula* for alfalfa).

### v5.1 - The Modernization Update
*Release Date: December 2025*
* **Targeted Workflow:** Tailored for **manually filtered gene sets** rather than raw expression filtering.
* **Visualization:** Introduced `ComplexHeatmap` for publication-quality module-trait correlations.
* **Automation:** Automatic parsing of sample names into biological factors (Stem Color, Photoperiod).
* **External Support:** Dedicated modules for **KEGG Mapper** and **KOBAS** integration.

## 🔬 Scientific Integrity Statement

**Strict Data Policy:** This pipeline adheres to the highest standards of scientific data integrity.

1. **No Synthetic Data:** All visualizations (heatmaps, dendrograms, networks) are generated directly from provided expression data.
2. **Enrichment Analysis:** Due to the lack of standardized `OrgDb` for *Medicago sativa*, this pipeline **does not** perform internal "guesswork" enrichment. It exports precise gene lists and provides guides for established external tools (AgriGO v2, EggNOG-mapper, KEGG Mapper) to ensure valid biological interpretation.
