# WGCNA Analysis Pipeline for *Medicago sativa*

**Current Version:** v5.0 (Targeted/Custom) | v2.7.1 (Unbiased/Legacy)  
**Author:** Jiaqing Li  
**Date:** December 8<sup>th</sup>, 2025

## Overview

This repository hosts a robust, publication-ready R pipeline for **Weighted Gene Co-expression Network Analysis (WGCNA)**. It is tailored for *Medicago sativa* (Alfalfa) and supports two distinct analytical strategies:

1.  **Targeted Analysis (v5.0):** An end-to-end automated pipeline designed for **manually pre-selected gene sets** (e.g., \~5,000 specific genes). It features modern visualization, automatic phenotype association (Stem Color/Photoperiod), and external enrichment integration.
2.  **Unbiased Discovery (v2.7.1):** A classic workflow for **raw transcriptome data**, featuring strict low-expression and low-variance filtering to discover novel modules from scratch. (Need  Basic Script WGCNA_pipeline.R v1.1.0 to perform Preliminary Analysis first.)

### Core Philosophy

  * **System Optimized:** Designed for 16GB RAM environments using memory-aware data handling and multi-threading.
  * **Scientific Integrity:** Enforces strict data handling standards. All visualizations are derived from real experimental data; no synthetic data is used.
  * **Species Specific:** Addresses the lack of specific annotation packages for *M. sativa* by facilitating rigorous external validation (KEGG/GO).

-----

## 🚀 Key Features

### New in v5.0 (Modern Visualization & Automation)

  * **Automatic Phenotype Mapping:** Automatically parses sample names (e.g., `LLGS`, `SLRS`) to create binary trait matrices for **Stem Color** (Red/Green) and **Photoperiod** (Long/Short).
  * **Next-Gen Visualization:**
      * **Module-Trait Heatmaps:** Utilizes **ComplexHeatmap** for dynamic row/column annotation and significance starring.
      * **Hub Gene Visualization:** Generates high-resolution expression heatmaps (`pheatmap`) and interaction networks (`igraph`) for key modules.
  * **Third-Party Integration:** Automatically formats and exports gene lists for **KEGG Mapper**, **KOBAS**, and **ShinyGO**, and includes a template script to visualize these external results back in R.

### Legacy Features (v1.1.0 - v2.7.1)

  * **Strict Quality Control:** Filters genes with \>50% zero expression and removes the bottom 25% low-variance genes.
  * **Robust Construction:** Uses `blockwiseModules` with multi-threading (6 threads).
  * **Dynamic Thresholding:** Uses quantile-based thresholding (Top 10%) for reliable Cytoscape export.

-----

## 🛠 Dependencies & Installation

This pipeline requires **R 4.5.2+**.

### 1\. Standard Dependencies

Run the following in R to install CRAN and Bioconductor packages:

```r
# Set Mirrors (Tsinghua)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioc/")

# Core Packages
install.packages(c("WGCNA", "dplyr", "stringr", "ggplot2", "RColorBrewer", 
                   "viridis", "pheatmap", "igraph", "patchwork", "ggrepel", 
                   "gridExtra", "corrplot", "scales", "data.table", "reshape2", "devtools"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GO.db", "impute", "preprocessCore"))
```

### 2\. ComplexHeatmap (Critical for v5.0)

The v5.0 pipeline relies on **ComplexHeatmap** for advanced plotting. Please install it directly from GitHub to ensure compatibility:

```r
library(devtools)
install_github("jokergoo/ComplexHeatmap")
```

-----

## 📖 Workflow Usage Guide

### Mode A: Targeted Analysis (v5.0)

*Best for: Analyzing a specific subset of genes (e.g., manually filtered lists) with known phenotype groupings.*

1.  **Input Data:** Load your data into RStudio as `Summary_of_manually_filtered_genes` (Expression) and `Transcriptome_grouping_summary` (Phenotype).
2.  **Run Script:** Execute `WGCNA_v5.0_Modern.R`.
3.  **Key Outputs:**
      * `04_Module_Trait_Correlations/`: The master heatmap showing relationships between modules and Red/Green or Long/Short traits.
      * `05_Hub_Genes/`: Detailed analysis of specific modules (e.g., Magenta, Darkgreen).
      * `09_Third_Party_Analysis/`: Gene lists ready for KEGG Mapper.

### Mode B: Unbiased Discovery (v1.1 + v2.7)

*Best for: Raw RNA-seq data where you need to filter noise and find novel patterns.*

1.  **Run Basic Script (v1.1.0):** Performs data cleaning (variance/expression filtering) and network construction.
2.  **Run Advanced Script (v2.7.1):** Generates visualizations and hub gene lists based on the calculated network.

-----

## 📂 Output Directory Structure (v5.0)

The v5.0 pipeline creates a comprehensive directory structure:

```text
Project_Root/
├── 01_Data_Clean/                  # Preprocessed expression & trait matrices
├── 02_Network_Topology/            # Soft-thresholding power analysis plots
├── 03_Modules/                     # Module dendrograms & statistics
├── 04_Module_Trait_Correlations/   # COMPLEX HEATMAPS (Key Result)
│   ├── Module_Trait_Correlation_Heatmap_Modern.pdf
│   └── Detailed_Correlation_Results_Full.csv
├── 05_Hub_Genes/                   # Hub gene networks & expression heatmaps
│   ├── magenta/
│   └── darkgreen/
├── 06_Enrichment_Analysis/         # Guides for external tools
├── 07_Cytoscape_Files/             # Edges/Nodes files for Cytoscape
├── 08_Summary_Reports/             # Executive Summary & Session Info
└── 09_Third_Party_Analysis/        # Ready-to-use inputs for KEGG/GO tools
    ├── KEGG_Input/
    └── Visualization_Template.R
```

-----

## 📝 Version History

### v5.0 - The Modernization Update (Current)

*Release Date: December 2025*

  * **New Workflow:** tailored for **manually filtered gene sets** rather than raw expression filtering.
  * **Visualization:** Introduced `ComplexHeatmap` for publication-quality module-trait correlations with significance annotations (`***`).
  * **Automation:** Added automatic parsing of sample names (e.g., `BPGS_1`) into biological factors (Stem Color, Photoperiod).
  * **External Support:** Added dedicated modules to prepare data for **KEGG Mapper** and **KOBAS**.

### v2.7.1 - The Robustness Patch

  * **Fix:** Addressed edge cases where identical variance across Module Eigengenes caused console warnings.
  * **Optimization:** Enhanced memory monitoring for 16GB RAM systems.

### v2.7 - The Integrity Update

  * **Scientific Integrity:** Removed synthetic data generation. Added workflows for *real* external enrichment.
  * **Network:** Implemented dynamic quantile thresholding (Top 10%) for Cytoscape export.

### v1.1.0 - Basic Pipeline

  * Initial release including 50% zero-expression filtering and 25% low-variance filtering.

-----

## 🔬 Scientific Integrity Statement

**Strict Data Policy:**
This pipeline is designed to adhere to the highest standards of scientific data integrity.

1.  **No Synthetic Data:** All visualizations (heatmaps, dendrograms, networks) are generated directly from your provided expression data.
2.  **Enrichment Analysis:** Due to the lack of a standardized `OrgDb` for *Medicago sativa*, this pipeline **does not** perform internal "guesswork" enrichment. Instead, it exports precise gene lists and provides a guide for using established external tools (AgriGO v2, EggNOG-mapper, KEGG Mapper) to ensure valid biological interpretation.
