# WGCNA Analysis Pipeline for *Medicago sativa*

**Current Version:** v5.1 (Targeted/Custom) | v9.3-Visual-Pro (Unbiased/Automated)  
**Author:** Jiaqing Li, Huanwei Lei  
**Date:** December 10<sup>th</sup>, 2025

## Overview

This repository hosts a robust, publication-ready R pipeline for **Weighted Gene Co-expression Network Analysis (WGCNA)**. Tailored for *Medicago sativa* (Alfalfa) and supports two distinct analytical strategies:

1.  **Targeted Analysis (v5.1):** An end-to-end automated pipeline designed for **manually pre-selected gene sets** (e.g., ~5,000 specific genes). Features modern visualization, automatic phenotype association (Stem Color/Photoperiod), and external enrichment integration.
2.  **Unbiased Discovery (v9.3-Visual-Pro):** A complete, single-script pipeline for **raw transcriptome data** with robust filtering, dynamic module merging, and **publication-grade visualization** optimized for 16GB RAM systems with performance enhancements.

### Core Philosophy

* **System Optimized:** v5.1 for 16GB RAM; v9.3-Visual-Pro specifically tuned for 16GB systems with **memory-aware blockwise computation** and automatic module merging.
* **Scientific Integrity:** Enforces strict data handling. All visualizations derived from real experimental data, with **complete CSV data export** for every plot.
* **Species Specific:** Addresses *M. sativa* annotation gaps by facilitating rigorous external validation (KEGG/GO) with automated gene list preparation.

## 🚀 Key Features

### New in v9.3-Visual-Pro (Publication-Grade Analysis & Performance)
* **Performance-Optimized Computation:** **Blockwise processing (maxBlockSize=9000)** prevents memory overflow on large datasets while maintaining full network integration.
* **Intelligent Module Management:** **Dynamic module merging** when >15 modules detected (threshold: 0.30) combined with **automatic block-wise dendrogram plotting**.
* **Enhanced Visualization Suite:**
  * **Publication-Ready Module-Trait Heatmaps:** **Color-coded module annotation bars** with **intelligent significance marking** (star color adapts to background brightness).
  * **Modern Aesthetic:** **Professional blue-white-red gradient** (colorRamp2) and **consistent styling** across all plots.
  * **Hub Gene Analysis:** **Enhanced heatmaps** with **Set2 color palette** for experimental groups and **network visualizations** with white node borders.
* **Complete Data Export:** **Every PDF visualization** is paired with its corresponding **CSV source data**, including:
  * Full module assignments with kME values
  * Correlation matrices and p-values
  * Hub gene rankings and expression matrices
  * Gene interaction edge lists
* **Robust Validation:** **Automatic output verification** ensures all critical files are generated successfully before completion.

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

### Mode B: Unbiased Discovery (v9.3-Visual-Pro)
*Best for: Raw RNA-seq data requiring noise filtering and novel pattern discovery with publication-grade output.*

1. **Input Data:** Requires `Expression_with_annotation` (raw expression with FPKM/TPM) and `Transcriptome_grouping_summary`.
2. **Run Script:** Execute `WGCNA_v9.3_Visual_Pro.R` – performs everything in one run with real-time progress display.
3. **Key Outputs:**
   * `01_Data_Clean/`: QC plots, outlier detection, filtered data, and **complete gene mapping tables**.
   * `03_Modules/`: Module assignments with kME values **and block-wise dendrograms**.
   * `04_Module_Trait_Correlations/`: **Publication-grade association heatmaps** with module color annotations.
   * `05_Hub_Genes/`: Module-specific hub gene analyses with **enhanced visual styling**.
   * `06_Enrichment_Materials/`: Gene lists for external KEGG/GO analysis.

## 📂 Output Directory Structure

### v9.3-Visual-Pro Structure
```text
Project_Root/
├── 01_Data_Clean/                  # QC plots, filtered data, outlier detection
│   ├── Gene_Name_Mapping_*.csv     # Complete gene ID mappings
│   ├── Sample_Distance_Matrix.csv  # Full distance data for QC plots
│   └── Sample_QC_Stats.csv         # Z-score and outlier statistics
├── 02_Network_Topology/            # Soft-thresholding power analysis
├── 03_Modules/                     # Module dendrograms & statistics
│   ├── Module_Assignment_Full.csv  # Complete module assignments + kME
│   ├── Module_Eigengenes.csv       # Module eigengene matrix
│   ├── Module_Size_Statistics.csv  # Module size distribution
│   └── Gene_Dendrogram_Final_All_Blocks.pdf  # Multi-block dendrograms
├── 04_Module_Trait_Correlations/   # Module-trait association heatmaps
│   ├── Module_Trait_Correlation.csv     # Raw correlation matrix
│   ├── Module_Trait_Pvalue.csv          # Significance matrix
│   ├── Module_Trait_Combined_Results.csv # Integrated results table
│   └── Module_Trait_Heatmap.pdf         # Publication-grade heatmap
├── 05_Hub_Genes/                   # Hub gene networks & expression heatmaps
│   ├── [module_color]/
│   │   ├── Hub_Gene_Ranking_[color].csv           # Top hub genes
│   │   ├── Top_Hub_Genes_Expression_[color].csv   # Expression matrix
│   │   ├── Gene_Interaction_Edges_[color].csv     # Network edge list
│   │   ├── Heatmap_[color].pdf                    # Enhanced heatmap
│   │   └── Network_[color].pdf                    # Network visualization
│   ├── magenta/
│   └── darkgreen/
├── 06_Enrichment_Materials/        # Gene lists for external KEGG/GO tools
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

### v9.3-Visual-Pro - Publication-Grade Visualization & Performance
*Release Date: December 2025*

* **Performance Breakthrough:** **Blockwise computation (maxBlockSize=9000)** solves memory overflow issues on large datasets while maintaining network integrity.
* **Visualization Enhancement:** **Color-coded module annotation bars** in heatmaps with **intelligent significance markers** that adapt to background colors.
* **Publication-Ready Outputs:** **Professional blue-white-red gradient** color scheme and **consistent styling** across all visualizations.
* **Complete Data Transparency:** **Every PDF plot** is paired with its **full CSV source data** for complete reproducibility.
* **Robust Module Management:** **Dynamic merging** of similar modules when count >15, plus **block-wise dendrogram plotting** for multi-block analyses.
* **Enhanced Hub Gene Analysis:** **Modernized heatmaps** with experimental group coloring and **network visualizations** with improved aesthetics.

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
2. **Complete Data Export:** Every visualization is accompanied by its raw data in CSV format, enabling independent verification and re-analysis.
3. **Enrichment Analysis:** Due to the lack of standardized `OrgDb` for *Medicago sativa*, this pipeline **does not** perform internal "guesswork" enrichment. It exports precise gene lists and provides guides for established external tools (AgriGO v2, EggNOG-mapper, KEGG Mapper) to ensure valid biological interpretation.
