[![DOI](https://zenodo.org/badge/1109997346.svg)](https://doi.org/10.5281/zenodo.17831030)
# WGCNA Analysis Pipeline 

**Current Versions:**

  * **Step 1: Preliminary Analysis:** v1.1.0
  * **Step 2: Advanced Analysis:** v2.4

## 1\. System Requirements

  * **OS:** Windows 10/11 (Tested). Compatible with Linux/macOS with minor path adjustments.
  * **R Version:** R 4.5.2 or higher is highly recommended.
  * **Hardware:**
      * **CPU:** Multi-core processor recommended (Scripts are optimized for 4-6 threads).
      * **RAM:** **16GB+ is strongly recommended.**
          * *Note:* The v2.4 Advanced script includes active memory management (garbage collection) to ensure stability on 16GB systems, but WGCNA network construction remains memory-intensive.

## 2\. R Dependencies

The pipeline is divided into two stages. Ensure all packages below are installed to run the full workflow.

### Core Analysis & Statistics (Step 1)

  * **`WGCNA`**: Main package for weighted correlation network analysis.
  * **`impute`** (Bioconductor): Required by WGCNA for handling missing data.
  * **`preprocessCore`** (Bioconductor): Required for WGCNA normalization functions.
  * **`fastcluster`**: Accelerates hierarchical clustering operations.
  * **`dynamicTreeCut`**: Adaptive branch pruning for module detection.

### Data Manipulation (Step 1 & 2)

  * **`dplyr`** / **`tidyverse`**: Advanced data frame manipulation.
  * **`stringr`**: String cleaning and RegEx operations.
  * **`reshape2`**: Data transformation (melting/casting).
  * **`data.table`**: High-performance file reading (Crucial for v2.4 speed).

### Advanced Visualization (Step 2 - New)

  * **`ggplot2`**: Base system for publication-quality static plotting.
  * **`pheatmap`**: Professional heatmaps with dynamic clustering.
  * **`igraph`**: Network topology analysis and visualization.
  * **`patchwork`**: Arranging multiple `ggplot` figures into composite figures.
  * **`ggrepel`**: Smart label placement to prevent text overlapping.
  * **`RColorBrewer`** / **`viridis`**: Scientific color palettes.

### Functional Enrichment (Step 2 - New)

  * **`clusterProfiler`**: Comprehensive enrichment analysis (KEGG/GO).
  * **`enrichplot`**: Visualization of enrichment results (Dot plots, etc.).
  * **`org.Hs.eg.db`**: Genome-wide annotation for Human (Replace with `org.Mm.eg.db` for Mouse, etc.).
  * **`GOSemSim`**: Semantic similarity measurement (Dependency for clusterProfiler).

## 3\. Installation Guide

To ensure both Script A (Preliminary) and Script B (Advanced) run smoothly, copy and paste the following code into your R console. This includes the "Smart Repair" logic to handle difficult Bioconductor packages.

```r
# ==============================================================================
# Unified Dependency Installer
# ==============================================================================

# 1. Install CRAN packages
cran_pkgs <- c("WGCNA", "dplyr", "ggplot2", "reshape2", "stringr", 
               "RColorBrewer", "viridis", "BiocManager", "fastcluster", 
               "dynamicTreeCut", "pheatmap", "igraph", "patchwork", 
               "ggrepel", "data.table")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# 2. Install Bioconductor packages
# Note: This step checks for dependencies that often cause installation failures
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install core annotation and enrichment packages
# CHANGE "org.Hs.eg.db" if you are not analyzing Human data
bioc_pkgs <- c("GO.db", "impute", "preprocessCore", "clusterProfiler", 
               "enrichplot", "GOSemSim", "AnnotationDbi", "org.Hs.eg.db")

BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

cat("Installation complete. Please restart your R session before running the pipeline.\n")
```

## 4\. Pipeline Workflow

### Part 1: Preliminary Analysis (Script v1.1.0)

**Focus:** Data Cleaning & Network Construction.

  * **Inputs:** Raw FPKM/Counts CSV file.
  * **Key Operations:**
      * Data cleaning (filtering low expression/low variance genes).
      * Outlier detection.
      * Soft-threshold selection ($\beta$).
      * Module detection via `blockwiseModules`.
  * **Outputs:** `03_Network/WGCNA_Network_Object.RData` (Critical input for Part 2).

### Part 2: Advanced Analysis (Script v2.4)

**Focus:** Mining, Visualization, and Biological Interpretation.

  * **Inputs:** The `.RData` file and Preprocessed Matrix from Part 1.
  * **Key Features:**
      * **Robustness:** Auto-detects data inconsistencies and fixes vector/dataframe errors.
      * **Hub Gene Mining:** Calculates Module Membership (kME) to identify top 20 driver genes per module.
      * **Visualization:** Generates publication-ready PDFs (300 DPI) for:
          * Hierarchical Clustering Dendrograms.
          * Network Topology Validation.
          * Hub Gene Interaction Networks (`igraph`).
          * Module-Module Correlation Heatmaps (`pheatmap` with dynamic breaks).
      * **Enrichment:** Automated KEGG and GO enrichment for top modules.
  * **Memory Optimization:** Uses aggressive garbage collection (`gc()`) to run safely on 16GB RAM machines.
