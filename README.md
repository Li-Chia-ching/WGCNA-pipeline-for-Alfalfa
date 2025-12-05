# WGCNA Analysis Pipeline (v1.1.0)

## 1. System Requirements
* **OS:** Windows 10 (Tested). Compatible with Linux/macOS with path adjustments.
* **R Version:** R 4.5.2 or higher.
* **Hardware:** * CPU: Multi-core processor recommended (Script is optimized for 6 threads).
    * RAM: 16GB+ recommended (WGCNA is memory-intensive).

## 2. R Dependencies
The script relies on a specific set of CRAN and Bioconductor packages. 

### Core Analysis & Statistics
* **`WGCNA`**: Main package for weighted correlation network analysis.
* **`impute`** (Bioconductor): Required by WGCNA for handling missing data.
* **`GO.db`** (Bioconductor): Annotation database for biological context.
* **`preprocessCore`** (Bioconductor): Often required for WGCNA normalization functions.

### Data Manipulation & Visualization
* **`dplyr`**: Data frame manipulation.
* **`stringr`**: String cleaning and RegEx operations.
* **`reshape2`**: Data transformation.
* **`ggplot2`**: Static plotting.
* **`RColorBrewer`**: Color palettes for modules.

### Performance & Utilities
* **`BiocManager`**: For managing Bioconductor installations.
* **`fastcluster`**: (Implicit) Highly recommended for accelerating WGCNA clustering.
* **`dynamicTreeCut`**: (Implicit) For adaptive branch pruning.

## 3. Installation Guide
Copy and paste the following code into your R console to install all dependencies:

```r
# 1. Install CRAN packages
install.packages(c("WGCNA", "dplyr", "ggplot2", "reshape2", 
                   "stringr", "RColorBrewer", "BiocManager", 
                   "fastcluster", "dynamicTreeCut"))

# 2. Install Bioconductor packages
BiocManager::install(c("GO.db", "impute", "preprocessCore"))
