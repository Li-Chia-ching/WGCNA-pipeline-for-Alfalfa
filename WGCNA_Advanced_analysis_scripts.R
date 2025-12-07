# ==============================================================================
# WGCNA ADVANCED ANALYSIS v2.7.1 (16GB RAM优化版)
# 功能：发表级可视化、富集分析、核心基因挖掘
# 物种：Medicago sativa (紫花苜蓿/苜蓿)
# 
# 主要改进：
#   1. 修复了网络动态阈值问题（移除固定阈值0.6）
#   2. 移除了随机生成的功能分类饼图，避免学术诚信问题
#   3. 添加了详细的外部富集分析指导
#   4. 修复了所有括号匹配错误
#   5. 针对16GB内存系统进行了全面优化
#   6. 计算相关性矩阵 module_cor 之后，插入一段温和的提示逻辑让代码更健壮
# ==============================================================================

#### 1. 初始化和包加载 ####
# ==============================================================================
cat("╔═══════════════════════════════════════════════════╗\n")
cat("║         WGCNA ADVANCED ANALYSIS v2.7              ║\n")
cat("║         16GB RAM Optimized Version               ║\n")
cat("║         Species: Medicago sativa                 ║\n")
cat("╚═══════════════════════════════════════════════════╝\n\n")

# 设置关键选项
options(stringsAsFactors = FALSE)
options(warn = 1)  # 立即显示警告
options(max.print = 100)  # 限制控制台输出行数
options(encoding = "UTF-8")  # 设置编码为UTF-8

# 函数：检查并安装缺失的包
check_and_install_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("Installing package:", pkg, "...\n"))
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
      cat(paste("  ✓", pkg, "installed\n"))
    }
  }
}

# 函数：监控内存使用情况（简化版）
monitor_memory <- function(step_name) {
  # 使用gc()来获取内存信息
  gc_info <- gc()
  mem_used <- sum(gc_info[, 2])  # 第二列是MB
  cat(paste("  Memory usage (", step_name, "): ", 
            round(mem_used, 1), " MB\n", sep = ""))
}

# 步骤1.1: 加载基础分析包
cat("[Step 1.1] Loading essential packages...\n")
essential_pkgs <- c("WGCNA", "ggplot2", "dplyr", "RColorBrewer", "viridis")
check_and_install_packages(essential_pkgs)

# 步骤1.2: 加载可视化包
cat("[Step 1.2] Loading visualization packages...\n")
viz_pkgs <- c("pheatmap", "igraph", "patchwork", "ggrepel")
check_and_install_packages(viz_pkgs)

# 步骤1.3: 加载富集分析包
cat("[Step 1.3] Loading enrichment analysis packages...\n")
enrichment_pkgs <- c("clusterProfiler", "enrichplot")
check_and_install_packages(enrichment_pkgs)

# 注意：由于Medicago sativa没有专门的OrgDb包，我们需要使用其他方法
cat("Note: Medicago sativa analysis - using plant-specific databases\n")
cat("For enrichment analysis, gene lists will be saved for external tools\n")

# 配置并行处理（内存感知）
cat("[Step 1.4] Configuring parallel processing...\n")
n_cores <- min(4, parallel::detectCores() - 2)  # 保留2个核心给系统
WGCNA::enableWGCNAThreads(nThreads = n_cores)
cat(paste("  Using", n_cores, "CPU cores for parallel processing\n"))

# 设置统一的绘图主题（发表级）
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.2),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
)

#### 2. 数据加载和验证 ####
# ==============================================================================
cat("\n[Step 2] Loading and validating Step1 analysis results...\n")

# 增强的数据加载函数，包含内存监控
load_and_validate_data <- function() {
  
  # 1. 检查并加载网络对象
  if (!file.exists("03_Network/WGCNA_Network_Object.RData")) {
    stop("❌ ERROR: Cannot find file 03_Network/WGCNA_Network_Object.RData")
  }
  
  cat("  Loading WGCNA network object...\n")
  load("03_Network/WGCNA_Network_Object.RData")
  monitor_memory("After loading network object")
  gc()  # 强制垃圾回收
  
  # 2. 加载表达矩阵（内存高效方式）
  if (!file.exists("01_InputData/Preprocessed_Expression_Matrix.csv")) {
    stop("❌ ERROR: Cannot find file 01_InputData/Preprocessed_Expression_Matrix.csv")
  }
  
  cat("  Loading expression matrix (memory-optimized)...\n")
  
  # 使用data.table进行快速、内存高效的加载（如果可用）
  if (requireNamespace("data.table", quietly = TRUE)) {
    library(data.table)
    datExpr <- fread("01_InputData/Preprocessed_Expression_Matrix.csv", 
                     showProgress = FALSE)
    gene_names <- datExpr[[1]]
    datExpr <- as.data.frame(datExpr[, -1])
    rownames(datExpr) <- gene_names
    datExpr <- as.data.frame(t(datExpr))
  } else {
    # 备选方案：使用base R，但限制行数以节省内存
    cat("  Note: Consider installing 'data.table' for faster file loading\n")
    datExpr <- read.csv("01_InputData/Preprocessed_Expression_Matrix.csv", 
                        row.names = 1, 
                        nrows = 50000)  # 限制读取行数
    datExpr <- as.data.frame(t(datExpr))
  }
  
  # 3. 修复module_colors的命名问题（如果存在）
  if (is.null(names(module_colors))) {
    if (length(module_colors) == ncol(datExpr)) {
      names(module_colors) <- colnames(datExpr)
      cat("  ✓ Fixed gene names in module_colors\n")
    } else {
      warning("⚠️ module_colors length doesn't match gene count. Attempting alignment...")
    }
  }
  
  # 4. 对齐数据
  common_genes <- intersect(colnames(datExpr), names(module_colors))
  if (length(common_genes) < 100) {
    warning("⚠️ Few common genes found between expression data and module assignments")
  }
  
  datExpr <- datExpr[, common_genes, drop = FALSE]
  module_colors <- module_colors[common_genes]
  
  # 输出统计信息
  cat(paste("  ✓ Data loaded successfully:\n"))
  cat(paste("    - Samples:", nrow(datExpr), "\n"))
  cat(paste("    - Genes:", ncol(datExpr), "\n"))
  
  # 计算模块统计
  module_table <- table(module_colors)
  non_grey_modules <- sum(names(module_table) != "grey")
  cat(paste("    - Total modules detected:", length(module_table), "\n"))
  cat(paste("    - Non-grey modules:", non_grey_modules, "\n"))
  
  # 显示前5大模块
  module_sorted <- sort(module_table, decreasing = TRUE)
  if (length(module_sorted) > 5) {
    cat("    - Top 5 largest modules:\n")
    for (i in 1:min(5, length(module_sorted))) {
      cat(paste("      ", names(module_sorted)[i], ": ", module_sorted[i], " genes\n", sep = ""))
    }
  }
  
  monitor_memory("Data loading completed")
  
  return(list(
    datExpr = datExpr, 
    module_colors = module_colors, 
    net = net, 
    MEs = MEs, 
    softPower = softPower
  ))
}

# 执行数据加载
cat("Starting data loading...\n")
data_list <- load_and_validate_data()

# 将数据加载到全局环境
list2env(data_list, envir = .GlobalEnv)

# 清理临时变量以释放内存
rm(data_list)
gc()

# 创建输出目录结构
cat("  Creating output directories...\n")
output_dirs <- c(
  "06_Advanced_Results",
  "06_Advanced_Results/Figures",
  "06_Advanced_Results/Tables", 
  "06_Advanced_Results/Enrichment",
  "06_Advanced_Results/Enrichment/Figures",  # 用于存放真实分析后的图片
  "06_Advanced_Results/Network_Data"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat(paste("    Created directory:", dir, "\n"))
  }
}

monitor_memory("Directory creation completed")

#### 3. 可视化：基因层次聚类树 (发表级) ####
# ==============================================================================
cat("\n[Step 3] Generating gene hierarchical clustering dendrogram...\n")

generate_dendrogram <- function() {
  cat("  Plotting gene dendrogram...\n")
  
  # 创建PDF版本 - 使用英文标题避免编码问题
  pdf("06_Advanced_Results/Figures/1_Gene_Dendrogram_and_Modules.pdf", 
      width = 14, height = 8)
  
  WGCNA::plotDendroAndColors(
    dendro = net$dendrograms[[1]], 
    colors = module_colors[net$blockGenes[[1]]],
    groupLabels = "Module Assignment",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    guideAll = FALSE,
    main = "Gene Hierarchical Clustering and WGCNA Module Detection",
    marAll = c(1, 5, 3, 1),
    cex.main = 1.2,
    cex.colorLabels = 0.7,
    cex.dendroLabels = 0.5
  )
  
  dev.off()
  
  # 创建PNG版本用于快速查看
  png("06_Advanced_Results/Figures/1_Gene_Dendrogram_and_Modules.png", 
      width = 1400, height = 800, res = 150)
  
  WGCNA::plotDendroAndColors(
    dendro = net$dendrograms[[1]], 
    colors = module_colors[net$blockGenes[[1]]],
    groupLabels = "Module Assignment",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    main = "Gene Hierarchical Clustering and WGCNA Module Detection"
  )
  
  dev.off()
  
  cat("  → Output files: Figures/1_Gene_Dendrogram_and_Modules.pdf/.png\n")
  cat("    Purpose: Shows how genes are assigned to different co-expression modules through hierarchical clustering\n")
}

# 执行绘图
generate_dendrogram()
gc()

#### 4. 可视化：网络拓扑分析 ####
# ==============================================================================
cat("\n[Step 4] Generating network topology analysis plots...\n")

generate_topology_plots <- function() {
  cat("  Calculating network topology metrics...\n")
  
  # 重新计算拓扑指标（轻量级操作）
  powers <- c(1:10, seq(12, 20, 2))
  sft <- WGCNA::pickSoftThreshold(datExpr, 
                                  powerVector = powers, 
                                  networkType = "unsigned", 
                                  verbose = 0)
  
  # 准备数据框用于绘图
  sft_df <- data.frame(
    Power = sft$fitIndices$Power,
    R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
    MeanK = sft$fitIndices$mean.k.
  )
  
  cat("  Creating publication-quality topology plots...\n")
  
  # 图1：无尺度拓扑拟合
  p1 <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = R2)) + 
    ggplot2::geom_point(size = 2.5, color = "#E41A1C", alpha = 0.8) + 
    ggplot2::geom_line(color = "#E41A1C", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0.85, linetype = "dashed", 
                        color = "gray40", linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = softPower, linetype = "dashed", 
                        color = "#377EB8", linewidth = 0.6) +
    ggplot2::annotate("text", x = softPower, y = max(sft_df$R2), 
                      label = paste("Selected:", softPower), 
                      vjust = -0.5, hjust = 1.2, 
                      color = "#377EB8", fontface = "bold", size = 3) +
    ggplot2::labs(
      title = "Scale-free Topology Model Fit",
      x = "Soft Threshold (Power)",
      y = expression("Scale Free Topology Model Fit (R"^2*")")
    ) +
    ggplot2::scale_x_continuous(breaks = powers) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # 图2：平均连接度
  p2 <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = MeanK)) + 
    ggplot2::geom_point(size = 2.5, color = "#377EB8", alpha = 0.8) + 
    ggplot2::geom_line(color = "#377EB8", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = softPower, linetype = "dashed", 
                        color = "#377EB8", linewidth = 0.6) +
    ggplot2::labs(
      title = "Mean Gene Connectivity",
      x = "Soft Threshold (Power)", 
      y = "Mean Connectivity (log10 scale)"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = powers) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.4),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # 合并图形
  combined_plot <- p1 + p2 + 
    patchwork::plot_annotation(
      title = "WGCNA Network Topology Analysis",
      subtitle = paste("Selected soft threshold power =", softPower),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
      )
    )
  
  # 保存高质量版本
  ggplot2::ggsave("06_Advanced_Results/Figures/2_Network_Topology_Analysis.pdf", 
                  combined_plot, width = 14, height = 6, dpi = 300)
  
  ggplot2::ggsave("06_Advanced_Results/Figures/2_Network_Topology_Analysis.png", 
                  combined_plot, width = 14, height = 6, dpi = 300)
  
  cat("  → Output files: Figures/2_Network_Topology_Analysis.pdf/.png\n")
  cat("    Purpose: Validates if the selected soft threshold power is appropriate, ensuring the network approximates scale-free topology\n")
}

# 执行绘图
generate_topology_plots()
gc()

#### 5. 功能富集分析 - Medicago sativa专用版（安全版本） ####
# ==============================================================================
cat("\n[Step 5] Preparing data for external functional enrichment analysis for Medicago sativa...\n")

perform_enrichment_analysis <- function() {
  
  # 识别最大的模块（排除grey模块）
  module_table <- table(module_colors)
  module_summary <- data.frame(
    Module = names(module_table),
    Size = as.numeric(module_table),
    stringsAsFactors = FALSE
  )
  
  # 排除grey模块
  module_summary <- module_summary[module_summary$Module != "grey", ]
  
  # 选择前5个最大的模块进行富集分析
  module_summary <- module_summary[order(-module_summary$Size), ]
  target_modules <- module_summary$Module[1:min(5, nrow(module_summary))]
  
  cat(paste("  Preparing gene lists for top", length(target_modules), "largest modules:\n"))
  cat(paste("    ", paste(target_modules, collapse = ", "), "\n"))
  
  # 为每个模块准备基因列表
  for (mod in target_modules) {
    module_genes <- names(module_colors)[module_colors == mod]
    cat(paste("    Processing module", mod, ":", length(module_genes), "genes\n"))
    
    if (length(module_genes) > 10) {  # 只有足够多基因的模块才保存
      # 保存基因列表
      write.table(module_genes,
                  paste0("06_Advanced_Results/Enrichment/Genes_", mod, ".txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      # 创建模块表达模式热图（真实数据）
      cat("      Creating expression pattern heatmap (real data)...\n")
      tryCatch({
        # 随机选择最多50个基因进行可视化（确保多样性）
        n_genes_to_plot <- min(50, length(module_genes))
        set.seed(123)  # 可重复的随机选择
        selected_genes <- sample(module_genes, n_genes_to_plot)
        
        expr_subset <- datExpr[, selected_genes, drop = FALSE]
        
        # Z-score标准化
        expr_scaled <- t(scale(expr_subset))
        
        # 设置颜色方案
        color_palette <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)
        
        # 创建热图 - 这是真实数据，不是随机生成的
        pheatmap::pheatmap(expr_scaled,
                           main = paste("Expression Pattern - Module", mod),
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           show_colnames = ncol(expr_scaled) <= 20,
                           show_rownames = FALSE,  # 基因太多不显示名称
                           color = color_palette,
                           fontsize_row = 6,
                           fontsize_col = 7,
                           border_color = NA,
                           filename = paste0("06_Advanced_Results/Enrichment/Figures/Expression_Pattern_", mod, ".pdf"),
                           width = 10,
                           height = 8)
        
        cat(paste("      Expression pattern heatmap created for module", mod, "\n"))
        
      }, error = function(e) {
        cat(paste("      Could not create expression pattern heatmap:", e$message, "\n"))
      })
    } else {
      cat(paste("      Skipping - too few genes for meaningful analysis\n"))
    }
    
    # 强制垃圾回收
    gc()
  }
  
  # 创建详细的富集分析说明文件
  enrichment_note <- c(
    "=================================================================================",
    "EXTERNAL FUNCTIONAL ENRICHMENT ANALYSIS GUIDE FOR MEDICAGO SATIVA",
    paste("Generated:", Sys.time()),
    "=================================================================================",
    "",
    "IMPORTANT ACADEMIC INTEGRITY NOTICE:",
    "  This script does NOT perform actual functional enrichment analysis due to database",
    "  limitations for Medicago sativa. It ONLY prepares gene lists for external tools.",
    "  Do NOT use any randomly generated functional category plots in publications.",
    "",
    "SPECIES INFORMATION:",
    "  - Organism: Medicago sativa (Alfalfa)",
    "  - Family: Fabaceae (Legumes)",
    "  - Common names: Alfalfa, Lucerne",
    "",
    "GENE LISTS PREPARED:",
    "  - Genes_*.txt files contain gene symbols for each co-expression module",
    "  - These are REAL gene names from your expression data (not simulated)",
    "",
    "RECOMMENDED EXTERNAL TOOLS FOR ENRICHMENT ANALYSIS:",
    "",
    "1. AgriGO v2.0 (http://systemsbiology.cau.edu.cn/agriGOv2/)",
    "   - Specifically designed for agricultural species",
    "   - Supports GO enrichment analysis",
    "   - Upload your Genes_*.txt files directly",
    "",
    "2. EggNOG-mapper (http://eggnog-mapper.embl.de/)",
    "   - Functional annotation and orthology assignments",
    "   - Use Medicago truncatula as reference species",
    "   - Download results and perform enrichment in R",
    "",
    "3. KEGG Mapper (https://www.kegg.jp/kegg/mapper/)",
    "   - Use Medicago truncatula (KEGG code: mtr) as reference",
    "   - Search for your genes in pathway maps",
    "",
    "4. PLAZA (https://bioinformatics.psb.ugent.be/plaza/)",
    "   - Plant comparative genomics platform",
    "   - Functional annotations and orthologous groups",
    "",
    "5. LegumeIP (http://plantgrn.noble.org/LegumeIP/)",
    "   - Legume-specific database",
    "   - Integrated omics data for legumes",
    "",
    "HOW TO CREATE FUNCTIONAL CATEGORY PLOTS AFTER EXTERNAL ANALYSIS:",
    "  1. Run enrichment analysis using one of the tools above",
    "  2. Download the results (usually CSV or Excel format)",
    "  3. In R, use the following code template to create bar plots or pie charts:",
    "",
    "  # Example R code for creating functional enrichment bar plot",
    "  library(ggplot2)",
    "  library(dplyr)",
    "  ",
    "  # Load your enrichment results",
    "  enrichment_results <- read.csv('your_enrichment_results.csv')",
    "  ",
    "  # Create bar plot (recommended over pie charts)",
    "  p <- ggplot(enrichment_results, aes(x = reorder(Term, -log10(PValue)), y = -log10(PValue))) +",
    "    geom_bar(stat = 'identity', fill = 'steelblue') +",
    "    coord_flip() +",
    "    labs(title = 'Functional Enrichment Results',",
    "         x = 'GO Term / Pathway',",
    "         y = '-log10(P-value)') +",
    "    theme_minimal()",
    "  ",
    "  ggsave('Functional_Enrichment_BarPlot.pdf', p, width = 10, height = 8)",
    "",
    "NEXT STEPS:",
    "  1. Upload Genes_*.txt files to external enrichment tools",
    "  2. Use Medicago truncatula annotations as reference where needed",
    "  3. Import results back to R for visualization",
    "  4. Create publication-quality figures with REAL enrichment data",
    "",
    "================================================================================="
  )
  
  writeLines(enrichment_note, "06_Advanced_Results/Enrichment/Enrichment_Analysis_Guide.txt")
  
  # 创建简化的R脚本模板用于结果可视化
  create_visualization_template <- function() {
    template <- c(
      "# ==============================================================================",
      "# R Script Template for Visualizing External Enrichment Results",
      "# Use this after running enrichment analysis on external tools",
      "# ==============================================================================",
      "",
      "library(ggplot2)",
      "library(dplyr)",
      "library(RColorBrewer)",
      "",
      "# 1. Load your enrichment results (modify file path as needed)",
      "enrichment_data <- read.csv('path/to/your/enrichment_results.csv')",
      "",
      "# 2. Basic data cleaning and filtering",
      "# Filter significant results (adjust p-value cutoff as needed)",
      "significant_results <- enrichment_data %>%",
      "  filter(PValue < 0.05) %>%",
      "  arrange(PValue) %>%",
      "  head(15)  # Show top 15 most significant terms",
      "",
      "# 3. Create bar plot (recommended for enrichment visualization)",
      "p_bar <- ggplot(significant_results, aes(x = reorder(Term, -log10(PValue)), y = -log10(PValue))) +",
      "  geom_bar(stat = 'identity', fill = 'steelblue', alpha = 0.8) +",
      "  coord_flip() +",
      "  labs(title = 'Top Enriched Functional Terms',",
      "       x = 'Functional Term',",
      "       y = '-log10(P-value)') +",
      "  theme_minimal(base_size = 11) +",
      "  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),",
      "        axis.text.y = element_text(size = 9))",
      "",
      "# Save bar plot",
      "ggsave('06_Advanced_Results/Enrichment/Figures/Real_Enrichment_BarPlot.pdf',",
      "       p_bar, width = 10, height = 7, dpi = 300)",
      "",
      "# 4. Alternative: Dot plot with gene counts and p-values",
      "if ('GeneCount' %in% colnames(significant_results) && 'TotalGenes' %in% colnames(significant_results)) {",
      "  p_dot <- ggplot(significant_results, aes(x = GeneCount/TotalGenes, y = reorder(Term, GeneCount/TotalGenes))) +",
      "    geom_point(aes(size = GeneCount, color = -log10(PValue))) +",
      "    scale_color_gradient(low = 'blue', high = 'red') +",
      "    labs(title = 'Enrichment Dot Plot',",
      "         x = 'Gene Ratio',",
      "         y = 'Functional Term',",
      "         size = 'Gene Count',",
      "         color = '-log10(P-value)') +",
      "    theme_minimal()",
      "  ",
      "  ggsave('06_Advanced_Results/Enrichment/Figures/Real_Enrichment_DotPlot.pdf',",
      "         p_dot, width = 10, height = 7, dpi = 300)",
      "}",
      "",
      "# 5. Create a summary table",
      "write.csv(significant_results,",
      "          '06_Advanced_Results/Enrichment/Significant_Enrichment_Results.csv',",
      "          row.names = FALSE)",
      "",
      "cat('Enrichment visualization completed. Check the Enrichment/Figures directory.\\n')",
      ""
    )
    
    writeLines(template, "06_Advanced_Results/Enrichment/Visualization_Template.R")
  }
  
  create_visualization_template()
  
  cat("  → Output files: Gene lists saved in '06_Advanced_Results/Enrichment/'\n")
  cat("  → Expression heatmaps (real data) saved in '06_Advanced_Results/Enrichment/Figures/'\n")
  cat("  → Detailed guide for external analysis: 'Enrichment_Analysis_Guide.txt'\n")
  cat("  → R template for result visualization: 'Visualization_Template.R'\n")
  cat("    Purpose: Prepares gene lists for external functional enrichment analysis\n")
  cat("    IMPORTANT: No random functional categories generated to ensure academic integrity\n")
}

# 执行富集分析准备
perform_enrichment_analysis()
gc()

#### 6. 可视化：核心基因分析 ####
# ==============================================================================
cat("\n[Step 6] Analyzing and visualizing hub genes...\n")

# 修复后的核心基因分析函数（使用动态阈值）
analyze_hub_genes <- function(target_module, datExpr, MEs, n_top = 20) {
  cat(paste("    Processing module", target_module, "...\n"))
  
  # 提取属于目标模块的基因
  module_genes <- names(module_colors)[module_colors == target_module]
  
  if (length(module_genes) == 0) {
    cat(paste("      Warning: No genes in the module\n"))
    return(NULL)
  }
  
  if (length(module_genes) < n_top) {
    cat(paste("      Warning: Module has only", length(module_genes), "genes\n"))
    n_top <- length(module_genes)
  }
  
  # 计算模块成员度 (kME) - 修复后的版本
  cat("      Calculating module membership (kME)...\n")
  tryCatch({
    # 确保MEs中有目标模块的特征基因
    me_name <- paste0("ME", target_module)
    if (!me_name %in% colnames(MEs)) {
      cat(paste("      Error: Cannot find eigengene", me_name, "in MEs\n"))
      return(NULL)
    }
    
    # 提取表达数据和特征基因
    expr_data <- datExpr[, module_genes, drop = FALSE]
    me_data <- MEs[, me_name, drop = FALSE]
    
    # 计算kME
    kME_values <- WGCNA::signedKME(expr_data, me_data)
    
    # 确保kME_values是矩阵或数据框
    if (is.null(kME_values) || nrow(kME_values) == 0) {
      cat("      Error: Failed to calculate kME, returned NULL\n")
      return(NULL)
    }
    
    # 转换为数据框 - 修复列选择错误
    kME <- data.frame(
      Gene = module_genes,
      kME = kME_values[, 1],  # 明确选择第一列
      stringsAsFactors = FALSE
    )
    
    # 按kME绝对值排序
    kME_sorted <- kME[order(-abs(kME$kME)), ]
    hub_genes <- kME_sorted[1:min(n_top, nrow(kME_sorted)), ]
    top_gene_names <- hub_genes$Gene
    
    cat(paste("      Identified", length(top_gene_names), "hub genes\n"))
    
    # 1. 创建核心基因表达热图
    if (length(top_gene_names) >= 3) {
      cat("      Creating expression heatmap...\n")
      
      # --- 检查方差齐性，并给出友好提示 ---
      if (length(unique(apply(MEs_clean, 2, var))) == 1) {
        cat("  ℹ️ 提示: 所有模块特征基因方差相同，这通常源于数据标准化步骤，不影响生物学解释。\n")
      }
      
      tryCatch({
        expr_subset <- datExpr[, top_gene_names, drop = FALSE]
        
        # Z-score标准化
        expr_scaled <- t(scale(expr_subset))
        
        # 设置颜色方案
        color_palette <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)
        
        # 创建热图
        pheatmap::pheatmap(expr_scaled,
                           main = paste("Hub Gene Expression - Module", target_module),
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           show_colnames = ncol(expr_scaled) <= 25,  # 样本太多时不显示列名
                           show_rownames = TRUE,
                           color = color_palette,
                           fontsize_row = 7,
                           fontsize_col = 6,
                           border_color = NA,
                           filename = paste0("06_Advanced_Results/Figures/3_Heatmap_Hub_", target_module, ".pdf"),
                           width = 9,
                           height = 8)
        
        # 保存PNG版本
        grDevices::png(paste0("06_Advanced_Results/Figures/3_Heatmap_Hub_", target_module, ".png"),
                       width = 900, height = 800, res = 150)
        
        pheatmap::pheatmap(expr_scaled,
                           main = paste("Hub Gene Expression - Module", target_module),
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           show_colnames = ncol(expr_scaled) <= 25,
                           color = color_palette)
        
        grDevices::dev.off()
        
      }, error = function(e) {
        cat(paste("      Error creating heatmap:", e$message, "\n"))
      })
    } else {
      cat("      Too few genes, skipping heatmap\n")
    }
    
    # 2. 创建核心基因互作网络（如果有足够基因）- 使用动态阈值
    if (length(top_gene_names) >= 5) {
      cat("      Creating interaction network with dynamic thresholding...\n")
      tryCatch({
        # 计算邻接矩阵（基于相关性）
        adj_matrix <- abs(cor(datExpr[, top_gene_names, drop = FALSE], 
                              use = "pairwise.complete.obs"))
        
        # 使用动态阈值而非固定阈值0.6
        # 保留最强的连接（前10%）
        upper_tri_values <- adj_matrix[upper.tri(adj_matrix)]
        
        if (length(upper_tri_values) > 0) {
          # 计算90%分位数作为阈值
          threshold <- quantile(upper_tri_values, probs = 0.90, na.rm = TRUE)
          
          # 如果计算出的阈值太小，使用最小值确保有连接
          if (threshold < 0.3 || is.na(threshold)) {
            threshold <- 0.3  # 安全阈值
            cat(paste("      Using safe threshold:", threshold, "\n"))
          } else {
            cat(paste("      Using dynamic threshold (90th percentile):", round(threshold, 3), "\n"))
          }
          
          # 应用动态阈值
          adj_matrix[adj_matrix < threshold] <- 0
        } else {
          cat("      Warning: No valid correlation values for threshold calculation\n")
          # 使用保守的默认阈值
          adj_matrix[adj_matrix < 0.5] <- 0
        }
        
        diag(adj_matrix) <- 0
        
        # 创建igraph对象
        g <- igraph::graph_from_adjacency_matrix(adj_matrix, 
                                                 mode = "undirected", 
                                                 weighted = TRUE,
                                                 diag = FALSE)
        
        # 移除孤立节点
        g <- igraph::delete_vertices(g, igraph::degree(g) == 0)
        
        if (igraph::vcount(g) > 1) {  # 只有在有网络时才绘图
          # 计算布局
          layout <- igraph::layout_with_fr(g)
          
          # 创建PDF
          grDevices::pdf(paste0("06_Advanced_Results/Figures/3_Network_Hub_", target_module, ".pdf"),
                         width = 9, height = 9)
          
          # 计算节点属性
          node_degree <- igraph::degree(g)
          node_colors <- grDevices::colorRampPalette(c("lightblue", "darkblue"))(100)
          degree_scaled <- as.integer(cut(node_degree, breaks = 100))
          
          # 创建增强版图形
          igraph::plot.igraph(g,
                              layout = layout,
                              vertex.size = 12 + 1.5 * node_degree,  # 大小与连接度成正比
                              vertex.color = node_colors[degree_scaled],
                              vertex.label.cex = 0.7,
                              vertex.label.color = "black",
                              vertex.label.dist = 0.8,
                              vertex.frame.color = "white",
                              edge.width = igraph::E(g)$weight * 3,
                              edge.color = "gray60",
                              main = paste("Hub Gene Interaction Network - Module", target_module))
          
          # 添加图例
          graphics::legend("topleft",
                           legend = c("High connectivity", "Low connectivity"),
                           fill = c("darkblue", "lightblue"),
                           bty = "n",
                           cex = 0.7)
          
          grDevices::dev.off()
          
          # 保存网络数据
          write.csv(igraph::as_data_frame(g, what = "edges"),
                    paste0("06_Advanced_Results/Network_Data/Network_Edges_", target_module, ".csv"),
                    row.names = FALSE)
          
          write.csv(hub_genes,
                    paste0("06_Advanced_Results/Network_Data/Hub_Genes_", target_module, ".csv"),
                    row.names = FALSE)
          
          cat(paste("      Network created:", igraph::vcount(g), "nodes,", 
                    igraph::ecount(g), "edges\n"))
        } else {
          cat("      Insufficient network connections, skipping network plot\n")
        }
        
      }, error = function(e) {
        cat(paste("      Error creating network:", e$message, "\n"))
      })
    } else {
      cat("      Too few genes, skipping network plot\n")
    }
    
    # 返回核心基因信息
    return(hub_genes)
    
  }, error = function(e) {
    cat(paste("      Error calculating kME:", e$message, "\n"))
    return(NULL)
  })
}

# 分析前3个最大的模块
module_sizes <- table(module_colors[module_colors != "grey"])
if (length(module_sizes) > 0) {
  top_modules <- names(sort(module_sizes, decreasing = TRUE))[
    1:min(3, length(module_sizes))
  ]
  
  hub_gene_results <- list()
  for (mod in top_modules) {
    result <- analyze_hub_genes(mod, datExpr, MEs, n_top = 20)
    if (!is.null(result)) {
      hub_gene_results[[mod]] <- result
    }
    gc()  # 每个模块处理后清理内存
  }
} else {
  cat("  Warning: No non-grey modules found\n")
  hub_gene_results <- list()
}

cat("  → Output files: Hub gene analysis results in '06_Advanced_Results/Figures/' and '06_Advanced_Results/Network_Data/'\n")
cat("    Purpose: Identifies and visualizes key regulatory genes within modules\n")

#### 7. 可视化：模块-模块特征基因相关性 ####
# ==============================================================================
cat("\n[Step 7] Generating module-module eigengene correlation heatmap...\n")

generate_module_correlation_heatmap <- function() {
  
  if (!exists("MEs") || ncol(MEs) < 3) {
    cat("  Warning: Insufficient modules for correlation analysis\n")
    return(NULL)
  }
  
  cat("  Preparing module eigengene data...\n")
  
  # 移除grey模块（如果存在）
  grey_cols <- grep("grey", colnames(MEs), ignore.case = TRUE)
  if (length(grey_cols) > 0) {
    MEs_clean <- MEs[, -grey_cols, drop = FALSE]
  } else {
    MEs_clean <- MEs
  }
  
  # 确保有足够的模块
  if (ncol(MEs_clean) < 3) {
    cat("  Warning: Too few modules after removing grey module\n")
    return(NULL)
  }
  
  # 计算相关性矩阵
  cat("  Calculating module-module correlations...\n")
  module_cor <- cor(MEs_clean, use = "pairwise.complete.obs")
  
  # ====== 修复：增强检查 ======
  # 检查相关性矩阵是否有有效值
  if (all(is.na(module_cor))) {
    cat("  ⚠️ 警告: 相关性矩阵全为NA值，跳过热图生成\n")
    cat("    可能原因:\n")
    cat("    1. 模块特征基因全为相同值（零方差）\n")
    cat("    2. 样本数太少（<2）\n")
    cat("    3. 数据格式异常\n")
    return(NULL)
  }
  
  # 提取所有非NA的相关性值
  cor_values <- as.vector(module_cor)
  cor_values <- cor_values[!is.na(cor_values)]
  
  # 检查是否有有效值
  if (length(cor_values) == 0) {
    cat("  ⚠️ 警告: 无有效相关性值，跳过热图生成\n")
    return(NULL)
  }
  
  # 检查是否所有值都相同（会导致热图颜色映射失败）
  if (var(cor_values) == 0) {
    cat("  ℹ️ 提示: 所有模块相关性值相同，创建简化热图\n")
    
    # 创建简化的热图（无颜色梯度）
    cor_df <- as.data.frame(module_cor)
    cor_df$Module <- rownames(cor_df)
    
    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = 1, y = Module)) +
      ggplot2::geom_tile(fill = "gray80") +
      ggplot2::geom_text(ggplot2::aes(label = round(module_cor[,1], 3)), 
                         size = 3) +
      ggplot2::labs(title = "模块相关性（所有值相同）",
                    x = "", y = "模块") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    
    ggplot2::ggsave("06_Advanced_Results/Figures/4_Module_Module_Correlation.pdf",
                    p, width = 8, height = 10)
    
    cat("  → 已保存简化版模块相关性图\n")
    return(NULL)
  }
  
  # ====== 原始代码继续（但有保护） ======
  
  # 确保范围在-1到1之间
  cor_min <- max(-1, min(cor_values, na.rm = TRUE))
  cor_max <- min(1, max(cor_values, na.rm = TRUE))
  
  # 如果所有值都相同，调整范围
  if (cor_min == cor_max) {
    cor_min <- -1
    cor_max <- 1
  }
  
  # 创建固定的breaks序列
  breaks_seq <- seq(cor_min, cor_max, length.out = 100)
  
  tryCatch({
    pheatmap::pheatmap(module_cor,
                       main = "Module-Module Eigengene Correlation",
                       subtitle = "Similarity between co-expression modules (non-trait analysis)",
                       color = color_palette,
                       breaks = breaks_seq,
                       display_numbers = TRUE,
                       number_format = "%.2f",
                       number_color = ifelse(abs(module_cor) > 0.7, "white", "black"),
                       fontsize_number = 6,
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       annotation_row = annotation_row,
                       annotation_col = annotation_row,
                       annotation_colors = list(
                         ModuleSize = grDevices::colorRampPalette(c("white", "darkred"))(100)
                       ),
                       fontsize_row = 8,
                       fontsize_col = 8,
                       border_color = NA,
                       cellwidth = 18,
                       cellheight = 18,
                       filename = "06_Advanced_Results/Figures/4_Module_Module_Correlation.pdf",
                       width = 12,
                       height = 11)
    
    # 创建PNG版本
    grDevices::png("06_Advanced_Results/Figures/4_Module_Module_Correlation.png",
                   width = 1100, height = 1000, res = 150)
    
    pheatmap::pheatmap(module_cor,
                       main = "Module-Module Eigengene Correlation",
                       color = color_palette,
                       breaks = breaks_seq,
                       display_numbers = FALSE,
                       cluster_rows = TRUE,
                       cluster_cols = TRUE)
    
    grDevices::dev.off()
    
    # 保存相关性矩阵
    write.csv(module_cor,
              "06_Advanced_Results/Tables/Module_Correlation_Matrix.csv",
              row.names = TRUE)
    
    cat("  → Output files: Figures/4_Module_Module_Correlation.pdf/.png\n")
    cat("    Purpose: Visualizes similarity between modules\n")
    cat("             Helps identify potentially redundant modules for merging\n")
    cat("    Note: This shows internal module structure, not trait associations\n")
    
  }, error = function(e) {
    cat(paste("  Error creating heatmap:", e$message, "\n"))
    cat("  Attempting alternative heatmap creation...\n")
    
    # 备用方案：使用ggplot2创建热图
    try({
      library(reshape2)
      cor_melt <- reshape2::melt(module_cor)
      colnames(cor_melt) <- c("Module1", "Module2", "Correlation")
      
      p <- ggplot2::ggplot(cor_melt, ggplot2::aes(x = Module1, y = Module2, fill = Correlation)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                                      midpoint = 0, limits = c(-1, 1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Module-Module Eigengene Correlation",
                      x = "Module", y = "Module")
      
      ggplot2::ggsave("06_Advanced_Results/Figures/4_Module_Module_Correlation_ggplot.pdf",
                      p, width = 10, height = 9)
      cat("  Created alternative heatmap using ggplot2\n")
    })
  })
}

generate_module_correlation_heatmap()
gc()

#### 8. 总结和导出 ####
# ==============================================================================
cat("\n[Step 8] Generating analysis summary and exporting final results...\n")

generate_summary_report <- function() {
  
  cat("  Compiling analysis summary...\n")
  
  # 模块统计
  module_table <- table(module_colors)
  module_stats <- data.frame(
    Module = names(module_table),
    GeneCount = as.numeric(module_table),
    stringsAsFactors = FALSE
  )
  
  # 按大小排序（排除grey）
  module_stats_no_grey <- module_stats[module_stats$Module != "grey", ]
  if (nrow(module_stats_no_grey) > 0) {
    module_stats_no_grey <- module_stats_no_grey[order(-module_stats_no_grey$GeneCount), ]
  }
  
  # 核心基因总结
  hub_summary <- data.frame()
  if (length(hub_gene_results) > 0) {
    for (mod in names(hub_gene_results)) {
      if (!is.null(hub_gene_results[[mod]]) && nrow(hub_gene_results[[mod]]) > 0) {
        top_hub <- head(hub_gene_results[[mod]][order(-abs(hub_gene_results[[mod]]$kME)), ], 3)
        hub_summary <- rbind(hub_summary,
                             data.frame(Module = mod,
                                        HubGenes = paste(top_hub$Gene, collapse = ", ")))
      }
    }
  }
  
  # 创建总结文本 - 包含Medicago sativa特定信息和学术诚信提醒
  summary_text <- c(
    "=================================================================================",
    "WGCNA ADVANCED ANALYSIS SUMMARY REPORT",
    paste("Species: Medicago sativa (Alfalfa)"),
    paste("Generated:", Sys.time()),
    "=================================================================================",
    "",
    "ACADEMIC INTEGRITY NOTICE:",
    "  IMPORTANT: No randomly generated functional enrichment plots were created.",
    "  This analysis only prepares gene lists for EXTERNAL functional enrichment tools.",
    "  All visualizations in this report are based on REAL experimental data.",
    "",
    "1. DATA OVERVIEW",
    paste("   Total genes analyzed:", ncol(datExpr)),
    paste("   Total samples:", nrow(datExpr)),
    paste("   Network soft threshold (power):", softPower),
    "",
    "2. MODULE DETECTION",
    paste("   Total modules detected:", length(unique(module_colors))),
    paste("   Non-grey modules:", sum(names(module_table) != "grey")),
    "",
    "3. MAIN MODULES (sorted by size)",
    if (nrow(module_stats_no_grey) > 0) {
      paste("   ", apply(head(module_stats_no_grey, 10), 1, 
                         function(x) paste(x[1], ": ", x[2], " genes", sep = "")), 
            collapse = "\n   ")
    } else {
      "   No non-grey modules detected"
    },
    "",
    "4. HUB GENE SUMMARY",
    if (nrow(hub_summary) > 0) {
      paste("   ", apply(hub_summary, 1, 
                         function(x) paste(x[1], ": ", x[2], sep = "")), 
            collapse = "\n   ")
    } else {
      "   No hub gene analysis performed or available"
    },
    "",
    "5. OUTPUT FILES",
    "   06_Advanced_Results/Figures/ (REAL DATA VISUALIZATIONS):",
    "     - 1_Gene_Dendrogram_and_Modules.pdf: Gene clustering and module assignment",
    "     - 2_Network_Topology_Analysis.pdf: Network topology validation",
    "     - 3_Heatmap_Hub_*.pdf: Hub gene expression heatmaps (real data)",
    "     - 3_Network_Hub_*.pdf: Hub gene interaction networks (dynamic thresholding)",
    "     - 4_Module_Module_Correlation.pdf: Inter-module similarity",
    "",
    "   06_Advanced_Results/Enrichment/ (PREPARATION FOR EXTERNAL ANALYSIS):",
    "     - Genes_*.txt: Gene lists for each module (upload to external tools)",
    "     - Enrichment_Analysis_Guide.txt: Detailed instructions for external analysis",
    "     - Visualization_Template.R: R script template for visualizing real results",
    "",
    "   06_Advanced_Results/Enrichment/Figures/ (REAL DATA):",
    "     - Expression_Pattern_*.pdf: Expression heatmaps of module genes",
    "     - [EMPTY - Use external tools then run Visualization_Template.R]",
    "",
    "   06_Advanced_Results/Tables/",
    "     - Module_Correlation_Matrix.csv: Module-module correlation matrix",
    "",
    "   06_Advanced_Results/Network_Data/",
    "     - Hub_Genes_*.csv: Hub gene lists",
    "     - Network_Edges_*.csv: Network edge lists (dynamic thresholding applied)",
    "",
    "6. CRITICAL NEXT STEPS FOR FUNCTIONAL ANALYSIS:",
    "   1. Upload Genes_*.txt files from Enrichment/ directory to:",
    "      - AgriGO v2.0 (http://systemsbiology.cau.edu.cn/agriGOv2/)",
    "      - EggNOG-mapper (http://eggnog-mapper.embl.de/)",
    "      - KEGG Mapper with Medicago truncatula (code: mtr)",
    "   2. Download enrichment results from external tools",
    "   3. Run Visualization_Template.R to create publication-quality figures",
    "   4. Do NOT use any randomly generated functional category plots",
    "",
    "7. TECHNICAL IMPROVEMENTS IN THIS VERSION:",
    "   - Dynamic network thresholding (no fixed 0.6 cutoff)",
    "   - Academic integrity: No random functional categories generated",
    "   - Complete bracket-matching error fixes",
    "   - Memory optimization for 16GB RAM systems",
    "",
    "=================================================================================",
    "ANALYSIS COMPLETE - PROCEED TO EXTERNAL FUNCTIONAL ENRICHMENT"
  )
  
  # 将总结写入文件
  writeLines(summary_text, "06_Advanced_Results/WGCNA_Advanced_Analysis_Summary.txt")
  
  # 同时输出到控制台
  cat(paste(summary_text, collapse = "\n"))
  cat("\n\n")
}

generate_summary_report()

# 最终内存清理
cat("  Performing final memory cleanup...\n")
rm(list = ls(pattern = "^temp_"))
gc()

monitor_memory("Final state")

cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║               ANALYSIS COMPLETE!                  ║\n")
cat("║  Publication-ready results saved to:              ║\n")
cat("║  06_Advanced_Results/                             ║\n")
cat("╚═══════════════════════════════════════════════════╝\n\n")

cat("SUMMARY OF KEY OUTPUTS:\n")
cat("1. GENE CLUSTERING: Figures/1_Gene_Dendrogram_and_Modules.pdf\n")
cat("   - Shows hierarchical clustering and module assignment of all genes\n\n")
cat("2. NETWORK VALIDATION: Figures/2_Network_Topology_Analysis.pdf\n")
cat("   - Validates scale-free topology assumption for selected power\n\n")
cat("3. HUB GENE ANALYSIS: Figures/3_Heatmap_Hub_* and 3_Network_Hub_*.pdf\n")
cat("   - Heatmaps show expression patterns of top hub genes (REAL DATA)\n")
cat("   - Network graphs with DYNAMIC THRESHOLDING (no fixed 0.6 cutoff)\n\n")
cat("4. MODULE CORRELATION: Figures/4_Module_Module_Correlation.pdf\n")
cat("   - Shows similarity between module eigengenes\n")
cat("   - Used to assess module independence (NOT trait association)\n\n")
cat("5. FUNCTIONAL ENRICHMENT PREPARATION (ACADEMIC INTEGRITY PROTECTED):\n")
cat("   - Enrichment/Genes_*.txt: REAL gene lists for external tools\n")
cat("   - Enrichment/Figures/Expression_Pattern_*.pdf: REAL expression heatmaps\n")
cat("   - NO RANDOMLY GENERATED FUNCTIONAL CATEGORY PLOTS\n")
cat("   - Use external tools then run Visualization_Template.R\n\n")
cat("6. ACADEMIC INTEGRITY REMINDER:\n")
cat("   - This script does NOT create fake functional enrichment results\n")
cat("   - All visualizations are based on real experimental data\n")
cat("   - Upload Genes_*.txt to external tools for proper enrichment analysis\n\n")
cat("Note: All figures are optimized for publication with consistent styling.\n")
cat("      Memory usage was optimized for 16GB RAM systems throughout.\n")

# 记录会话信息
sink("06_Advanced_Results/R_Session_Info.txt")
sessionInfo()
sink()

