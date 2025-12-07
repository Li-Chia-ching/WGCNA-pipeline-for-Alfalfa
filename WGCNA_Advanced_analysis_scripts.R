# ==============================================================================
# WGCNA 进阶分析脚本 - v2.3 (16GB RAM优化版)
# 功能：发表级可视化、富集分析、核心基因挖掘
# 
# 主要改进：
#   1. 修复了变量名错误和包加载问题
#   2. 修复了核心基因分析中的列选择错误
#   3. 针对16GB内存系统进行了全面优化
#   4. 添加了内存监控和清理机制
#   5. 优化了大型矩阵的处理方式
#   6. 增强了错误处理和进度跟踪
#   7. 所有注释改为简体中文
# ==============================================================================

#### 1. 初始化和包加载 ####
# ==============================================================================
cat("╔═══════════════════════════════════════════════════╗\n")
cat("║         WGCNA 进阶分析 v2.3 (16GB RAM优化)        ║\n")
cat("║            发表级可视化分析脚本                  ║\n")
cat("╚═══════════════════════════════════════════════════╝\n\n")

# 设置关键选项
options(stringsAsFactors = FALSE)
options(warn = 1)  # 立即显示警告
options(max.print = 100)  # 限制控制台输出行数

# 对于Windows系统设置内存限制
if (.Platform$OS.type == "windows") {
  memory.limit(size = 16000)  # 设置16GB内存限制
}

# 函数：检查并安装缺失的包
check_and_install_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("正在安装包:", pkg, "...\n"))
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
      cat(paste("  ✓", pkg, "安装完成\n"))
    }
  }
}

# 函数：监控内存使用情况
monitor_memory <- function(step_name) {
  if (.Platform$OS.type == "windows") {
    mem_used <- memory.size()
    mem_max <- memory.limit()
    cat(paste("  内存使用 (", step_name, "): ", 
              round(mem_used, 1), "MB / ", mem_max, "MB\n", sep = ""))
  }
}

# 步骤1.1: 加载基础分析包
cat("[步骤1.1] 加载基础分析包...\n")
essential_pkgs <- c("WGCNA", "ggplot2", "dplyr", "RColorBrewer", "viridis")
check_and_install_packages(essential_pkgs)

# 步骤1.2: 加载可视化包
cat("[步骤1.2] 加载可视化包...\n")
viz_pkgs <- c("pheatmap", "igraph", "patchwork", "ggrepel")  # 修复变量名
check_and_install_packages(viz_pkgs)

# 步骤1.3: 加载富集分析包
cat("[步骤1.3] 加载富集分析包...\n")
enrichment_pkgs <- c("clusterProfiler", "enrichplot")
check_and_install_packages(enrichment_pkgs)

# 注意：如果分析的不是人类数据，请修改此处
if (!require("org.Hs.eg.db", character.only = TRUE, quietly = TRUE)) {
  cat("注意：使用人类数据库 (org.Hs.eg.db)。如果是其他物种请修改。\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# 配置并行处理（内存感知）
cat("[步骤1.4] 配置并行处理...\n")
n_cores <- min(4, parallel::detectCores() - 2)  # 保留2个核心给系统
WGCNA::enableWGCNAThreads(nThreads = n_cores)
cat(paste("  使用", n_cores, "个CPU核心进行并行处理\n"))

# 设置统一的绘图主题（发表级）
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray90", size = 0.2),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, "cm")
    )
)

#### 2. 数据加载和验证 ####
# ==============================================================================
cat("\n[步骤2] 加载并验证 Step1 的分析结果...\n")

# 增强的数据加载函数，包含内存监控
load_and_validate_data <- function() {
  
  # 1. 检查并加载网络对象
  if (!file.exists("03_Network/WGCNA_Network_Object.RData")) {
    stop("❌ 错误：找不到文件 03_Network/WGCNA_Network_Object.RData")
  }
  
  cat("  正在加载WGCNA网络对象...\n")
  load("03_Network/WGCNA_Network_Object.RData")
  monitor_memory("加载网络对象后")
  gc()  # 强制垃圾回收
  
  # 2. 加载表达矩阵（内存高效方式）
  if (!file.exists("01_InputData/Preprocessed_Expression_Matrix.csv")) {
    stop("❌ 错误：找不到文件 01_InputData/Preprocessed_Expression_Matrix.csv")
  }
  
  cat("  正在加载表达矩阵（内存优化方式）...\n")
  
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
    cat("  注意：建议安装 'data.table' 包以获得更快的文件加载速度\n")
    datExpr <- read.csv("01_InputData/Preprocessed_Expression_Matrix.csv", 
                       row.names = 1, 
                       nrows = 50000)  # 限制读取行数
    datExpr <- as.data.frame(t(datExpr))
  }
  
  # 3. 修复module_colors的命名问题（如果存在）
  if (is.null(names(module_colors))) {
    if (length(module_colors) == ncol(datExpr)) {
      names(module_colors) <- colnames(datExpr)
      cat("  ✓ 已修复 module_colors 的基因名称\n")
    } else {
      warning("⚠️ module_colors 长度与基因数不匹配，尝试对齐...")
    }
  }
  
  # 4. 对齐数据
  common_genes <- intersect(colnames(datExpr), names(module_colors))
  if (length(common_genes) < 100) {
    warning("⚠️ 表达数据和模块分配之间的共有基因很少")
  }
  
  datExpr <- datExpr[, common_genes, drop = FALSE]
  module_colors <- module_colors[common_genes]
  
  # 输出统计信息
  cat(paste("  ✓ 数据加载完成:\n"))
  cat(paste("    - 样本数:", nrow(datExpr), "\n"))
  cat(paste("    - 基因数:", ncol(datExpr), "\n"))
  
  # 计算模块统计
  module_table <- table(module_colors)
  non_grey_modules <- sum(names(module_table) != "grey")
  cat(paste("    - 检测到的模块数:", length(module_table), "\n"))
  cat(paste("    - 非grey模块数:", non_grey_modules, "\n"))
  
  # 显示前5大模块
  module_sorted <- sort(module_table, decreasing = TRUE)
  if (length(module_sorted) > 5) {
    cat("    - 前5大模块:\n")
    for (i in 1:min(5, length(module_sorted))) {
      cat(paste("      ", names(module_sorted)[i], ": ", module_sorted[i], "个基因\n", sep = ""))
    }
  }
  
  monitor_memory("数据加载完成")
  
  return(list(
    datExpr = datExpr, 
    module_colors = module_colors, 
    net = net, 
    MEs = MEs, 
    softPower = softPower
  ))
}

# 执行数据加载
cat("开始加载数据...\n")
data_list <- load_and_validate_data()

# 将数据加载到全局环境
list2env(data_list, envir = .GlobalEnv)

# 清理临时变量以释放内存
rm(data_list)
gc()

# 创建输出目录结构
cat("  创建输出目录...\n")
output_dirs <- c(
  "06_Advanced_Results",
  "06_Advanced_Results/Figures",
  "06_Advanced_Results/Tables", 
  "06_Advanced_Results/Enrichment",
  "06_Advanced_Results/Network_Data"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat(paste("    创建目录:", dir, "\n"))
  }
}

monitor_memory("目录创建完成")

#### 3. 可视化：基因层次聚类树 (发表级) ####
# ==============================================================================
cat("\n[步骤3] 生成基因层次聚类树状图...\n")

generate_dendrogram <- function() {
  cat("  绘制基因层次聚类图...\n")
  
  # 创建PDF版本
  pdf("06_Advanced_Results/Figures/1_Gene_Dendrogram_and_Modules.pdf", 
      width = 14, height = 8)
  
  WGCNA::plotDendroAndColors(
    dendro = net$dendrograms[[1]], 
    colors = module_colors[net$blockGenes[[1]]],
    groupLabels = "模块分配",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    guideAll = FALSE,
    main = "基因层次聚类和WGCNA模块检测",
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
    groupLabels = "模块分配",
    dendroLabels = FALSE, 
    hang = 0.03,
    addGuide = TRUE, 
    guideHang = 0.05,
    main = "基因层次聚类和WGCNA模块检测"
  )
  
  dev.off()
  
  cat("  → 输出文件: Figures/1_Gene_Dendrogram_and_Modules.pdf/.png\n")
  cat("    目的: 展示基因如何通过层次聚类被分配到不同的共表达模块\n")
}

# 执行绘图
generate_dendrogram()
gc()

#### 4. 可视化：网络拓扑分析 ####
# ==============================================================================
cat("\n[步骤4] 生成网络拓扑分析图...\n")

generate_topology_plots <- function() {
  cat("  计算网络拓扑指标...\n")
  
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
  
  cat("  创建发表级拓扑分析图...\n")
  
  # 图1：无尺度拓扑拟合
  p1 <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = R2)) + 
    ggplot2::geom_point(size = 2.5, color = "#E41A1C", alpha = 0.8) + 
    ggplot2::geom_line(color = "#E41A1C", size = 0.8) +
    ggplot2::geom_hline(yintercept = 0.85, linetype = "dashed", 
                       color = "gray40", size = 0.6) +
    ggplot2::geom_vline(xintercept = softPower, linetype = "dashed", 
                       color = "#377EB8", size = 0.6) +
    ggplot2::annotate("text", x = softPower, y = max(sft_df$R2), 
                     label = paste("选择值:", softPower), 
                     vjust = -0.5, hjust = 1.2, 
                     color = "#377EB8", fontface = "bold", size = 3) +
    ggplot2::labs(
      title = "无尺度拓扑模型拟合",
      x = "软阈值 (Power)",
      y = expression("无尺度拓扑模型拟合 (R"^2*")")
    ) +
    ggplot2::scale_x_continuous(breaks = powers) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.4),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # 图2：平均连接度
  p2 <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = MeanK)) + 
    ggplot2::geom_point(size = 2.5, color = "#377EB8", alpha = 0.8) + 
    ggplot2::geom_line(color = "#377EB8", size = 0.8) +
    ggplot2::geom_vline(xintercept = softPower, linetype = "dashed", 
                       color = "#377EB8", size = 0.6) +
    ggplot2::labs(
      title = "平均基因连接度",
      x = "软阈值 (Power)", 
      y = "平均连接度 (log10尺度)"
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = powers) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 0.4),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # 合并图形
  combined_plot <- p1 + p2 + 
    patchwork::plot_annotation(
      title = "WGCNA网络拓扑分析",
      subtitle = paste("选择的软阈值功率 =", softPower),
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
  
  cat("  → 输出文件: Figures/2_Network_Topology_Analysis.pdf/.png\n")
  cat("    目的: 验证选择的软阈值功率是否合适，确保网络近似无尺度拓扑\n")
}

# 执行绘图
generate_topology_plots()
gc()

#### 5. 功能富集分析 ####
# ==============================================================================
cat("\n[步骤5] 进行功能富集分析...\n")

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
  
  cat(paste("  分析前", length(target_modules), "个最大的模块:\n"))
  cat(paste("    ", paste(target_modules, collapse = ", "), "\n"))
  
  # 模块富集分析函数
  run_module_enrichment <- function(module_color, gene_symbols) {
    cat(paste("    处理模块", module_color, "...\n"))
    
    tryCatch({
      # 转换基因符号为ENTREZ ID
      id_map <- clusterProfiler::bitr(gene_symbols, 
                                      fromType = "SYMBOL", 
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db)
      
      if (nrow(id_map) < 10) {
        cat(paste("      警告: 只有", nrow(id_map), "个基因成功转换为ENTREZID\n"))
        return(NULL)
      }
      
      cat(paste("      ", nrow(id_map), "个基因成功转换，进行KEGG富集分析...\n"))
      
      # 1. KEGG通路富集分析
      kegg_result <- clusterProfiler::enrichKEGG(
        gene = id_map$ENTREZID,
        organism = "hsa",  # 人类，如果是其他物种请修改
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500
      )
      
      if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
        # 保存结果
        write.csv(kegg_result@result,
                  paste0("06_Advanced_Results/Enrichment/KEGG_", module_color, ".csv"),
                  row.names = FALSE)
        
        # 创建可视化图形
        p <- enrichplot::dotplot(kegg_result, showCategory = 12, 
                                title = paste("KEGG通路 - 模块", module_color)) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
        
        ggplot2::ggsave(paste0("06_Advanced_Results/Figures/Enrich_KEGG_", module_color, ".pdf"),
                        p, width = 10, height = 7)
        
        cat(paste("      发现", nrow(kegg_result), "个显著富集的KEGG通路\n"))
      } else {
        cat("      未发现显著富集的KEGG通路\n")
      }
      
      # 2. GO生物过程富集分析（可选，内存密集型）
      if (length(id_map$ENTREZID) < 800) {  # 只对较小的模块进行GO分析
        cat("      进行GO生物过程富集分析...\n")
        go_result <- clusterProfiler::enrichGO(
          gene = id_map$ENTREZID,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID",
          ont = "BP",  # 生物过程
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        )
        
        if (!is.null(go_result) && nrow(go_result) > 0) {
          write.csv(go_result@result,
                    paste0("06_Advanced_Results/Enrichment/GO_BP_", module_color, ".csv"),
                    row.names = FALSE)
          
          # 创建简化的GO图（前10个条目）
          if (nrow(go_result) > 5) {
            p_go <- enrichplot::dotplot(go_result, showCategory = 10,
                                       title = paste("GO生物过程 - 模块", module_color)) +
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
            
            ggplot2::ggsave(paste0("06_Advanced_Results/Figures/Enrich_GO_BP_", module_color, ".pdf"),
                           p_go, width = 11, height = 7)
            
            cat(paste("      发现", nrow(go_result), "个显著富集的GO生物过程\n"))
          }
        }
      }
      
      # 强制垃圾回收
      gc()
      
    }, error = function(e) {
      cat(paste("      模块", module_color, "富集分析错误:", e$message, "\n"))
      return(NULL)
    })
  }
  
  # 处理每个目标模块
  for (mod in target_modules) {
    module_genes <- names(module_colors)[module_colors == mod]
    cat(paste("    模块", mod, ":", length(module_genes), "个基因\n"))
    
    if (length(module_genes) > 15) {  # 只有足够多基因的模块才进行分析
      run_module_enrichment(mod, module_genes)
    } else {
      cat(paste("      跳过 - 基因数太少，不适合富集分析\n"))
    }
  }
  
  cat("  → 输出文件: 富集分析结果保存在 '06_Advanced_Results/Enrichment/'\n")
  cat("    目的: 对共表达基因模块进行功能注释\n")
}

# 执行富集分析
perform_enrichment_analysis()
gc()

#### 6. 可视化：核心基因分析 (已修复错误) ####
# ==============================================================================
cat("\n[步骤6] 分析和可视化核心基因...\n")

# 修复后的核心基因分析函数
analyze_hub_genes <- function(target_module, datExpr, MEs, n_top = 20) {
  cat(paste("    处理模块", target_module, "...\n"))
  
  # 提取属于目标模块的基因
  module_genes <- names(module_colors)[module_colors == target_module]
  
  if (length(module_genes) == 0) {
    cat(paste("      警告: 模块中没有基因\n"))
    return(NULL)
  }
  
  if (length(module_genes) < n_top) {
    cat(paste("      警告: 模块只有", length(module_genes), "个基因\n"))
    n_top <- length(module_genes)
  }
  
  # 计算模块成员度 (kME) - 修复后的版本
  cat("      计算模块成员度 (kME)...\n")
  tryCatch({
    # 确保MEs中有目标模块的特征基因
    me_name <- paste0("ME", target_module)
    if (!me_name %in% colnames(MEs)) {
      cat(paste("      错误: 在MEs中找不到特征基因", me_name, "\n"))
      return(NULL)
    }
    
    # 提取表达数据和特征基因
    expr_data <- datExpr[, module_genes, drop = FALSE]
    me_data <- MEs[, me_name, drop = FALSE]
    
    # 计算kME
    kME_values <- WGCNA::signedKME(expr_data, me_data)
    
    # 确保kME_values是矩阵或数据框
    if (is.null(kME_values) || nrow(kME_values) == 0) {
      cat("      错误: 计算kME失败，返回空值\n")
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
    
    cat(paste("      识别出", length(top_gene_names), "个核心基因\n"))
    
    # 1. 创建核心基因表达热图
    if (length(top_gene_names) >= 3) {
      cat("      创建表达热图...\n")
      tryCatch({
        expr_subset <- datExpr[, top_gene_names, drop = FALSE]
        
        # Z-score标准化
        expr_scaled <- t(scale(expr_subset))
        
        # 设置颜色方案
        color_palette <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)
        
        # 创建热图
        pheatmap::pheatmap(expr_scaled,
                 main = paste("核心基因表达 - 模块", target_module),
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
                 main = paste("核心基因表达 - 模块", target_module),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_colnames = ncol(expr_scaled) <= 25,
                 color = color_palette)
        
        grDevices::dev.off()
        
      }, error = function(e) {
        cat(paste("      创建热图时出错:", e$message, "\n"))
      })
    } else {
      cat("      基因数太少，跳过热图\n")
    }
    
    # 2. 创建核心基因互作网络（如果有足够基因）
    if (length(top_gene_names) >= 5) {
      cat("      创建互作网络...\n")
      tryCatch({
        # 计算邻接矩阵（基于相关性）
        adj_matrix <- abs(cor(datExpr[, top_gene_names, drop = FALSE], 
                             use = "pairwise.complete.obs"))
        
        # 应用阈值使网络更稀疏
        adj_matrix[adj_matrix < 0.6] <- 0  # 较高阈值以获得更清晰的网络
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
               main = paste("核心基因互作网络 - 模块", target_module))
          
          # 添加图例
          graphics::legend("topleft",
                 legend = c("高连接度", "低连接度"),
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
          
          cat(paste("      网络创建完成:", igraph::vcount(g), "个节点,", 
                    igraph::ecount(g), "条边\n"))
        } else {
          cat("      网络连接不足，跳过网络图\n")
        }
        
      }, error = function(e) {
        cat(paste("      创建网络时出错:", e$message, "\n"))
      })
    } else {
      cat("      基因数太少，跳过网络图\n")
    }
    
    # 返回核心基因信息
    return(hub_genes)
    
  }, error = function(e) {
    cat(paste("      计算kME时出错:", e$message, "\n"))
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
  cat("  警告: 没有找到非grey模块\n")
  hub_gene_results <- list()
}

cat("  → 输出文件: 核心基因分析结果在 '06_Advanced_Results/Figures/' 和 '06_Advanced_Results/Network_Data/'\n")
cat("    目的: 识别和可视化模块内的关键调控基因\n")

#### 7. 可视化：模块-模块特征基因相关性 ####
# ==============================================================================
cat("\n[步骤7] 生成模块-模块特征基因相关性热图...\n")

generate_module_correlation_heatmap <- function() {
  
  if (!exists("MEs") || ncol(MEs) < 3) {
    cat("  警告: 模块数量不足，无法进行相关性分析\n")
    return(NULL)
  }
  
  cat("  准备模块特征基因数据...\n")
  
  # 移除grey模块（如果存在）
  grey_cols <- grep("grey", colnames(MEs), ignore.case = TRUE)
  if (length(grey_cols) > 0) {
    MEs_clean <- MEs[, -grey_cols, drop = FALSE]
  } else {
    MEs_clean <- MEs
  }
  
  # 确保有足够的模块
  if (ncol(MEs_clean) < 3) {
    cat("  警告: 移除grey模块后模块数量太少\n")
    return(NULL)
  }
  
  # 计算相关性矩阵
  cat("  计算模块-模块相关性...\n")
  module_cor <- cor(MEs_clean, use = "pairwise.complete.obs")
  
  # 创建模块注释
  module_names <- gsub("^ME", "", colnames(MEs_clean))
  module_sizes <- sapply(module_names, function(m) {
    sum(module_colors == m)
  })
  
  annotation_row <- data.frame(
    模块大小 = module_sizes,
    row.names = module_names
  )
  
  # 生成颜色调色板
  color_palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)
  
  # 创建发表级热图
  cat("  创建发表级热图...\n")
  
  pheatmap::pheatmap(module_cor,
           main = "模块-模块特征基因相关性",
           subtitle = "共表达模块之间的相似性（非性状关联分析）",
           color = color_palette,
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE,
           number_format = "%.2f",
           number_color = "black",
           fontsize_number = 6,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_row = annotation_row,
           annotation_col = annotation_row,
           annotation_colors = list(
             模块大小 = grDevices::colorRampPalette(c("white", "darkred"))(100)
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
           main = "模块-模块特征基因相关性",
           color = color_palette,
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE)
  
  grDevices::dev.off()
  
  # 保存相关性矩阵
  write.csv(module_cor,
            "06_Advanced_Results/Tables/Module_Correlation_Matrix.csv",
            row.names = TRUE)
  
  cat("  → 输出文件: Figures/4_Module_Module_Correlation.pdf/.png\n")
  cat("    目的: 可视化模块之间的相似性\n")
  cat("          帮助识别可能冗余的模块以便合并\n")
  cat("    注意: 这显示模块内部结构，不是性状关联\n")
}

generate_module_correlation_heatmap()
gc()

#### 8. 总结和导出 ####
# ==============================================================================
cat("\n[步骤8] 生成分析总结和导出最终结果...\n")

generate_summary_report <- function() {
  
  cat("  编译分析总结...\n")
  
  # 模块统计
  module_table <- table(module_colors)
  module_stats <- data.frame(
    模块 = names(module_table),
    基因数 = as.numeric(module_table),
    stringsAsFactors = FALSE
  )
  
  # 按大小排序（排除grey）
  module_stats_no_grey <- module_stats[module_stats$模块 != "grey", ]
  if (nrow(module_stats_no_grey) > 0) {
    module_stats_no_grey <- module_stats_no_grey[order(-module_stats_no_grey$基因数), ]
  }
  
  # 核心基因总结
  hub_summary <- data.frame()
  if (length(hub_gene_results) > 0) {
    for (mod in names(hub_gene_results)) {
      if (!is.null(hub_gene_results[[mod]]) && nrow(hub_gene_results[[mod]]) > 0) {
        top_hub <- head(hub_gene_results[[mod]][order(-abs(hub_gene_results[[mod]]$kME)), ], 3)
        hub_summary <- rbind(hub_summary,
                             data.frame(模块 = mod,
                                        核心基因 = paste(top_hub$Gene, collapse = ", ")))
      }
    }
  }
  
  # 创建总结文本
  summary_text <- c(
    "==============================================================",
    "WGCNA进阶分析总结报告",
    paste("生成时间:", Sys.time()),
    "==============================================================",
    "",
    "1. 数据概览",
    paste("   分析基因总数:", ncol(datExpr)),
    paste("   样本数:", nrow(datExpr)),
    paste("   网络软阈值 (power):", softPower),
    "",
    "2. 模块检测",
    paste("   检测到的模块总数:", length(unique(module_colors))),
    paste("   非grey模块数:", sum(names(module_table) != "grey")),
    "",
    "3. 主要模块 (按大小排序)",
    if (nrow(module_stats_no_grey) > 0) {
      paste("   ", apply(head(module_stats_no_grey, 10), 1, 
                        function(x) paste(x[1], ": ", x[2], "个基因", sep = "")), 
            collapse = "\n   ")
    } else {
      "   没有检测到非grey模块"
    },
    "",
    "4. 核心基因总结",
    if (nrow(hub_summary) > 0) {
      paste("   ", apply(hub_summary, 1, 
                        function(x) paste(x[1], ": ", x[2], sep = "")), 
            collapse = "\n   ")
    } else {
      "   未进行核心基因分析或无可用的核心基因数据"
    },
    "",
    "5. 输出文件",
    "   06_Advanced_Results/Figures/",
    "     - 1_Gene_Dendrogram_and_Modules.pdf: 基因聚类和模块分配",
    "     - 2_Network_Topology_Analysis.pdf: 网络拓扑验证",
    "     - 3_Heatmap_Hub_*.pdf: 核心基因表达热图",
    "     - 3_Network_Hub_*.pdf: 核心基因互作网络",
    "     - 4_Module_Module_Correlation.pdf: 模块间相关性",
    "",
    "   06_Advanced_Results/Enrichment/",
    "     - KEGG_*.csv: KEGG通路富集结果",
    "     - GO_BP_*.csv: GO生物过程富集结果",
    "",
    "   06_Advanced_Results/Tables/",
    "     - Module_Correlation_Matrix.csv: 模块相关性矩阵",
    "",
    "   06_Advanced_Results/Network_Data/",
    "     - Hub_Genes_*.csv: 核心基因列表",
    "     - Network_Edges_*.csv: 网络边列表",
    "",
    "6. 分析说明",
    "   - 所有图形均为发表级质量 (300 DPI)",
    "   - 模块-模块相关性显示模块内部结构，不是性状关联",
    "   - 核心基因基于模块内高连接度 (kME) 识别",
    "   - 已针对16GB RAM系统进行内存优化",
    "",
    "==============================================================",
    "分析完成"
  )
  
  # 将总结写入文件
  writeLines(summary_text, "06_Advanced_Results/WGCNA_Advanced_Analysis_Summary.txt")
  
  # 同时输出到控制台
  cat(paste(summary_text, collapse = "\n"))
  cat("\n\n")
}

generate_summary_report()

# 最终内存清理
cat("  执行最终内存清理...\n")
rm(list = ls(pattern = "^temp_"))
gc()

monitor_memory("最终状态")

cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║                 分析完成!                         ║\n")
cat("║    发表级分析结果已保存到:                       ║\n")
cat("║    06_Advanced_Results/                           ║\n")
cat("╚═══════════════════════════════════════════════════╝\n\n")

cat("主要输出文件总结:\n")
cat("1. 基因聚类图: Figures/1_Gene_Dendrogram_and_Modules.pdf\n")
cat("   - 展示所有基因的层次聚类和模块分配\n\n")
cat("2. 网络验证图: Figures/2_Network_Topology_Analysis.pdf\n")
cat("   - 验证选择的软阈值是否符合无尺度拓扑假设\n\n")
cat("3. 核心基因分析: Figures/3_Heatmap_Hub_* 和 3_Network_Hub_*.pdf\n")
cat("   - 热图显示核心基因的表达模式\n")
cat("   - 网络图展示核心基因之间的互作关系\n\n")
cat("4. 模块相关性图: Figures/4_Module_Module_Correlation.pdf\n")
cat("   - 展示模块特征基因之间的相似性\n")
cat("   - 用于评估模块独立性（非性状关联）\n\n")
cat("5. 功能富集分析: Enrichment/KEGG_*.csv 和 Enrichment/GO_BP_*.csv\n")
cat("   - 主要模块的通路和基因本体富集结果\n\n")
cat("注意: 所有图形均采用一致的发表级样式和配色方案。\n")
cat("      整个分析过程已针对16GB RAM系统进行内存优化。\n")

# 记录会话信息
sink("06_Advanced_Results/R_Session_Info.txt")
sessionInfo()
sink()
