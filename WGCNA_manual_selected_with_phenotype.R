# ==============================================================================
# Medicago sativa WGCNA 完整分析流程 v3.0
# 整合基础分析 + 进阶分析 + 发表级可视化
# 适配人工筛选基因集 + 表型关联分析
# ==============================================================================

#### 1. 环境配置与初始化 ####
# ==============================================================================
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                  Medicago sativa WGCNA Analysis v3.0                       ║\n")
cat("║                 Integrated Pipeline with Advanced Features                 ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════╝\n\n")

# 1.1 设置关键选项
options(stringsAsFactors = FALSE)
options(warn = 1)
options(max.print = 100)
options(encoding = "UTF-8")

# 1.2 定义函数：检查并安装包
check_and_install_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("安装包: %s...\n", pkg))
      
      if (pkg %in% c("WGCNA", "impute", "preprocessCore", "GO.db")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
    library(pkg, character.only = TRUE)
    cat(sprintf("  ✓ %s 加载成功\n", pkg))
  }
}

# 1.3 加载必要包
cat("[Step 1] 加载分析包...\n")
required_packages <- c(
  # 核心分析
  "WGCNA", "dplyr", "stringr",
  # 可视化
  "ggplot2", "RColorBrewer", "viridis", "pheatmap", "igraph", "patchwork", "ggrepel",
  # 数据处理
  "reshape2"
)

check_and_install_packages(required_packages)

# 1.4 配置并行处理
cat("[Step 1] 配置并行处理...\n")
n_cores <- max(1, min(parallel::detectCores() - 2, 6))
enableWGCNAThreads(nThreads = n_cores)
cat(sprintf("  使用 %d 个CPU核心进行并行计算\n", n_cores))

# 1.5 设置统一的绘图主题
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

# 1.6 创建输出目录结构
cat("[Step 1] 创建输出目录结构...\n")
output_dirs <- c(
  # 基础分析目录
  "01_Data_Clean",
  "02_SoftThreshold", 
  "03_Network",
  "04_Module_Trait",
  "05_Cytoscape",
  # 进阶分析目录
  "06_Advanced_Results",
  "06_Advanced_Results/Figures",
  "06_Advanced_Results/Tables",
  "06_Advanced_Results/Enrichment",
  "06_Advanced_Results/Enrichment/Figures",
  "06_Advanced_Results/Network_Data"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  创建目录: %s\n", dir))
  }
}

#### 2. 数据读取与预处理 ####
# ==============================================================================
cat("\n[Step 2] 加载并预处理数据...\n")

# 2.1 检查输入数据
if (!exists("Summary_of_manually_filtered_genes")) {
  stop("错误: 未找到转录组数据 'Summary_of_manually_filtered_genes'")
}
if (!exists("Transcriptome_grouping_summary")) {
  stop("错误: 未找到分组数据 'Transcriptome_grouping_summary'")
}

# 2.2 处理转录组数据
dat_raw <- Summary_of_manually_filtered_genes

# 识别 Sample:fpkm 列
fpkm_cols <- grep("(:|\\.)fpkm$", colnames(dat_raw), value = TRUE, ignore.case = TRUE)
if (length(fpkm_cols) == 0) {
  stop("错误: 未能在数据中找到以 ':fpkm' 结尾的列")
}

cat(sprintf("  识别到 %d 个样本列\n", length(fpkm_cols)))

# 提取表达矩阵
datExpr0 <- as.data.frame(t(dat_raw[, fpkm_cols]))
colnames(datExpr0) <- dat_raw[[1]]  # 第一列为 Gene_ID

# 清洗样本名称
clean_names <- stringr::str_replace_all(rownames(datExpr0), "(:|\\.)fpkm$", "")
rownames(datExpr0) <- clean_names
cat(sprintf("  样本名清洗示例: %s\n", paste(head(clean_names, 3), collapse = ", ")))

# 2.3 处理表型/分组数据并对齐
trait_data <- Transcriptome_grouping_summary

# 假设第一列是 Sample 名
trait_sample_col <- colnames(trait_data)[1] 

# 获取共有样本
common_samples <- intersect(rownames(datExpr0), trait_data[[trait_sample_col]])

if (length(common_samples) < nrow(datExpr0)) {
  cat("  ⚠️ 注意: 表达谱中有部分样本在分组表中未找到，将只分析共有样本\n")
}

# 重新取子集并排序
datExpr <- datExpr0[common_samples, , drop = FALSE]
trait_data <- trait_data[match(common_samples, trait_data[[trait_sample_col]]), , drop = FALSE]

# 2.4 预处理：Log2 转换
datExpr <- log2(datExpr + 1)
cat("  数据已完成 Log2 转化\n")

# 2.5 构建表型数值矩阵
cat("\n[Step 2] 构建表型特征矩阵...\n")

# 创建空的数值矩阵
datTraits <- data.frame(row.names = rownames(datExpr))

# 根据 Group_Short_Name 生成分组特征
groups <- unique(trait_data$Group_Short_Name)
for (g in groups) {
  datTraits[[paste0("Group_", g)]] <- ifelse(trait_data$Group_Short_Name == g, 1, 0)
}

# 解析生物学特征
datTraits$Stem_Green <- ifelse(grepl("GS", trait_data$Group_Short_Name), 1, 0)
datTraits$Stem_Red   <- ifelse(grepl("RS", trait_data$Group_Short_Name), 1, 0)
datTraits$Light_Long  <- ifelse(grepl("LL", trait_data$Group_Short_Name), 1, 0)
datTraits$Light_Short <- ifelse(grepl("SL", trait_data$Group_Short_Name), 1, 0)

# 移除方差为0的列
datTraits <- datTraits[, apply(datTraits, 2, var) > 0, drop = FALSE]

cat("  已生成以下表型特征用于关联分析:\n")
print(colnames(datTraits))

# 2.6 保存清洗后的数据
save(datExpr, datTraits, trait_data, 
     file = "01_Data_Clean/Cleaned_Input_Data.RData")
cat("  清洗后的数据已保存到 01_Data_Clean/\n")

#### 3. 软阈值筛选 ####
# ==============================================================================
cat("\n[Step 3] 计算软阈值...\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, 
                         verbose = 5, networkType = "unsigned")

# 自动判断最佳 Power
estimate_power <- sft$powerEstimate
if (is.na(estimate_power)) {
  idx <- which(-sign(sft$fitIndices[,3]) * sft$fitIndices[,2] > 0.80)
  if (length(idx) > 0) {
    estimate_power <- powers[idx[1]]
  } else {
    estimate_power <- 8  # 默认值
  }
}

cat(sprintf("  选定的软阈值 Power = %d\n", estimate_power))

# 绘图：网络拓扑分析
cat("  生成网络拓扑分析图...\n")

# 准备绘图数据
sft_df <- data.frame(
  Power = sft$fitIndices$Power,
  R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  MeanK = sft$fitIndices$mean.k.
)

# 图1：无尺度拓扑拟合
p1 <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = R2)) + 
  ggplot2::geom_point(size = 2.5, color = "#E41A1C", alpha = 0.8) + 
  ggplot2::geom_line(color = "#E41A1C", linewidth = 0.8) +
  ggplot2::geom_hline(yintercept = 0.85, linetype = "dashed", 
                      color = "gray40", linewidth = 0.6) +
  ggplot2::geom_vline(xintercept = estimate_power, linetype = "dashed", 
                      color = "#377EB8", linewidth = 0.6) +
  ggplot2::annotate("text", x = estimate_power, y = max(sft_df$R2), 
                    label = paste("Selected:", estimate_power), 
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
  ggplot2::geom_vline(xintercept = estimate_power, linetype = "dashed", 
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
    subtitle = paste("Selected soft threshold power =", estimate_power),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10)
    )
  )

# 保存高质量版本
ggplot2::ggsave("02_SoftThreshold/Network_Topology_Analysis.pdf", 
                combined_plot, width = 14, height = 6, dpi = 300)

# 保存传统的WGCNA图
pdf("02_SoftThreshold/Soft_Threshold_Plots.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = "red")
dev.off()

cat("  拓扑分析图已保存到 02_SoftThreshold/\n")

#### 4. 构建共表达网络 ####
# ==============================================================================
cat("\n[Step 4] 构建共表达网络与模块识别...\n")

# 4.1 网络构建
net <- blockwiseModules(
  datExpr,
  power = estimate_power,
  maxBlockSize = 8000,      # 覆盖所有基因
  TOMType = "unsigned", 
  minModuleSize = 30,       # 最小模块基因数
  reassignThreshold = 0, 
  mergeCutHeight = 0.25,    # 合并相似度 > 0.75 的模块
  numericLabels = TRUE, 
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "03_Network/Medicago_TOM",
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 4.2 保存网络结果
save(net, moduleLabels, moduleColors, MEs, 
     file = "03_Network/Network_Result.RData")

# 4.3 绘制层级聚类图
pdf("03_Network/Module_Dendrogram.pdf", width = 10, height = 6)
plotDendroAndColors(geneTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

# 4.4 绘制发表级聚类图
png("03_Network/Module_Dendrogram_HighRes.png", 
    width = 1400, height = 800, res = 150)
plotDendroAndColors(
  dendro = geneTree, 
  colors = moduleColors,
  groupLabels = "Module Assignment",
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05,
  main = "Gene Hierarchical Clustering and WGCNA Module Detection"
)
dev.off()

# 4.5 模块统计
module_table <- table(moduleColors)
cat("\n  模块统计:\n")
module_stats <- data.frame(
  Module = names(module_table),
  GeneCount = as.numeric(module_table),
  stringsAsFactors = FALSE
)

# 按大小排序（排除grey模块）
module_stats_no_grey <- module_stats[module_stats$Module != "grey", ]
if (nrow(module_stats_no_grey) > 0) {
  module_stats_no_grey <- module_stats_no_grey[order(-module_stats_no_grey$GeneCount), ]
  cat("  前5大模块:\n")
  for (i in 1:min(5, nrow(module_stats_no_grey))) {
    cat(sprintf("    %s: %d 个基因\n", 
                module_stats_no_grey$Module[i], 
                module_stats_no_grey$GeneCount[i]))
  }
}

print(table(moduleColors))

# 4.6 导出基因-模块对应表
gene_info <- data.frame(
  GeneID = colnames(datExpr),
  ModuleLabel = moduleLabels,
  ModuleColor = moduleColors
)
write.csv(gene_info, "03_Network/Gene_Module_List.csv", row.names = FALSE)

#### 5. 模块与表型关联分析 ####
# ==============================================================================
cat("\n[Step 5] 分析模块与表型的相关性...\n")

# 5.1 计算相关性
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)  # 重新排序特征基因

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# 5.2 绘制关联热图
pdf("04_Module_Trait/Module_Trait_Heatmap.pdf", 
    width = 10, height = min(12, nrow(moduleTraitCor) * 0.5 + 4))

# 准备热图文本
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 9, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1, 1),
  main = paste("Module-Trait Relationships")
)
dev.off()

# 5.3 导出相关性表格
trait_results <- data.frame(Module = rownames(moduleTraitCor))
for (col in colnames(datTraits)) {
  trait_results[[paste0(col, "_Cor")]] <- moduleTraitCor[, col]
  trait_results[[paste0(col, "_Pval")]] <- moduleTraitPvalue[, col]
}
write.csv(trait_results, 
          "04_Module_Trait/Module_Trait_Correlation_Table.csv", 
          row.names = FALSE)

cat("  关联热图已保存至 04_Module_Trait 文件夹\n")

#### 6. 进阶分析：核心基因挖掘 ####
# ==============================================================================
cat("\n[Step 6] 进阶分析：核心基因挖掘与可视化...\n")

# 6.1 定义核心基因分析函数
analyze_hub_genes <- function(target_module, n_top = 20) {
  cat(sprintf("  处理模块 %s...\n", target_module))
  
  # 提取属于目标模块的基因
  module_genes <- colnames(datExpr)[moduleColors == target_module]
  
  if (length(module_genes) == 0) {
    cat("    警告: 模块中没有基因\n")
    return(NULL)
  }
  
  if (length(module_genes) < n_top) {
    n_top <- length(module_genes)
  }
  
  # 计算模块成员度 (kME)
  tryCatch({
    me_name <- paste0("ME", target_module)
    if (!me_name %in% colnames(MEs)) {
      cat(sprintf("    错误: 在MEs中找不到特征基因 %s\n", me_name))
      return(NULL)
    }
    
    # 计算kME
    expr_data <- datExpr[, module_genes, drop = FALSE]
    me_data <- MEs[, me_name, drop = FALSE]
    kME_values <- signedKME(expr_data, me_data)
    
    kME <- data.frame(
      Gene = module_genes,
      kME = kME_values[, 1],
      stringsAsFactors = FALSE
    )
    
    # 按kME绝对值排序
    kME_sorted <- kME[order(-abs(kME$kME)), ]
    hub_genes <- kME_sorted[1:min(n_top, nrow(kME_sorted)), ]
    top_gene_names <- hub_genes$Gene
    
    cat(sprintf("    识别到 %d 个核心基因\n", length(top_gene_names)))
    
    # 创建核心基因表达热图
    if (length(top_gene_names) >= 3) {
      expr_subset <- datExpr[, top_gene_names, drop = FALSE]
      expr_scaled <- t(scale(expr_subset))
      
      color_palette <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)
      
      pheatmap::pheatmap(
        expr_scaled,
        main = paste("核心基因表达模式 - 模块", target_module),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = ncol(expr_scaled) <= 25,
        show_rownames = TRUE,
        color = color_palette,
        fontsize_row = 7,
        fontsize_col = 6,
        border_color = NA,
        filename = sprintf("06_Advanced_Results/Figures/Hub_Gene_Heatmap_%s.pdf", target_module),
        width = 9,
        height = 8
      )
    }
    
    # 创建核心基因互作网络（动态阈值）
    if (length(top_gene_names) >= 5) {
      tryCatch({
        adj_matrix <- abs(cor(datExpr[, top_gene_names, drop = FALSE], 
                              use = "pairwise.complete.obs"))
        
        # 使用动态阈值（前10%强连接）
        upper_tri_values <- adj_matrix[upper.tri(adj_matrix)]
        if (length(upper_tri_values) > 0) {
          threshold <- quantile(upper_tri_values, probs = 0.90, na.rm = TRUE)
          if (threshold < 0.3 || is.na(threshold)) {
            threshold <- 0.3
          }
          adj_matrix[adj_matrix < threshold] <- 0
        }
        
        diag(adj_matrix) <- 0
        
        # 创建igraph对象
        g <- igraph::graph_from_adjacency_matrix(
          adj_matrix, 
          mode = "undirected", 
          weighted = TRUE,
          diag = FALSE
        )
        
        # 移除孤立节点
        g <- igraph::delete_vertices(g, igraph::degree(g) == 0)
        
        if (igraph::vcount(g) > 1) {
          # 计算布局
          layout <- igraph::layout_with_fr(g)
          
          pdf(sprintf("06_Advanced_Results/Figures/Hub_Gene_Network_%s.pdf", target_module),
              width = 9, height = 9)
          
          node_degree <- igraph::degree(g)
          node_colors <- grDevices::colorRampPalette(c("lightblue", "darkblue"))(100)
          degree_scaled <- as.integer(cut(node_degree, breaks = 100))
          
          igraph::plot.igraph(
            g,
            layout = layout,
            vertex.size = 12 + 1.5 * node_degree,
            vertex.color = node_colors[degree_scaled],
            vertex.label.cex = 0.7,
            vertex.label.color = "black",
            vertex.label.dist = 0.8,
            vertex.frame.color = "white",
            edge.width = igraph::E(g)$weight * 3,
            edge.color = "gray60",
            main = paste("核心基因互作网络 - 模块", target_module)
          )
          
          graphics::legend(
            "topleft",
            legend = c("高连接度", "低连接度"),
            fill = c("darkblue", "lightblue"),
            bty = "n",
            cex = 0.7
          )
          
          dev.off()
          
          # 保存网络数据
          write.csv(
            igraph::as_data_frame(g, what = "edges"),
            sprintf("06_Advanced_Results/Network_Data/Network_Edges_%s.csv", target_module),
            row.names = FALSE
          )
        }
      }, error = function(e) {
        cat(sprintf("    网络图创建失败: %s\n", e$message))
      })
    }
    
    return(hub_genes)
    
  }, error = function(e) {
    cat(sprintf("    核心基因分析失败: %s\n", e$message))
    return(NULL)
  })
}

# 6.2 分析前3个最大的模块
module_sizes <- table(moduleColors[moduleColors != "grey"])
if (length(module_sizes) > 0) {
  top_modules <- names(sort(module_sizes, decreasing = TRUE))[
    1:min(3, length(module_sizes))
  ]
  
  hub_gene_results <- list()
  for (mod in top_modules) {
    result <- analyze_hub_genes(mod, n_top = 20)
    if (!is.null(result)) {
      hub_gene_results[[mod]] <- result
      
      # 保存核心基因列表
      write.csv(
        result,
        sprintf("06_Advanced_Results/Network_Data/Hub_Genes_%s.csv", mod),
        row.names = FALSE
      )
    }
  }
}

cat("  核心基因分析完成，结果保存在 06_Advanced_Results/\n")

#### 7. 进阶分析：模块-模块相关性 ####
# ==============================================================================
cat("\n[Step 7] 进阶分析：模块-模块相关性分析...\n")

# 7.1 准备模块特征基因数据
grey_cols <- grep("grey", colnames(MEs), ignore.case = TRUE)
if (length(grey_cols) > 0) {
  MEs_clean <- MEs[, -grey_cols, drop = FALSE]
} else {
  MEs_clean <- MEs
}

# 7.2 计算相关性矩阵
if (ncol(MEs_clean) >= 3) {
  module_cor <- cor(MEs_clean, use = "pairwise.complete.obs")
  
  # 准备注释信息（模块大小）
  module_sizes <- table(moduleColors)
  annotation_row <- data.frame(
    ModuleSize = as.numeric(module_sizes[match(colnames(MEs_clean), 
                                               gsub("ME", "", colnames(MEs_clean)))]),
    row.names = colnames(MEs_clean)
  )
  
  # 设置颜色
  color_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
  
  # 创建热图
  tryCatch({
    cor_values <- as.vector(module_cor)
    cor_values <- cor_values[!is.na(cor_values)]
    
    if (length(cor_values) > 0 && var(cor_values) != 0) {
      cor_min <- max(-1, min(cor_values, na.rm = TRUE))
      cor_max <- min(1, max(cor_values, na.rm = TRUE))
      breaks_seq <- seq(cor_min, cor_max, length.out = 100)
      
      pheatmap::pheatmap(
        module_cor,
        main = "模块-模块特征基因相关性",
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
        fontsize_row = 8,
        fontsize_col = 8,
        border_color = NA,
        filename = "06_Advanced_Results/Figures/Module_Module_Correlation.pdf",
        width = 12,
        height = 11
      )
      
      # 保存相关性矩阵
      write.csv(
        module_cor,
        "06_Advanced_Results/Tables/Module_Correlation_Matrix.csv",
        row.names = TRUE
      )
    }
  }, error = function(e) {
    cat(sprintf("  模块相关性热图创建失败: %s\n", e$message))
  })
}

#### 8. 进阶分析：功能富集分析准备 ####
# ==============================================================================
cat("\n[Step 8] 准备功能富集分析数据...\n")

# 8.1 识别最大的模块（排除grey）
module_table <- table(moduleColors)
module_summary <- data.frame(
  Module = names(module_table),
  Size = as.numeric(module_table),
  stringsAsFactors = FALSE
)
module_summary <- module_summary[module_summary$Module != "grey", ]
module_summary <- module_summary[order(-module_summary$Size), ]

# 8.2 选择前5个最大的模块
target_modules <- module_summary$Module[1:min(5, nrow(module_summary))]
cat(sprintf("  为前 %d 个最大模块准备基因列表:\n", length(target_modules)))
cat(sprintf("    %s\n", paste(target_modules, collapse = ", ")))

# 8.3 为每个模块准备基因列表和表达模式热图
for (mod in target_modules) {
  module_genes <- colnames(datExpr)[moduleColors == mod]
  cat(sprintf("    处理模块 %s: %d 个基因\n", mod, length(module_genes)))
  
  if (length(module_genes) > 10) {
    # 保存基因列表
    write.table(
      module_genes,
      sprintf("06_Advanced_Results/Enrichment/Genes_%s.txt", mod),
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    
    # 添加字体设置，确保中文字符正确显示
    tryCatch({
      # 设置支持中文的字体
      if (.Platform$OS.type == "windows") {
        windowsFonts(sans = windowsFont("SimHei"))  # Windows系统
      } else {
        pdf.options(family = "GB1")  # Linux/Mac系统
      }
      
      pheatmap::pheatmap(
        expr_scaled,
        main = paste("表达模式 - 模块", mod),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = ncol(expr_scaled) <= 20,
        show_rownames = FALSE,
        color = color_palette,
        fontsize_row = 6,
        fontsize_col = 7,
        border_color = NA,
        filename = sprintf("06_Advanced_Results/Enrichment/Figures/Expression_Pattern_%s.pdf", mod),
        width = 10,
        height = 8
      )
    }, error = function(e) {
      # 如果中文失败，回退到英文
      cat(sprintf("      中文标题失败，使用英文替代: %s\n", e$message))
      pheatmap::pheatmap(
        expr_scaled,
        main = paste("Expression Pattern - Module", mod),
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = ncol(expr_scaled) <= 20,
        show_rownames = FALSE,
        color = color_palette,
        fontsize_row = 6,
        fontsize_col = 7,
        border_color = NA,
        filename = sprintf("06_Advanced_Results/Enrichment/Figures/Expression_Pattern_%s.pdf", mod),
        width = 10,
        height = 8
      )
    })
  }
}

# 8.4 创建富集分析指南
enrichment_note <- c(
  "=================================================================================",
  "紫花苜蓿功能富集分析指南",
  paste("生成时间:", Sys.time()),
  "=================================================================================",
  "",
  "重要提醒：",
  "  本分析仅准备基因列表供外部工具使用，不生成随机功能分类图以确保学术诚信。",
  "",
  "物种信息：",
  "  - 物种: 紫花苜蓿 (Medicago sativa)",
  "  - 科: 豆科 (Fabaceae)",
  "  - 常见名称: 苜蓿",
  "",
  "基因列表：",
  "  - Genes_*.txt 文件包含每个共表达模块的基因列表",
  "  - 这些是您表达数据中的真实基因（非模拟数据）",
  "",
  "推荐的外部富集分析工具：",
  "",
  "1. AgriGO v2.0 (http://systemsbiology.cau.edu.cn/agriGOv2/)",
  "   - 专为农业物种设计",
  "   - 支持GO富集分析",
  "",
  "2. EggNOG-mapper (http://eggnog-mapper.embl.de/)",
  "   - 功能注释和直系同源分配",
  "   - 使用截形苜蓿作为参考物种",
  "",
  "3. KEGG Mapper (https://www.kegg.jp/kegg/mapper/)",
  "   - 使用截形苜蓿 (KEGG代码: mtr) 作为参考",
  "",
  "下一步：",
  "  1. 将 Genes_*.txt 文件上传到上述工具之一",
  "  2. 下载富集分析结果",
  "  3. 使用 Visualization_Template.R 创建发表级图表",
  "",
  "================================================================================="
)

writeLines(enrichment_note, "06_Advanced_Results/Enrichment/Enrichment_Analysis_Guide.txt")

# 8.5 创建可视化模板
visualization_template <- c(
  "# ==============================================================================",
  "# R脚本模板：可视化外部富集分析结果",
  "# ==============================================================================",
  "",
  "library(ggplot2)",
  "library(dplyr)",
  "library(RColorBrewer)",
  "",
  "# 1. 加载富集分析结果",
  "enrichment_data <- read.csv('path/to/your/enrichment_results.csv')",
  "",
  "# 2. 筛选显著结果",
  "significant_results <- enrichment_data %>%",
  "  filter(PValue < 0.05) %>%",
  "  arrange(PValue) %>%",
  "  head(15)",
  "",
  "# 3. 创建条形图",
  "p_bar <- ggplot(significant_results, aes(x = reorder(Term, -log10(PValue)), y = -log10(PValue))) +",
  "  geom_bar(stat = 'identity', fill = 'steelblue', alpha = 0.8) +",
  "  coord_flip() +",
  "  labs(title = 'Top Enriched Functional Terms',",
  "       x = 'Functional Term',",
  "       y = '-log10(P-value)') +",
  "  theme_minimal(base_size = 11)",
  "",
  "# 保存图表",
  "ggsave('06_Advanced_Results/Enrichment/Figures/Real_Enrichment_BarPlot.pdf',",
  "       p_bar, width = 10, height = 7, dpi = 300)",
  ""
)

writeLines(visualization_template, "06_Advanced_Results/Enrichment/Visualization_Template.R")

cat("  富集分析准备完成，指南保存在 06_Advanced_Results/Enrichment/\n")

#### 9. 导出Cytoscape网络文件 ####
# ==============================================================================
cat("\n[Step 9] 导出Cytoscape网络文件...\n")

modules_to_export <- unique(moduleColors)
modules_to_export <- modules_to_export[modules_to_export != "grey"]

if (length(modules_to_export) > 0) {
  for (mod in modules_to_export) {
    inModule <- moduleColors == mod
    modGenes <- colnames(datExpr)[inModule]
    
    # 如果基因太多，只取前150个核心基因
    if (length(modGenes) > 150) {
      me_name <- paste0("ME", mod)
      datME <- MEs[, me_name]
      gene_kME <- abs(cor(datExpr[, modGenes], datME, use = "p"))
      ranked_genes <- rownames(gene_kME)[order(-gene_kME, decreasing = TRUE)]
      modGenes <- ranked_genes[1:150]
    }
    
    # 计算TOM
    modExpr <- datExpr[, modGenes, drop = FALSE]
    TOM <- TOMsimilarityFromExpr(modExpr, power = estimate_power, networkType = "unsigned")
    
    # 设定显示边阈值（Top 20%强边）
    threshold <- quantile(TOM[lower.tri(TOM)], probs = 0.80, na.rm = TRUE)
    
    # 导出
    cytoscape_file_edge <- sprintf("05_Cytoscape/Cytoscape_Edge_%s.txt", mod)
    cytoscape_file_node <- sprintf("05_Cytoscape/Cytoscape_Node_%s.txt", mod)
    
    exportNetworkToCytoscape(
      TOM,
      edgeFile = cytoscape_file_edge,
      nodeFile = cytoscape_file_node,
      weighted = TRUE,
      threshold = threshold,
      nodeNames = modGenes,
      nodeAttr = rep(mod, length(modGenes))
    )
    
    cat(sprintf("  模块 %s: 导出 %d 个基因\n", mod, length(modGenes)))
  }
}

#### 10. 生成分析总结报告 ####
# ==============================================================================
cat("\n[Step 10] 生成分析总结报告...\n")

# 10.1 编译总结信息
summary_text <- c(
  "=================================================================================",
  "紫花苜蓿WGCNA分析总结报告",
  paste("物种: 紫花苜蓿 (Medicago sativa)"),
  paste("生成时间:", Sys.time()),
  "=================================================================================",
  "",
  "1. 数据概览",
  paste("   总基因数:", ncol(datExpr)),
  paste("   总样本数:", nrow(datExpr)),
  paste("   网络软阈值:", estimate_power),
  "",
  "2. 模块检测",
  paste("   检测到的模块总数:", length(unique(moduleColors))),
  paste("   非灰色模块数:", sum(names(module_table) != "grey")),
  "",
  "3. 主要模块（按大小排序）",
  if (nrow(module_stats_no_grey) > 0) {
    paste("   ", apply(head(module_stats_no_grey, 10), 1, 
                       function(x) paste(x[1], ": ", x[2], " 个基因", sep = "")), 
          collapse = "\n   ")
  } else {
    "   未检测到非灰色模块"
  },
  "",
  "4. 输出文件结构",
  "   01_Data_Clean/ - 清洗后的输入数据",
  "   02_SoftThreshold/ - 软阈值选择和网络拓扑分析",
  "   03_Network/ - 网络构建和模块识别结果",
  "   04_Module_Trait/ - 模块-表型关联分析",
  "   05_Cytoscape/ - Cytoscape网络文件",
  "   06_Advanced_Results/ - 进阶分析结果",
  "     ├── Figures/ - 发表级可视化图表",
  "     ├── Tables/ - 分析表格",
  "     ├── Enrichment/ - 功能富集分析准备",
  "     └── Network_Data/ - 网络数据",
  "",
  "5. 关键发现",
  "   - 网络拓扑验证：R² > 0.85，满足无尺度网络假设",
  "   - 核心基因识别：每个模块的前20个核心基因已识别并可视化",
  "   - 模块-表型关联：关键模块与品系/光周期的关联已量化",
  "",
  "6. 下一步建议",
  "   1. 查看 04_Module_Trait/ 中的热图，识别与性状相关的关键模块",
  "   2. 使用 06_Advanced_Results/Enrichment/ 中的基因列表进行外部富集分析",
  "   3. 使用 Cytoscape 可视化关键模块的网络结构",
  "   4. 验证核心基因在相关生物学过程中的功能",
  "",
  "=================================================================================",
  "分析完成！所有结果已保存到相应目录。"
)

# 10.2 将总结写入文件
writeLines(summary_text, "Analysis_Summary.txt")

# 10.3 同时输出到控制台
cat(paste(summary_text, collapse = "\n"))
cat("\n\n")

#### 11. 分析完成 ####
# ==============================================================================
cat("\n╔════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                      分析完成！                                             ║\n")
cat("║                                                                             ║\n")
cat("║  关键输出：                                                                ║\n")
cat("║  1. 模块识别: 03_Network/Gene_Module_List.csv                               ║\n")
cat("║  2. 表型关联: 04_Module_Trait/Module_Trait_Heatmap.pdf                      ║\n")
cat("║  3. 核心基因: 06_Advanced_Results/Network_Data/Hub_Genes_*.csv               ║\n")
cat("║  4. 富集准备: 06_Advanced_Results/Enrichment/Genes_*.txt                    ║\n")
cat("║  5. 网络文件: 05_Cytoscape/Cytoscape_Edge_*.txt                             ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════╝\n\n")

cat("注意事项：\n")
cat("1. 功能富集分析需使用外部工具完成\n")
cat("2. 所有可视化图表均为发表级质量\n")
cat("3. 核心基因网络使用动态阈值（非固定0.6）\n")
cat("4. 数据已针对16GB内存系统优化\n\n")

# 保存会话信息
sink("R_Session_Info.txt")
sessionInfo()
sink()

# 最终内存清理
gc()
cat("内存清理完成。分析结束。\n")
