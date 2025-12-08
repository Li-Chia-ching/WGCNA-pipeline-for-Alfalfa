# ==============================================================================
# Medicago sativa WGCNA 完整分析流程 v5.0
# 现代化配色方案 + 关键可视化增强 + 第三方分析支持
# ==============================================================================

#### 1. 环境配置与初始化 ####
# ==============================================================================
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════════════╗\n")
cat("║          Medicago sativa WGCNA Analysis v5.0: Modern Visualization         ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════╝\n\n")

# 1.1 设置关键选项
options(stringsAsFactors = FALSE)
options(warn = 1)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 1.2 加载必要包
cat("[Step 1] Loading required packages...\n")
required_packages <- c(
  "WGCNA", "dplyr", "stringr", "ggplot2", "RColorBrewer", "viridis",
  "pheatmap", "igraph", "patchwork", "ggrepel", "gridExtra", "corrplot",
  "ComplexHeatmap", "circlize", "scales", "data.table", "reshape2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# 1.3 配置并行处理
n_cores <- max(1, min(parallel::detectCores() - 2, 8))
enableWGCNAThreads(nThreads = n_cores)
cat(sprintf("  Using %d CPU cores\n", n_cores))

# 1.4 定义现代化配色方案
modern_colors <- list(
  # 热图配色
  heatmap_sequential = c("#2C3E50", "#3498DB", "#1ABC9C", "#F1C40F", "#E74C3C"),
  heatmap_diverging = c("#2980B9", "#6DD5FA", "#FFFFFF", "#FFB347", "#C0392B"),
  viridis_heatmap = viridis(100),
  plasma_heatmap = plasma(100),
  magma_heatmap = magma(100),
  
  # 网络图配色
  network_main = "#2C3E50",
  network_hub = "#E74C3C",
  network_edge = "#3498DB",
  network_highlight = "#F1C40F",
  
  # 模块配色
  module_palette = c("#E74C3C", "#3498DB", "#2ECC71", "#F1C40F", "#9B59B6", 
                     "#1ABC9C", "#34495E", "#E67E22", "#95A5A6", "#D35400")
)

# 1.5 创建输出目录结构
cat("[Step 1] Creating output directories...\n")
output_dirs <- c(
  "01_Data_Clean",
  "02_Network_Topology", 
  "03_Modules",
  "04_Module_Trait_Correlations",
  "05_Hub_Genes",
  "05_Hub_Genes/magenta",
  "05_Hub_Genes/darkgreen",
  "06_Enrichment_Analysis",
  "07_Cytoscape_Files",
  "08_Summary_Reports",
  "09_Third_Party_Analysis"
)

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  Created: %s\n", dir))
  }
}

#### 2. 数据预处理 ####
# ==============================================================================
cat("\n[Step 2] Data preprocessing...\n")

# 2.1 检查输入数据
if (!exists("Summary_of_manually_filtered_genes")) {
  stop("ERROR: Missing transcriptome data 'Summary_of_manually_filtered_genes'")
}
if (!exists("Transcriptome_grouping_summary")) {
  stop("ERROR: Missing grouping data 'Transcriptome_grouping_summary'")
}

# 2.2 处理表达数据
dat_raw <- Summary_of_manually_filtered_genes
fpkm_cols <- grep("(:|\\.)fpkm$", colnames(dat_raw), value = TRUE, ignore.case = TRUE)

if (length(fpkm_cols) == 0) {
  stop("ERROR: No columns with ':fpkm' suffix found")
}

datExpr0 <- as.data.frame(t(dat_raw[, fpkm_cols]))
colnames(datExpr0) <- dat_raw[[1]]
rownames(datExpr0) <- gsub("(:|\\.)fpkm$", "", rownames(datExpr0))

# 2.3 处理表型数据
trait_data <- Transcriptome_grouping_summary
trait_sample_col <- colnames(trait_data)[1]

# 对齐样本
common_samples <- intersect(rownames(datExpr0), trait_data[[trait_sample_col]])
datExpr <- datExpr0[common_samples, ]
trait_data <- trait_data[match(common_samples, trait_data[[trait_sample_col]]), ]

# 2.4 数据转换
datExpr <- log2(datExpr + 1)

# 2.5 构建表型矩阵
datTraits <- data.frame(row.names = rownames(datExpr))

# 处理分组信息
groups <- unique(trait_data$Group_Short_Name)
for (g in groups) {
  datTraits[[paste0("Group_", g)]] <- ifelse(trait_data$Group_Short_Name == g, 1, 0)
}

# 添加品系和光周期特征
datTraits$Stem_Green <- ifelse(grepl("GS", trait_data$Group_Short_Name), 1, 0)
datTraits$Stem_Red <- ifelse(grepl("RS", trait_data$Group_Short_Name), 1, 0)
datTraits$Light_Long <- ifelse(grepl("LL", trait_data$Group_Short_Name), 1, 0)
datTraits$Light_Short <- ifelse(grepl("SL", trait_data$Group_Short_Name), 1, 0)

# 移除零方差列
datTraits <- datTraits[, apply(datTraits, 2, var) > 0]

# 保存预处理数据
save(datExpr, datTraits, file = "01_Data_Clean/Processed_Data.RData")
cat("  Data preprocessing completed\n")

#### 3. 网络构建与模块识别 ####
# ==============================================================================
cat("\n[Step 3] Network construction and module identification...\n")

# 3.1 软阈值选择
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned")

# 自动选择阈值
estimate_power <- sft$powerEstimate
if (is.na(estimate_power)) {
  idx <- which(-sign(sft$fitIndices[,3]) * sft$fitIndices[,2] > 0.80)
  estimate_power <- ifelse(length(idx) > 0, powers[idx[1]], 6)
}
cat(sprintf("  Selected soft threshold: %d\n", estimate_power))

# 3.2 网络拓扑图（现代化配色）
pdf("02_Network_Topology/Network_Topology_Analysis.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# 图1：无尺度拓扑拟合
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, R²",
     type = "n", main = "Scale Independence",
     cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
grid(col = "gray90", lty = 2)
points(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       pch = 21, bg = modern_colors$heatmap_sequential[2], col = "white", 
       cex = 1.5, lwd = 2)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.8, col = modern_colors$network_main, pos = 3)
abline(h = 0.85, col = modern_colors$network_hub, lwd = 2, lty = 2)
abline(v = estimate_power, col = modern_colors$network_highlight, lwd = 2, lty = 2)
legend("bottomright", legend = paste("Selected power:", estimate_power),
       col = modern_colors$network_highlight, lwd = 2, bty = "n")

# 图2：平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity",
     cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
grid(col = "gray90", lty = 2)
points(sft$fitIndices[,1], sft$fitIndices[,5],
       pch = 21, bg = modern_colors$heatmap_sequential[3], col = "white",
       cex = 1.5, lwd = 2)
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = 0.8, col = modern_colors$network_main, pos = 3)
abline(v = estimate_power, col = modern_colors$network_highlight, lwd = 2, lty = 2)

dev.off()

# 3.3 构建共表达网络
net <- blockwiseModules(
  datExpr,
  power = estimate_power,
  maxBlockSize = 15000,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "03_Modules/Medicago_TOM",
  verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 3.4 保存模块结果
save(net, moduleLabels, moduleColors, MEs, 
     file = "03_Modules/Module_Results.RData")

# 3.5 模块统计
module_table <- table(moduleColors)
module_stats <- data.frame(
  Module = names(module_table),
  Size = as.numeric(module_table),
  stringsAsFactors = FALSE
)
module_stats <- module_stats[order(-module_stats$Size), ]

write.csv(module_stats, "03_Modules/Module_Statistics.csv", row.names = FALSE)

# 3.6 模块层级图（现代化配色）
pdf("03_Modules/Module_Dendrogram_Modern.pdf", width = 14, height = 8)
plotDendroAndColors(geneTree, moduleColors,
                    "Module Colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene Clustering and Module Assignment",
                    marAll = c(1, 5, 3, 1))
dev.off()

#### 4. 模块与处理条件相关性热图（重点图1 - 稳健版）####
# ==============================================================================
cat("\n[Step 4] Module-trait correlation analysis (Key Figure 1)...\n")

# 4.1 计算相关性
MEs_ordered <- orderMEs(MEs)
moduleTraitCor <- cor(MEs_ordered, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 4.2 动态准备注释信息（关键修正点）
# --- 自动判断每个性状的类型 ---
trait_names <- colnames(datTraits)
trait_types <- sapply(trait_names, function(name) {
  if (grepl("^Group_", name)) return("Group Effect")
  if (grepl("Stem_", name)) return("Stem Color")
  if (grepl("Light_", name)) return("Photoperiod")
  return("Other") # 兜底分类
})

# --- 创建与热图列数严格对应的顶部注释 ---
ha_column <- HeatmapAnnotation(
  Trait_Type = trait_types,
  col = list(Trait_Type = c("Group Effect" = "#4DAF4A", 
                            "Stem Color" = "#984EA3", 
                            "Photoperiod" = "#FF7F00",
                            "Other" = "#CCCCCC")), # 其他类型用灰色
  annotation_name_side = "left",
  show_legend = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# --- 创建行注释（模块大小）---
module_sizes <- table(moduleColors)
# 安全匹配模块名，防止因模块缺失导致长度错误
module_size_vec <- sapply(rownames(moduleTraitCor), function(me_name) {
  mod_name <- gsub("ME", "", me_name)
  as.numeric(module_sizes[mod_name])
})

ha_row <- rowAnnotation(
  Module_Size = anno_barplot(module_size_vec,
                             gp = gpar(fill = "#377EB8", border = NA),
                             axis_param = list(gp = gpar(fontsize = 9)),
                             width = unit(2, "cm")), # 控制条形图宽度
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# 4.3 创建现代化热图（使用ComplexHeatmap）
pdf("04_Module_Trait_Correlations/Module_Trait_Correlation_Heatmap_Modern.pdf", 
    width = 14, height = 10) # 调整了宽高比，更适合横向性状

library(ComplexHeatmap)
library(circlize)

# 定义颜色映射（蓝-白-红）
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

# 绘制主热图
ht <- Heatmap(moduleTraitCor,
              name = "Correlation\n(r)", # 图例名称更清晰
              col = col_fun,
              
              # --- 矩阵样式 ---
              rect_gp = gpar(col = "white", lwd = 0.5),
              na_col = "grey90",
              
              # --- 聚类设置 ---
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",
              show_column_dend = TRUE,
              show_row_dend = TRUE,
              
              # --- 行列名称设置 ---
              row_names_gp = gpar(fontsize = 11, fontface = "bold"),
              column_names_gp = gpar(fontsize = 11),
              column_names_rot = 45,
              column_names_centered = FALSE,
              
              # --- 标题 ---
              column_title = "Module-Trait Relationships",
              column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "#2C3E50"),
              column_title_side = "top",
              
              # --- 添加注释（关键：确保注释与矩阵对齐）---
              top_annotation = ha_column,
              right_annotation = ha_row,
              
              # --- 图例设置 ---
              heatmap_legend_param = list(
                title = "Pearson r",
                title_gp = gpar(fontsize = 11, fontface = "bold"),
                labels_gp = gpar(fontsize = 10),
                legend_height = unit(3, "cm"),
                grid_width = unit(0.5, "cm"),
                at = c(-1, -0.5, 0, 0.5, 1)
              ),
              
              # --- 单元格内添加显著性星号（优化位置）---
              cell_fun = function(j, i, x, y, width, height, fill) {
                pval <- moduleTraitPvalue[i, j]
                if(!is.na(pval) && pval < 0.05) {
                  star <- ifelse(pval < 0.001, "***", 
                                 ifelse(pval < 0.01, "**", "*"))
                  # 将星号放在单元格上半部分，避免与数字重叠（如果以后加数字）
                  grid.text(star, x, y - height*0.25, 
                            gp = gpar(fontsize = 9, fontface = "bold", col = "black"))
                }
              },
              
              # --- 其他优化 ---
              show_heatmap_legend = TRUE,
              use_raster = FALSE # 数据量大时可设为TRUE加速渲染
)

# 4.4 绘制并保存
draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     padding = unit(c(4, 2, 4, 2), "mm") # 上右下左边距
)

dev.off()

# 4.5 同时生成一个PNG版本（用于快速查看或演示）
png("04_Module_Trait_Correlations/Module_Trait_Correlation_Heatmap_Modern.png", 
    width = 1400, height = 1000, res = 150)
draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     padding = unit(c(4, 2, 4, 2), "mm"))
dev.off()

cat("  ✓ 现代化相关性热图已保存 (PDF & PNG)。\n")
cat("  ✓ 顶部注释基于性状名称自动生成：", paste(unique(trait_types), collapse=", "), "\n")

# 4.6 保存详细的、带标注的相关性结果表格
cor_results <- data.frame(
  Module = rownames(moduleTraitCor),
  Module_Color = gsub("ME", "", rownames(moduleTraitCor)),
  Module_Size = module_size_vec,
  stringsAsFactors = FALSE
)

for (i in seq_along(trait_names)) {
  trait <- trait_names[i]
  type <- trait_types[i]
  
  cor_results[[paste0(trait, "_r")]] <- round(moduleTraitCor[, trait], 4)
  cor_results[[paste0(trait, "_pval")]] <- signif(moduleTraitPvalue[, trait], 4)
  cor_results[[paste0(trait, "_FDR")]] <- signif(p.adjust(moduleTraitPvalue[, trait], method = "fdr"), 4)
  
  # 添加一列直接显示显著性与相关系数，便于阅读
  cor_results[[paste0(trait, "_summary")]] <- sapply(1:nrow(moduleTraitCor), function(idx) {
    r_val <- moduleTraitCor[idx, trait]
    p_val <- moduleTraitPvalue[idx, trait]
    if (is.na(p_val) || p_val >= 0.05) return(sprintf("%.3f", r_val))
    star <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))
    return(sprintf("%.3f%s", r_val, star))
  })
}

# 保存两份：一份完整，一份精简（仅含显著结果）
write.csv(cor_results, 
          "04_Module_Trait_Correlations/Detailed_Correlation_Results_Full.csv", 
          row.names = FALSE)

# 精简版：找出至少与一个性状显著相关(p<0.05)的模块
significant_modules <- unique(unlist(lapply(trait_names, function(trait) {
  rownames(moduleTraitCor)[moduleTraitPvalue[, trait] < 0.05]
})))
if (length(significant_modules) > 0) {
  sig_cor_results <- cor_results[cor_results$Module %in% significant_modules, ]
  write.csv(sig_cor_results,
            "04_Module_Trait_Correlations/Significant_Correlation_Results.csv",
            row.names = FALSE)
  cat("  ✓ 已保存显著相关结果。共", length(significant_modules), "个模块与至少一个性状显著相关。\n")
} else {
  cat("  ⚠️  未发现显著性关联 (p < 0.05)。\n")
}

#### 5. 核心基因分析（重点模块） ####
# ==============================================================================
cat("\n[Step 5] Hub gene analysis for key modules...\n")

# 5.1 核心基因分析函数（增强版）
analyze_hub_genes_enhanced <- function(target_module, n_top = 25, save_dir = "05_Hub_Genes") {
  cat(sprintf("  Analyzing hub genes for module: %s\n", target_module))
  
  # 提取模块基因
  module_genes <- colnames(datExpr)[moduleColors == target_module]
  
  if (length(module_genes) < 10) {
    cat(sprintf("    Warning: Module %s has only %d genes, skipping\n", 
                target_module, length(module_genes)))
    return(NULL)
  }
  
  # 计算kME
  me_name <- paste0("ME", target_module)
  if (!me_name %in% colnames(MEs)) {
    cat(sprintf("    Error: Module eigengene not found for %s\n", target_module))
    return(NULL)
  }
  
  kME_values <- signedKME(datExpr[, module_genes], MEs[, me_name, drop = FALSE])
  
  kME_df <- data.frame(
    Gene = module_genes,
    kME = kME_values[, 1],
    stringsAsFactors = FALSE
  )
  
  kME_df <- kME_df[order(-abs(kME_df$kME)), ]
  hub_genes <- head(kME_df, n_top)
  
  return(list(
    all_genes = kME_df,
    hub_genes = hub_genes,
    module_name = target_module
  ))
}

# 5.2 分析magenta模块（重点模块）
cat("  Processing magenta module...\n")
magenta_results <- analyze_hub_genes_enhanced("magenta", n_top = 25, save_dir = "05_Hub_Genes/magenta")

if (!is.null(magenta_results)) {
  # 保存基因列表
  write.csv(magenta_results$hub_genes, 
            "05_Hub_Genes/magenta/Magenta_Hub_Genes.csv", 
            row.names = FALSE)
  write.csv(magenta_results$all_genes, 
            "05_Hub_Genes/magenta/Magenta_All_Genes_with_kME.csv", 
            row.names = FALSE)
  
  # 5.2.1 Magenta核心基因表达热图（重点图2）
  cat("    Creating magenta hub gene expression heatmap...\n")
  
  # 选择前20个核心基因
  top_genes <- head(magenta_results$hub_genes$Gene, 20)
  expr_data <- datExpr[, top_genes, drop = FALSE]
  
  # Z-score标准化
  expr_scaled <- t(scale(expr_data))
  
  # 创建现代化热图
  pdf("05_Hub_Genes/magenta/Magenta_Hub_Gene_Heatmap_Modern.pdf", 
      width = 14, height = 10)
  
  # 创建样本注释（按处理条件）
  annotation_col <- data.frame(
    Stem_Color = ifelse(grepl("GS", trait_data$Group_Short_Name), "Green", "Red"),
    Light_Condition = ifelse(grepl("LL", trait_data$Group_Short_Name), "Long", "Short"),
    row.names = rownames(datExpr)
  )
  
  # 定义注释颜色
  ann_colors <- list(
    Stem_Color = c(Green = "#2ECC71", Red = "#E74C3C"),
    Light_Condition = c(Long = "#3498DB", Short = "#F1C40F")
  )
  
  # 使用pheatmap创建现代化热图
  pheatmap(expr_scaled,
           main = "Magenta Module: Hub Gene Expression Patterns",
           color = modern_colors$viridis_heatmap,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = ncol(expr_scaled) <= 30,
           show_rownames = TRUE,
           fontsize_row = 9,
           fontsize_col = 8,
           fontsize = 11,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           border_color = NA,
           cellwidth = ifelse(ncol(expr_scaled) <= 30, 12, NA),
           cellheight = 12)
  
  dev.off()
  
  # 5.2.2 Magenta核心基因互作网络图（重点图3）
  cat("    Creating magenta hub gene interaction network...\n")
  
  # 创建互作网络
  create_gene_network <- function(gene_list, expr_data, module_color, output_dir) {
    if (length(gene_list) < 5) return(NULL)
    
    # 计算相关性矩阵
    cor_matrix <- cor(expr_data[, gene_list, drop = FALSE], use = "pairwise.complete.obs")
    
    # 转换为邻接矩阵（动态阈值）
    threshold <- quantile(abs(cor_matrix[upper.tri(cor_matrix)]), 0.85, na.rm = TRUE)
    adj_matrix <- abs(cor_matrix)
    adj_matrix[adj_matrix < threshold] <- 0
    diag(adj_matrix) <- 0
    
    # 创建igraph对象
    g <- graph_from_adjacency_matrix(adj_matrix, 
                                     mode = "undirected", 
                                     weighted = TRUE,
                                     diag = FALSE)
    
    # 移除孤立节点
    g <- delete_vertices(g, degree(g) == 0)
    
    if (vcount(g) < 3) return(NULL)
    
    # 计算节点属性
    node_degree <- degree(g)
    betweenness_centrality <- betweenness(g)
    
    # 创建现代化网络图
    pdf(sprintf("%s/%s_Hub_Gene_Network_Modern.pdf", output_dir, module_color),
        width = 12, height = 10)
    
    # 设置布局
    layout <- layout_with_fr(g)
    
    # 设置节点颜色（按度大小）
    node_colors <- colorRampPalette(c("#ECF0F1", modern_colors$network_main))(100)
    degree_scaled <- as.integer(cut(node_degree, breaks = 100))
    
    # 设置节点大小（按中心性）
    node_size <- 10 + (betweenness_centrality / max(betweenness_centrality)) * 20
    
    # 设置边宽度和颜色
    edge_width <- E(g)$weight * 3
    edge_colors <- colorRampPalette(c("#BDC3C7", modern_colors$network_edge))(100)
    edge_weight_scaled <- as.integer(cut(E(g)$weight, breaks = 100))
    
    # 绘制网络
    plot(g,
         layout = layout,
         vertex.size = node_size,
         vertex.color = node_colors[degree_scaled],
         vertex.frame.color = "white",
         vertex.label.cex = 0.8,
         vertex.label.color = modern_colors$network_main,
         vertex.label.family = "sans",
         edge.width = edge_width,
         edge.color = edge_colors[edge_weight_scaled],
         edge.curved = 0.2,
         main = sprintf("%s Module: Hub Gene Interaction Network", toupper(module_color)),
         main.cex = 1.2)
    
    # 添加图例
    legend("bottomright",
           legend = c("High Connectivity", "Low Connectivity"),
           fill = c(node_colors[100], node_colors[1]),
           title = "Connectivity",
           bty = "n",
           cex = 0.8)
    
    legend("bottomleft",
           legend = c("Strong Interaction", "Weak Interaction"),
           col = c(edge_colors[100], edge_colors[1]),
           lwd = 3,
           title = "Edge Strength",
           bty = "n",
           cex = 0.8)
    
    dev.off()
    
    # 保存网络数据
    write.csv(as_data_frame(g, what = "edges"),
              sprintf("%s/%s_Network_Edges.csv", output_dir, module_color),
              row.names = FALSE)
    
    write.csv(data.frame(Gene = V(g)$name,
                         Degree = node_degree,
                         Betweenness = betweenness_centrality),
              sprintf("%s/%s_Network_Node_Attributes.csv", output_dir, module_color),
              row.names = FALSE)
    
    return(g)
  }
  
  # 创建magenta网络图
  magenta_network <- create_gene_network(
    top_genes,
    datExpr,
    "magenta",
    "05_Hub_Genes/magenta"
  )
}

# 5.3 分析darkgreen模块（重点模块）
cat("  Processing darkgreen module...\n")
darkgreen_results <- analyze_hub_genes_enhanced("darkgreen", n_top = 25, save_dir = "05_Hub_Genes/darkgreen")

if (!is.null(darkgreen_results)) {
  # 保存基因列表
  write.csv(darkgreen_results$hub_genes, 
            "05_Hub_Genes/darkgreen/Darkgreen_Hub_Genes.csv", 
            row.names = FALSE)
  write.csv(darkgreen_results$all_genes, 
            "05_Hub_Genes/darkgreen/Darkgreen_All_Genes_with_kME.csv", 
            row.names = FALSE)
  
  # 5.3.1 Darkgreen核心基因表达热图（重点图4）
  cat("    Creating darkgreen hub gene expression heatmap...\n")
  
  top_genes <- head(darkgreen_results$hub_genes$Gene, 20)
  expr_data <- datExpr[, top_genes, drop = FALSE]
  expr_scaled <- t(scale(expr_data))
  
  pdf("05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Heatmap_Modern.pdf", 
      width = 14, height = 10)
  
  # 样本注释
  annotation_col <- data.frame(
    Stem_Color = ifelse(grepl("GS", trait_data$Group_Short_Name), "Green", "Red"),
    Light_Condition = ifelse(grepl("LL", trait_data$Group_Short_Name), "Long", "Short"),
    row.names = rownames(datExpr)
  )
  
  ann_colors <- list(
    Stem_Color = c(Green = "#2ECC71", Red = "#E74C3C"),
    Light_Condition = c(Long = "#3498DB", Short = "#F1C40F")
  )
  
  pheatmap(expr_scaled,
           main = "Darkgreen Module: Hub Gene Expression Patterns",
           color = modern_colors$plasma_heatmap,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = ncol(expr_scaled) <= 30,
           show_rownames = TRUE,
           fontsize_row = 9,
           fontsize_col = 8,
           fontsize = 11,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           border_color = NA,
           cellwidth = ifelse(ncol(expr_scaled) <= 30, 12, NA),
           cellheight = 12)
  
  dev.off()
  
  # 5.3.2 Darkgreen核心基因互作网络图（重点图5）
  cat("    Creating darkgreen hub gene interaction network...\n")
  
  darkgreen_network <- create_gene_network(
    top_genes,
    datExpr,
    "darkgreen",
    "05_Hub_Genes/darkgreen"
  )
}

#### 6. 富集分析准备（KEGG第三方分析） ####
# ==============================================================================
cat("\n[Step 6] Preparing files for external KEGG enrichment analysis...\n")

# 6.1 识别关键模块（基于相关性）
identify_key_modules <- function(cor_results, trait_patterns, top_n = 3) {
  # 筛选与特定性状高度相关的模块
  key_modules <- list()
  
  for (pattern in trait_patterns) {
    # 找到与当前性状最相关的模块
    trait_cols <- grep(pattern, colnames(cor_results), value = TRUE)
    cor_cols <- trait_cols[grep("_Cor$", trait_cols)]
    
    if (length(cor_cols) > 0) {
      max_cor <- 0
      best_module <- NULL
      
      for (cor_col in cor_cols) {
        max_idx <- which.max(abs(cor_results[[cor_col]]))
        current_max <- abs(cor_results[max_idx, cor_col])
        
        if (current_max > max_cor) {
          max_cor <- current_max
          best_module <- cor_results$Module[max_idx]
        }
      }
      
      if (!is.null(best_module)) {
        module_name <- gsub("ME", "", best_module)
        key_modules[[pattern]] <- module_name
      }
    }
  }
  
  return(key_modules)
}

# 识别关键模块
trait_patterns <- c("Stem_Green", "Stem_Red", "Light_Long", "Light_Short")
key_modules <- identify_key_modules(cor_results, trait_patterns)

cat("  Key modules identified for enrichment analysis:\n")
for (trait in names(key_modules)) {
  cat(sprintf("    %s: %s module\n", trait, key_modules[[trait]]))
}

# 6.2 为关键模块准备KEGG分析文件
cat("  Preparing gene lists for external KEGG analysis...\n")

# 创建富集分析指南
enrichment_guide <- c(
  "=================================================================================",
  "EXTERNAL KEGG ENRICHMENT ANALYSIS GUIDE",
  paste("Generated:", Sys.time()),
  "=================================================================================",
  "",
  "IMPORTANT: Medicago sativa (Alfalfa) has limited direct support in R KEGG packages.",
  "This analysis prepares gene lists for external enrichment tools.",
  "",
  "KEY MODULES IDENTIFIED FOR ENRICHMENT ANALYSIS:",
  ""
)

# 添加模块信息
for (trait in names(key_modules)) {
  module_name <- key_modules[[trait]]
  enrichment_guide <- c(enrichment_guide,
                        sprintf("- %s: Strongly associated with %s", module_name, trait))
}

enrichment_guide <- c(enrichment_guide,
                      "",
                      "RECOMMENDED EXTERNAL TOOLS:",
                      "",
                      "1. KEGG Mapper (https://www.kegg.jp/kegg/mapper/)",
                      "   - Upload gene lists in '09_Third_Party_Analysis/KEGG_Input/'",
                      "   - Use 'Medicago truncatula' as reference species (KEGG code: mtr)",
                      "   - Select 'Search Against: KEGG Genes Database'",
                      "",
                      "2. KOBAS (http://kobas.cbi.pku.edu.cn/)",
                      "   - Supports Medicago truncatula annotations",
                      "   - Provides KEGG and GO enrichment",
                      "",
                      "3. ShinyGO (http://bioinformatics.sdstate.edu/go/)",
                      "   - User-friendly web interface",
                      "   - Supports multiple plant species",
                      "",
                      "ANALYSIS STEPS:",
                      "1. Upload 'Genes_[module].txt' files to your chosen tool",
                      "2. Select Medicago truncatula as reference (for KEGG)",
                      "3. Run enrichment analysis with default parameters",
                      "4. Download results and import into R for visualization",
                      "",
                      "VISUALIZATION IN R (after external analysis):",
                      "1. Load your enrichment results into R",
                      "2. Use the template in '09_Third_Party_Analysis/Visualization_Template.R'",
                      "3. Create publication-quality enrichment plots",
                      "",
                      "FILE STRUCTURE:",
                      "09_Third_Party_Analysis/",
                      "├── KEGG_Input/          # Gene lists for external tools",
                      "├── Expected_Output/     # Expected result formats",
                      "└── Visualization_Template.R  # R code for result visualization",
                      "",
                      "================================================================================="
)

# 保存指南
writeLines(enrichment_guide, "06_Enrichment_Analysis/KEGG_Enrichment_Guide.txt")

# 6.3 创建KEGG分析输入文件
kegg_input_dir <- "09_Third_Party_Analysis/KEGG_Input"
dir.create(kegg_input_dir, recursive = TRUE, showWarnings = FALSE)

# 为所有模块准备基因列表
all_modules <- unique(moduleColors[moduleColors != "grey"])
for (module in all_modules) {
  module_genes <- colnames(datExpr)[moduleColors == module]
  
  # 保存基因列表
  write.table(module_genes,
              file.path(kegg_input_dir, sprintf("Genes_%s.txt", module)),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 保存带kME值的列表
  if (module %in% c("magenta", "darkgreen")) {
    # 为关键模块保存详细信息
    if (module == "magenta" && !is.null(magenta_results)) {
      write.csv(magenta_results$all_genes,
                file.path(kegg_input_dir, sprintf("Genes_%s_with_kME.csv", module)),
                row.names = FALSE)
    }
    if (module == "darkgreen" && !is.null(darkgreen_results)) {
      write.csv(darkgreen_results$all_genes,
                file.path(kegg_input_dir, sprintf("Genes_%s_with_kME.csv", module)),
                row.names = FALSE)
    }
  }
}

# 6.4 创建可视化模板
visualization_template <- c(
  '# ==============================================================================',
  '# R VISUALIZATION TEMPLATE FOR EXTERNAL ENRICHMENT RESULTS',
  '# Use this after running KEGG/GO enrichment on external platforms',
  '# ==============================================================================',
  '',
  'library(ggplot2)',
  'library(dplyr)',
  'library(tidyr)',
  'library(RColorBrewer)',
  '',
  '# 1. LOAD YOUR ENRICHMENT RESULTS',
  '# Example: Load KEGG enrichment results (modify file path as needed)',
  'enrichment_data <- read.csv("path/to/your/kegg_results.csv", stringsAsFactors = FALSE)',
  '',
  '# 2. DATA PREPARATION',
  '# Filter significant results',
  'significant_results <- enrichment_data %>%',
  '  filter(p.adjust < 0.05) %>%',
  '  arrange(p.adjust) %>%',
  '  head(15)',
  '',
  '# 3. CREATE MODERN BAR PLOT',
  'p_bar <- ggplot(significant_results, aes(x = reorder(Description, -log10(p.adjust)), ',
  '                                          y = -log10(p.adjust))) +',
  '  geom_bar(stat = "identity", fill = "#3498DB", alpha = 0.8, width = 0.7) +',
  '  coord_flip() +',
  '  labs(title = "Top Enriched KEGG Pathways",',
  '       x = "Pathway",',
  '       y = "-log10(Adjusted P-value)") +',
  '  theme_minimal(base_size = 12) +',
  '  theme(',
  '    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),',
  '    axis.title = element_text(face = "bold"),',
  '    panel.grid.major = element_line(color = "gray90"),',
  '    panel.grid.minor = element_blank()',
  '  )',
  '',
  '# 4. CREATE DOT PLOT WITH GENE RATIO',
  'if ("GeneRatio" %in% colnames(significant_results)) {',
  '  p_dot <- ggplot(significant_results, aes(x = GeneRatio, ',
  '                                            y = reorder(Description, GeneRatio))) +',
  '    geom_point(aes(size = Count, color = -log10(p.adjust))) +',
  '    scale_color_gradient(low = "#3498DB", high = "#E74C3C") +',
  '    scale_size(range = c(3, 8)) +',
  '    labs(title = "Enrichment Dot Plot",',
  '         x = "Gene Ratio",',
  '         y = "Pathway",',
  '         size = "Gene Count",',
  '         color = "-log10(P-value)") +',
  '    theme_minimal() +',
  '    theme(plot.title = element_text(face = "bold", hjust = 0.5))',
  '}',
  '',
  '# 5. CREATE PATHWAY CATEGORY PLOT',
  '# If you have pathway categories in your data',
  'if ("Category" %in% colnames(significant_results)) {',
  '  p_category <- ggplot(significant_results, aes(x = Category, y = -log10(p.adjust))) +',
  '    geom_boxplot(fill = "#2ECC71", alpha = 0.6) +',
  '    geom_jitter(width = 0.2, size = 2, color = "#2C3E50") +',
  '    labs(title = "Enrichment by Pathway Category",',
  '         x = "Category",',
  '         y = "-log10(Adjusted P-value)") +',
  '    theme_minimal() +',
  '    theme(axis.text.x = element_text(angle = 45, hjust = 1))',
  '}',
  '',
  '# 6. SAVE PLOTS',
  'ggsave("09_Third_Party_Analysis/KEGG_Enrichment_BarPlot.pdf", p_bar, width = 12, height = 8)',
  'if (exists("p_dot")) {',
  '  ggsave("09_Third_Party_Analysis/KEGG_Enrichment_DotPlot.pdf", p_dot, width = 12, height = 8)',
  '}',
  '',
  '# 7. CREATE SUMMARY TABLE',
  'write.csv(significant_results,',
  '          "09_Third_Party_Analysis/Significant_Enrichment_Results.csv",',
  '          row.names = FALSE)',
  '',
  'cat("Enrichment visualization completed. Check 09_Third_Party_Analysis/ directory.\\n")',
  ''
)

writeLines(visualization_template, "09_Third_Party_Analysis/Visualization_Template.R")

#### 7. Cytoscape文件导出 ####
# ==============================================================================
cat("\n[Step 7] Exporting files for Cytoscape visualization...\n")

export_cytoscape_files <- function(target_module, max_genes = 200) {
  cat(sprintf("  Exporting %s module for Cytoscape...\n", target_module))
  
  # 获取模块基因
  module_genes <- colnames(datExpr)[moduleColors == target_module]
  
  if (length(module_genes) > max_genes) {
    # 如果基因太多，选择前max_genes个核心基因
    me_name <- paste0("ME", target_module)
    if (me_name %in% colnames(MEs)) {
      kME_values <- abs(cor(datExpr[, module_genes], MEs[, me_name], use = "p"))
      ranked_genes <- rownames(kME_values)[order(-kME_values, decreasing = TRUE)]
      module_genes <- ranked_genes[1:max_genes]
    }
  }
  
  if (length(module_genes) >= 10) {
    # 计算TOM
    modExpr <- datExpr[, module_genes, drop = FALSE]
    TOM <- TOMsimilarityFromExpr(modExpr, power = estimate_power, networkType = "unsigned")
    
    # 设定阈值
    threshold <- quantile(TOM[lower.tri(TOM)], probs = 0.85, na.rm = TRUE)
    
    # 导出Cytoscape文件
    exportNetworkToCytoscape(
      TOM,
      edgeFile = sprintf("07_Cytoscape_Files/%s_EdgeList.txt", target_module),
      nodeFile = sprintf("07_Cytoscape_Files/%s_NodeList.txt", target_module),
      weighted = TRUE,
      threshold = threshold,
      nodeNames = module_genes,
      nodeAttr = rep(target_module, length(module_genes))
    )
    
    return(TRUE)
  }
  return(FALSE)
}

# 导出关键模块
modules_to_export <- c("magenta", "darkgreen", 
                       module_stats$Module[1:min(5, nrow(module_stats))])

for (module in unique(modules_to_export)) {
  if (module != "grey") {
    success <- export_cytoscape_files(module)
    if (success) {
      cat(sprintf("    Exported: %s module\n", module))
    }
  }
}

#### 8. 生成综合总结报告 ####
# ==============================================================================
cat("\n[Step 8] Generating comprehensive summary report...\n")

# 8.1 创建详细总结
summary_report <- c(
  "=================================================================================",
  "MEDICAGO SATIVA WGCNA ANALYSIS - COMPREHENSIVE SUMMARY REPORT",
  paste("Analysis completed:", Sys.time()),
  "=================================================================================",
  "",
  "1. ANALYSIS OVERVIEW",
  paste("   Total genes analyzed:", ncol(datExpr)),
  paste("   Total samples:", nrow(datExpr)),
  paste("   Experimental design: 2×2 factorial (Stem Color × Light Condition)"),
  paste("   Network power:", estimate_power),
  "",
  "2. MODULE DISCOVERY",
  paste("   Total modules detected:", length(unique(moduleColors))),
  paste("   Non-grey modules:", sum(moduleColors != "grey")),
  "",
  "3. KEY VISUALIZATIONS GENERATED",
  "",
  "   [FIGURE 1] Module-Trait Correlation Heatmap",
  "   Location: 04_Module_Trait_Correlations/Module_Trait_Correlation_Heatmap_Modern.pdf",
  "   Description: Shows correlation between co-expression modules and experimental",
  "                conditions (stem color and light period).",
  "   Key findings: Identifies modules strongly associated with specific treatments.",
  "",
  "   [FIGURE 2] Magenta Module Hub Gene Expression Heatmap",
  "   Location: 05_Hub_Genes/magenta/Magenta_Hub_Gene_Heatmap_Modern.pdf",
  "   Description: Expression patterns of top 20 hub genes in the magenta module.",
  "   Key findings: Shows how key regulatory genes respond to different conditions.",
  "",
  "   [FIGURE 3] Magenta Module Hub Gene Interaction Network",
  "   Location: 05_Hub_Genes/magenta/Magenta_Hub_Gene_Network_Modern.pdf",
  "   Description: Network visualization of interactions between hub genes.",
  "   Key findings: Reveals core regulatory network structure.",
  "",
  "   [FIGURE 4] Darkgreen Module Hub Gene Expression Heatmap",
  "   Location: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Heatmap_Modern.pdf",
  "   Description: Expression patterns of top hub genes in darkgreen module.",
  "",
  "   [FIGURE 5] Darkgreen Module Hub Gene Interaction Network",
  "   Location: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Network_Modern.pdf",
  "   Description: Network structure of darkgreen module hub genes.",
  "",
  "4. HUB GENE ANALYSIS",
  "   Key modules with biological relevance:",
  ""
)

# 添加模块特异性信息
if (!is.null(magenta_results)) {
  summary_report <- c(summary_report,
                      paste("   - Magenta module: ", nrow(magenta_results$hub_genes), 
                            " hub genes identified"),
                      paste("     Top 5 hub genes: ", 
                            paste(head(magenta_results$hub_genes$Gene, 5), collapse = ", ")))
}

if (!is.null(darkgreen_results)) {
  summary_report <- c(summary_report,
                      paste("   - Darkgreen module: ", nrow(darkgreen_results$hub_genes), 
                            " hub genes identified"),
                      paste("     Top 5 hub genes: ", 
                            paste(head(darkgreen_results$hub_genes$Gene, 5), collapse = ", ")))
}

summary_report <- c(summary_report,
                    "",
                    "5. FUNCTIONAL ENRICHMENT ANALYSIS PREPARATION",
                    "   IMPORTANT: Medicago sativa requires external tools for KEGG enrichment.",
                    "   Prepared files are available in: 09_Third_Party_Analysis/KEGG_Input/",
                    "",
                    "   Recommended external tools:",
                    "   1. KEGG Mapper (use Medicago truncatula as reference)",
                    "   2. KOBAS (supports M. truncatula annotations)",
                    "   3. ShinyGO (user-friendly web interface)",
                    "",
                    "   Analysis steps:",
                    "   1. Upload gene lists from 09_Third_Party_Analysis/KEGG_Input/",
                    "   2. Select appropriate reference species (M. truncatula for KEGG)",
                    "   3. Run enrichment analysis",
                    "   4. Use Visualization_Template.R for result visualization",
                    "",
                    "6. OUTPUT DIRECTORY STRUCTURE",
                    "",
                    "   01_Data_Clean/                    # Preprocessed data",
                    "   02_Network_Topology/              # Network topology analysis",
                    "   03_Modules/                       # Module identification results",
                    "   04_Module_Trait_Correlations/     # Module-trait relationships (Key Figure 1)",
                    "   05_Hub_Genes/                     # Hub gene analysis",
                    "       ├── magenta/                  # Magenta module analysis (Figures 2-3)",
                    "       └── darkgreen/                # Darkgreen module analysis (Figures 4-5)",
                    "   06_Enrichment_Analysis/           # Enrichment analysis guidance",
                    "   07_Cytoscape_Files/               # Network files for Cytoscape",
                    "   08_Summary_Reports/               # This report and additional summaries",
                    "   09_Third_Party_Analysis/          # Files for external KEGG analysis",
                    "",
                    "7. NEXT STEPS RECOMMENDED",
                    "",
                    "   1. Examine module-trait correlations to identify key modules",
                    "   2. Validate hub genes with qPCR or other experimental methods",
                    "   3. Perform KEGG enrichment using external tools",
                    "   4. Use Cytoscape for network visualization and exploration",
                    "   5. Integrate findings with existing literature on Medicago sativa",
                    "",
                    "8. TECHNICAL NOTES",
                    "",
                    "   - All visualizations use modern color schemes optimized for publication",
                    "   - Network files are compatible with Cytoscape 3.8+",
                    "   - Gene lists are formatted for major enrichment analysis platforms",
                    "   - Analysis is reproducible with provided R session information",
                    "",
                    "=================================================================================",
                    "ANALYSIS COMPLETED SUCCESSFULLY",
                    "================================================================================="
)

# 保存总结报告
writeLines(summary_report, "08_Summary_Reports/Comprehensive_Analysis_Summary.txt")

# 8.2 创建简短的执行总结
exec_summary <- c(
  "╔════════════════════════════════════════════════════════════════════════════╗",
  "║                     WGCNA ANALYSIS EXECUTIVE SUMMARY                       ║",
  "╚════════════════════════════════════════════════════════════════════════════╝",
  "",
  paste("Analysis Date:   ", Sys.time()),
  paste("Genes Analyzed:  ", format(ncol(datExpr), big.mark = ",")),
  paste("Samples:         ", nrow(datExpr)),
  paste("Key Modules:     ", ifelse(!is.null(magenta_results), "magenta", ""), 
        ifelse(!is.null(darkgreen_results), "darkgreen", "")),
  "",
  "KEY OUTPUTS:",
  "✓ Modern module-trait correlation heatmap",
  "✓ Magenta module hub gene analysis (heatmap + network)",
  "✓ Darkgreen module hub gene analysis (heatmap + network)",
  "✓ External KEGG enrichment analysis files",
  "✓ Cytoscape network visualization files",
  "",
  "NEXT ACTIONS REQUIRED:",
  "1. Review module-trait correlations (04_Module_Trait_Correlations/)",
  "2. Perform KEGG enrichment using external tools (09_Third_Party_Analysis/)",
  "3. Validate key hub genes experimentally",
  "",
  "For detailed results, see the comprehensive summary report."
)

writeLines(exec_summary, "08_Summary_Reports/Executive_Summary.txt")

#### 9. 生成综合总结报告 ####
# ==============================================================================
cat("\n[Step 9] Generating comprehensive summary report...\n")

# 9.1 创建详细总结 (完全重构，避免复杂拼接)
# 使用安全的字符串构建方式
divider_line <- paste(rep("=", 80), collapse = "")

summary_report <- c(
  divider_line,
  "MEDICAGO SATIVA WGCNA ANALYSIS - COMPREHENSIVE SUMMARY REPORT",
  paste("Analysis completed:", Sys.time()),
  divider_line,
  "",
  "1. ANALYSIS OVERVIEW",
  paste("   Total genes analyzed:", ncol(datExpr)),
  paste("   Total samples:", nrow(datExpr)),
  paste("   Experimental design: 2x2 factorial (Stem Color x Light Condition)"),
  paste("   Network power:", estimate_power),
  "",
  "2. MODULE DISCOVERY",
  paste("   Total modules detected:", length(unique(moduleColors))),
  paste("   Non-grey modules:", sum(moduleColors != "grey")),  # 移除转义符，更清晰
  "",
  "3. KEY VISUALIZATIONS GENERATED",
  "",
  "   [FIGURE 1] Module-Trait Correlation Heatmap",
  "   Location: 04_Module_Trait_Correlations/Module_Trait_Correlation_Heatmap_Modern.pdf",
  "   Description: Shows correlation between co-expression modules and experimental",
  "                conditions (stem color and light period).",
  "",
  "   [FIGURE 2] Magenta Module Hub Gene Expression Heatmap",
  "   Location: 05_Hub_Genes/magenta/Magenta_Hub_Gene_Heatmap_Modern.pdf",
  "",
  "   [FIGURE 3] Magenta Module Hub Gene Interaction Network",
  "   Location: 05_Hub_Genes/magenta/Magenta_Hub_Gene_Network_Modern.pdf",
  "",
  "   [FIGURE 4] Darkgreen Module Hub Gene Expression Heatmap",
  "   Location: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Heatmap_Modern.pdf",
  "",
  "   [FIGURE 5] Darkgreen Module Hub Gene Interaction Network",
  "   Location: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Network_Modern.pdf",
  "",
  "4. HUB GENE ANALYSIS",
  "   Key modules with biological relevance:",
  ""
)

# 动态添加模块信息 - 使用独立的逻辑块，确保括号闭合
if (exists("magenta_results") && !is.null(magenta_results)) {
  magenta_info <- c(
    paste("   - Magenta module:", nrow(magenta_results$hub_genes), "hub genes identified"),
    paste("     Top 5 hub genes:", paste(head(magenta_results$hub_genes$Gene, 5), collapse = ", "))
  )
  summary_report <- c(summary_report, magenta_info, "")
}

if (exists("darkgreen_results") && !is.null(darkgreen_results)) {
  darkgreen_info <- c(
    paste("   - Darkgreen module:", nrow(darkgreen_results$hub_genes), "hub genes identified"),
    paste("     Top 5 hub genes:", paste(head(darkgreen_results$hub_genes$Gene, 5), collapse = ", "))
  )
  summary_report <- c(summary_report, darkgreen_info, "")
}

# 继续添加报告的其余部分 (结构清晰，逐段添加)
report_footer <- c(
  "",
  "5. FUNCTIONAL ENRICHMENT ANALYSIS PREPARATION",
  "   IMPORTANT: Medicago sativa requires external tools for KEGG enrichment.",
  "   Prepared files are available in: 09_Third_Party_Analysis/KEGG_Input/",
  "",
  "6. OUTPUT DIRECTORY STRUCTURE",
  "   01_Data_Clean/                    # Preprocessed data",
  "   02_Network_Topology/              # Network topology analysis",
  "   03_Modules/                       # Module identification results",
  "   04_Module_Trait_Correlations/     # Module-trait relationships",
  "   05_Hub_Genes/                     # Hub gene analysis",
  "   06_Enrichment_Analysis/           # Enrichment analysis guidance",
  "   07_Cytoscape_Files/               # Network files for Cytoscape",
  "   08_Summary_Reports/               # This report and additional summaries",
  "   09_Third_Party_Analysis/          # Files for external KEGG analysis",
  "",
  "7. NEXT STEPS RECOMMENDED",
  "   1. Examine module-trait correlations to identify key modules",
  "   2. Validate hub genes with qPCR or other experimental methods",
  "   3. Perform KEGG enrichment using external tools",
  "   4. Use Cytoscape for network visualization and exploration",
  "",
  divider_line,
  "ANALYSIS COMPLETED SUCCESSFULLY",
  divider_line
)

summary_report <- c(summary_report, report_footer)

# 保存总结报告
writeLines(summary_report, "08_Summary_Reports/Comprehensive_Analysis_Summary.txt")
cat("  [完成] 综合总结报告已生成\n")

# 9.2 创建简短的执行总结
exec_divider <- paste(rep("═", 70), collapse = "")

exec_summary <- c(
  exec_divider,
  "               WGCNA ANALYSIS EXECUTIVE SUMMARY               ",
  exec_divider,
  paste("Analysis Date:  ", Sys.time()),
  paste("Genes Analyzed: ", format(ncol(datExpr), big.mark = ",")),
  paste("Samples:        ", nrow(datExpr)),
  paste("Network Power:  ", estimate_power),
  paste("Key Modules:    ", ifelse(exists("magenta_results"), "magenta ", ""),
        ifelse(exists("darkgreen_results"), "darkgreen", "")),
  "",
  "KEY OUTPUTS:",
  "• Module-trait correlation heatmap",
  "• Magenta module hub gene analysis (heatmap + network)",
  "• Darkgreen module hub gene analysis (heatmap + network)",
  "• External KEGG enrichment analysis files",
  "• Cytoscape network visualization files",
  "",
  "NEXT ACTIONS:",
  "1. Review module-trait correlations in 04_Module_Trait_Correlations/",
  "2. Use external tools for KEGG enrichment (see 06_Enrichment_Analysis/)",
  "3. Validate key hub genes experimentally",
  "",
  exec_divider
)

writeLines(exec_summary, "08_Summary_Reports/Executive_Summary.txt")
cat("  [完成] 执行摘要已生成\n")

#### 10. 最终输出 ####
# ==============================================================================
cat("\n")
cat(paste(rep("═", 70), collapse = ""), "\n")
cat("               ANALYSIS SUCCESSFULLY COMPLETED               \n")
cat(paste(rep("═", 70), collapse = ""), "\n\n")

cat("SUMMARY OF KEY OUTPUTS:\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("1. MODULE-TRAIT CORRELATION HEATMAP (Figure 1)\n")
cat("   File: 04_Module_Trait_Correlations/Module_Trait_Correlation_Heatmap_Modern.pdf\n")
cat("   Description: Modern visualization of relationships between modules and treatments\n\n")

cat("2. MAGENTA MODULE ANALYSIS (Figures 2-3)\n")
cat("   Expression Heatmap: 05_Hub_Genes/magenta/Magenta_Hub_Gene_Heatmap_Modern.pdf\n")
cat("   Interaction Network: 05_Hub_Genes/magenta/Magenta_Hub_Gene_Network_Modern.pdf\n")
cat("   Hub Gene List: 05_Hub_Genes/magenta/Magenta_Hub_Genes.csv\n\n")

cat("3. DARKGREEN MODULE ANALYSIS (Figures 4-5)\n")
cat("   Expression Heatmap: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Heatmap_Modern.pdf\n")
cat("   Interaction Network: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Gene_Network_Modern.pdf\n")
cat("   Hub Gene List: 05_Hub_Genes/darkgreen/Darkgreen_Hub_Genes.csv\n\n")

cat("4. KEGG ENRICHMENT ANALYSIS PREPARATION\n")
cat("   Gene Lists: 09_Third_Party_Analysis/KEGG_Input/\n")
cat("   Analysis Guide: 06_Enrichment_Analysis/KEGG_Enrichment_Guide.txt\n")
cat("   Visualization Template: 09_Third_Party_Analysis/Visualization_Template.R\n")
cat("   Note: Use external tools (KEGG Mapper, KOBAS) for actual enrichment analysis\n\n")

cat("5. NETWORK VISUALIZATION FILES\n")
cat("   Cytoscape Files: 07_Cytoscape_Files/\n")
cat("   Compatible with Cytoscape 3.8+ for interactive network exploration\n\n")

cat("6. SUMMARY REPORTS\n")
cat("   Comprehensive Report: 08_Summary_Reports/Comprehensive_Analysis_Summary.txt\n")
cat("   Executive Summary: 08_Summary_Reports/Executive_Summary.txt\n\n")

cat(paste0(rep("=", 70), collapse = ""), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
cat(paste0(rep("=", 70), collapse = ""), "\n\n")

# 保存会话信息
sink("08_Summary_Reports/R_Session_Info.txt")
cat("WGCNA Analysis Session Information\n")
cat("Analysis completed:", Sys.time(), "\n\n")
sessionInfo()
sink()
cat("  ✓ R会话信息已保存: 08_Summary_Reports/R_Session_Info.txt\n")

# 最终清理
gc()
cat("\n所有分析文件已保存完毕。请查阅上述目录中的结果进行后续分析和论文撰写。\n")
