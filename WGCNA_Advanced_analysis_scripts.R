# ==============================================================================
# WGCNA 进阶分析脚本：现代化可视化、富集分析与核心基因挖掘
# 版本：1.0 | 优化日期：2025-12-06
# 特点：现代化图形、交互式选项、发表级输出
# ==============================================================================

# 1. 加载必要的包
options(stringsAsFactors = FALSE)
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(igraph)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(viridis)          # 现代化颜色方案
library(ggsci)           # 科学期刊颜色方案
library(ggrepel)         # 智能标签避免重叠
library(patchwork)       # 图形拼接
library(scales)          # 图形缩放

# --- 富集分析包 ---
library(clusterProfiler)
library(enrichplot)      # 富集分析可视化增强
library(org.Hs.eg.db)

# 设置主题为现代化科学图表
theme_set(theme_minimal(base_size = 11) + 
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
                  axis.line = element_line(color = "black", linewidth = 0.5),
                  axis.title = element_text(face = "bold", size = 12),
                  plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(color = "grey40", hjust = 0.5),
                  legend.position = "right",
                  legend.background = element_rect(fill = "white", color = NA)))

# 设置颜色方案
modern_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", 
                   "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
                   "#7E6148FF", "#B09C85FF")

# 多线程设置，将获取核心数的标准函数替换为 detectCores()
enableWGCNAThreads(nThreads = ifelse(parallel::detectCores() > 4, 6, 4))

# ==============================================================================
# 2. 加载之前分析的保存数据（优化版）
# ==============================================================================
cat("📂 正在加载 Step1 的分析结果...\n")

# 使用tryCatch确保数据加载安全
load_and_validate_data <- function() {
  
  # ================== 1. 加载原始数据 ==================
  cat("📂 正在加载 Step1 的分析结果...\n")
  
  # 加载网络对象 (net, module_colors, MEs, softPower)
  if(!file.exists("03_Network/WGCNA_Network_Object.RData")) {
    stop("❌ 找不到 03_Network/WGCNA_Network_Object.RData，请先运行第一步脚本。")
  }
  load("03_Network/WGCNA_Network_Object.RData")
  
  # 加载预处理后的表达矩阵 (用于热图和计算kME)
  if(!file.exists("01_InputData/Preprocessed_Expression_Matrix.csv")) {
    stop("❌ 找不到 01_InputData/Preprocessed_Expression_Matrix.csv。")
  }
  datExpr <- read.csv("01_InputData/Preprocessed_Expression_Matrix.csv", row.names = 1)
  datExpr <- as.data.frame(t(datExpr)) # 转置为 WGCNA 格式 (行=样本, 列=基因)
  
  # ================== 2. 核心修复：确保 module_colors 有正确的基因名 ==================
  cat("🔧 检查并修复 module_colors 的基因名...\n")
  
  # 情况1: module_colors 完全没名字 (这正是你遇到的情况)
  if (is.null(names(module_colors))) {
    cat("   警告：module_colors 缺少基因名。\n")
    # 假设其顺序与 datExpr 的列名（基因）完全一致
    if (length(module_colors) == ncol(datExpr)) {
      names(module_colors) <- colnames(datExpr)
      cat("   ✅ 已按顺序为其赋予 datExpr 的列名。\n")
    } else {
      stop(paste("错误：module_colors 长度 (", length(module_colors), 
                 ") 与 datExpr 基因数 (", ncol(datExpr), ") 不匹配。请检查脚本-A。"))
    }
  } 
  # 情况2: module_colors 有名字，但与 datExpr 的基因名不匹配 (备用逻辑)
  else if (length(intersect(colnames(datExpr), names(module_colors))) < 100) {
    cat("   警告：module_colors 与 datExpr 的基因名交集很少。\n")
    # 再次检查长度是否一致，如果一致则直接替换名字（假设顺序一致）
    if (length(module_colors) == ncol(datExpr)) {
      cat("   ✅ 长度一致，正在将 module_colors 的名称同步为 datExpr 的列名（假设顺序一致）...\n")
      names(module_colors) <- colnames(datExpr)
    } else {
      warning("⚠️ 无法自动修复命名问题。共同基因数可能仍会很少。")
    }
  }
  
  # ================== 3. 进行交集匹配 ==================
  common_genes <- intersect(colnames(datExpr), names(module_colors))
  cat(paste("   共同基因数:", length(common_genes), "\n"))
  
  if(length(common_genes) < 100) {
    warning("⚠️ 共同基因数较少，请检查数据一致性。")
  }
  
  datExpr <- datExpr[, common_genes, drop = FALSE]
  module_colors <- module_colors[common_genes]
  
  # ================== 4. 最终验证 ==================
  cat(paste("   ✅ 最终 datExpr 维度:", nrow(datExpr), "个样本 x", ncol(datExpr), "个基因\n"))
  if (nrow(datExpr) < 3) {
    stop("❌ 错误：样本数不足3个，无法进行WGCNA分析。")
  }
  if (ncol(datExpr) < 5000) {
    warning("⚠️  基因数较少可能影响网络分析。")
  }
  
  # ================== 5. 返回所有必要数据 ==================
  return(list(datExpr = datExpr, module_colors = module_colors, net = net, MEs = MEs))
}

# ================== 执行函数并分配结果到全局环境 ==================
data_list <- load_and_validate_data()
list2env(data_list, envir = .GlobalEnv)

# 创建输出目录（分级目录）
dir.create("06_Advanced_Results", showWarnings = FALSE)
dir.create("06_Advanced_Results/Figures", showWarnings = FALSE)
dir.create("06_Advanced_Results/Tables", showWarnings = FALSE)
dir.create("06_Advanced_Results/Enrichment", showWarnings = FALSE)

# ==============================================================================
# 3. 现代化基因聚类树和模块分割可视化
# ==============================================================================
cat("🎨 绘制现代化基因聚类树与模块颜色...\n")

# 保存为高质量PDF和PNG
pdf("06_Advanced_Results/Figures/Gene_Dendrogram_with_Modules.pdf", 
    width = 14, height = 8, useDingbats = FALSE)

# 创建自定义颜色标签
color_labels <- module_colors[net$blockGenes[[1]]]
unique_colors <- unique(color_labels)
color_legend <- unique_colors[!is.na(unique_colors)]

plotDendroAndColors(net$dendrograms[[1]], 
                    colors = color_labels,
                    groupLabels = "Module Colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    guideAll = FALSE,
                    main = "Gene Dendrogram and Module Assignment",
                    cex.main = 1.5,
                    cex.colorLabels = 0.8,
                    marAll = c(1, 5, 3, 1))

# 添加图例
legend("topright", 
       legend = paste("Module", color_legend),
       fill = color_legend,
       border = NA,
       bty = "n",
       cex = 0.8,
       title = "Module Legend")

dev.off()

# 同时保存高分辨率PNG
png("06_Advanced_Results/Figures/Gene_Dendrogram_with_Modules.png", 
    width = 2800, height = 1600, res = 300)
plotDendroAndColors(net$dendrograms[[1]], 
                    colors = color_labels,
                    groupLabels = "Module Colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene Dendrogram and Module Assignment",
                    cex.main = 1.5)
dev.off()

cat("✅ 已保存: 06_Advanced_Results/Figures/Gene_Dendrogram_with_Modules.pdf\n")

# ==============================================================================
# 4. 现代化软阈值网络拓扑分析
# ==============================================================================
cat("📊 绘制现代化软阈值网络拓扑图...\n")

# ======= 调试代码开始 =======
cat("\n[调试] 检查 datExpr 维度以定位问题:\n")
cat(paste("  行数 (样本数):", nrow(datExpr), "\n"))
cat(paste("  列数 (基因数):", ncol(datExpr), "\n"))
cat(paste("  前5个样本名:", paste(head(rownames(datExpr), 5), collapse=", "), "\n"))
cat(paste("  前5个基因名:", paste(head(colnames(datExpr), 5), collapse=", "), "\n"))
# ======= 调试代码结束 =======

# 重新计算软阈值统计量
powers <- c(1:10, seq(12, 20, 2))
# 允许WGCNA使用多个CPU核心，启用多线程并行计算，适应低内存电脑
enableWGCNAThreads(nThreads = 4) # 可设为4或6，不要超过电脑的物理核心数

# 然后，在调用 pickSoftThreshold 时，设置更明确的参数
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         networkType = "unsigned",
                         verbose = 5, # 显示详细进度
                         blockSize = 2000) # 可尝试调整（如1000, 2000, 3000）以找到内存和速度最佳平衡点

# 转换为数据框用于ggplot
sft_df <- data.frame(
  Power = sft$fitIndices$Power,
  SFT_R.sq = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  MeanConnectivity = sft$fitIndices$mean.k.,
  MedianConnectivity = sft$fitIndices$median.k.
)

# 创建现代化双面板图
p1 <- ggplot(sft_df, aes(x = Power, y = SFT_R.sq)) +
  geom_line(color = "#3C5488FF", linewidth = 1.2) +
  geom_point(size = 3, color = "#E64B35FF") +
  geom_text(aes(label = Power), vjust = -1, size = 3) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = softPower, linetype = "dashed", color = "blue", alpha = 0.7) +
  labs(title = "Scale-Free Topology Fit",
       subtitle = paste("Selected soft threshold (β):", softPower),
       x = "Soft Threshold (Power)",
       y = expression(paste("Scale Free Topology Model Fit (signed R"^2*")"))) +
  scale_x_continuous(breaks = powers) +
  annotate("text", x = max(powers)*0.8, y = 0.87, 
           label = "R² = 0.85 cutoff", color = "red", size = 3.5)

p2 <- ggplot(sft_df, aes(x = Power, y = MeanConnectivity)) +
  geom_line(color = "#00A087FF", linewidth = 1.2) +
  geom_point(size = 3, color = "#4DBBD5FF") +
  geom_text(aes(label = Power), vjust = -1, size = 3) +
  geom_vline(xintercept = softPower, linetype = "dashed", color = "blue", alpha = 0.7) +
  labs(title = "Mean Connectivity",
       x = "Soft Threshold (Power)",
       y = "Mean Connectivity") +
  scale_x_continuous(breaks = powers) +
  scale_y_log10() +  # 对数尺度更好展示
  annotation_logticks(sides = "l")

# 使用patchwork组合图形
combined_plot <- p1 + p2 + 
  plot_annotation(title = "Network Topology Analysis for Soft Threshold Selection",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave("06_Advanced_Results/Figures/Soft_Threshold_Analysis.pdf", 
       plot = combined_plot, width = 14, height = 6, dpi = 300)
ggsave("06_Advanced_Results/Figures/Soft_Threshold_Analysis.png", 
       plot = combined_plot, width = 14, height = 6, dpi = 300)

cat("✅ 已保存现代化软阈值分析图\n")

# ==============================================================================
# 5. 增强型富集分析（支持KEGG和GO）
# ==============================================================================
cat("🔬 进行增强型富集分析...\n")

enhanced_enrichment_analysis <- function(module_color, gene_list, 
                                         organism_db, org_code, 
                                         top_n = 15) {
  
  cat(paste0("  📈 分析模块: ", module_color, " (", length(gene_list), " 个基因)\n"))
  
  results_list <- list()
  
  tryCatch({
    # ID转换
    gene_entrez <- bitr(gene_list, fromType = "SYMBOL", 
                        toType = c("ENTREZID", "ENSEMBL", "GENENAME"), 
                        OrgDb = organism_db)
    
    if(nrow(gene_entrez) < 10) {
      cat(paste0("    ⚠️  基因数太少 (", nrow(gene_entrez), ")，跳过富集分析\n"))
      return(NULL)
    }
    
    # 1. KEGG富集分析
    kk <- enrichKEGG(gene = gene_entrez$ENTREZID,
                     organism = org_code,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     minGSSize = 5,
                     maxGSSize = 500)
    
    if(!is.null(kk) && nrow(kk) > 0) {
      # 保存KEGG结果
      write.csv(kk@result, 
                paste0("06_Advanced_Results/Enrichment/KEGG_", module_color, ".csv"))
      
      # 创建增强型气泡图
      if(nrow(kk) > 1) {
        p_kegg <- dotplot(kk, showCategory = top_n, 
                          color = "qvalue", 
                          size = "Count",
                          label_format = 40) +  # 标签换行
          scale_color_viridis_c(direction = -1) +
          labs(title = paste("KEGG Enrichment:", module_color),
               subtitle = paste(length(gene_entrez$ENTREZID), "genes with Entrez ID")) +
          theme(axis.text.y = element_text(size = 10))
        
        ggsave(paste0("06_Advanced_Results/Figures/KEGG_", module_color, ".pdf"), 
               plot = p_kegg, width = 10, height = 8)
      }
      results_list$KEGG <- kk
    }
    
    # 2. GO富集分析（生物过程）
    ego_bp <- enrichGO(gene = gene_entrez$ENTREZID,
                       OrgDb = organism_db,
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
    
    if(!is.null(ego_bp) && nrow(ego_bp) > 0) {
      write.csv(ego_bp@result, 
                paste0("06_Advanced_Results/Enrichment/GO_BP_", module_color, ".csv"))
      
      # 简化GO术语
      ego_bp_simplify <- simplify(ego_bp, cutoff = 0.7)
      
      if(nrow(ego_bp_simplify) > 0) {
        p_go <- dotplot(ego_bp_simplify, showCategory = top_n,
                        color = "qvalue", size = "Count") +
          scale_color_viridis_c(direction = -1) +
          labs(title = paste("GO Biological Process:", module_color))
        
        ggsave(paste0("06_Advanced_Results/Figures/GO_BP_", module_color, ".pdf"), 
               plot = p_go, width = 12, height = 9)
      }
      results_list$GO_BP <- ego_bp
    }
    
    return(results_list)
    
  }, error = function(e) {
    cat(paste0("    ❌ 富集分析失败: ", e$message, "\n"))
    return(NULL)
  })
}

# 自动识别所有显著模块进行分析
module_summary <- as.data.frame(table(module_colors)) %>%
  arrange(desc(Freq)) %>%
  filter(Freq >= 30)  # 只分析基因数大于30的模块

cat(paste("📋 将分析", nrow(module_summary), "个模块\n"))

for(i in 1:nrow(module_summary)) {
  module_color <- as.character(module_summary$module_colors[i])
  gene_count <- module_summary$Freq[i]
  
  genes_in_module <- names(module_colors)[module_colors == module_color]
  
  cat(paste0("  [", i, "/", nrow(module_summary), "] ", 
             module_color, " (", gene_count, " genes)\n"))
  
  enhanced_enrichment_analysis(module_color, genes_in_module, 
                               org.Hs.eg.db, "hsa")
}

# ==============================================================================
# 6. 现代化核心基因互作网络与热图分析
# ==============================================================================
cat("🕸️ 绘制现代化核心基因网络与热图...\n")

analyze_module_hub_modern <- function(target_module, datExpr, MEs, 
                                      n_top_genes = 25, 
                                      correlation_cutoff = 0.6) {
  
  cat(paste0("  🔍 深入分析模块: ", target_module, "\n"))
  
  # 1. 筛选模块基因
  inModule <- module_colors == target_module
  modGenes <- names(module_colors)[inModule]
  
  if(length(modGenes) < 10) {
    cat("    ⚠️ 基因数量太少，跳过\n")
    return(NULL)
  }
  
  # 2. 计算kME
  ME_name <- paste0("ME", target_module)
  if(!ME_name %in% colnames(MEs)) {
    cat("    ℹ️ 重新计算模块特征基因...\n")
    ME_temp <- moduleEigengenes(datExpr[, modGenes], 
                                colors = rep(target_module, length(modGenes)))$eigengenes
    datKME <- signedKME(datExpr[, modGenes], ME_temp, outputColumnName = "kME")
  } else {
    datKME <- signedKME(datExpr[, modGenes], MEs[, ME_name, drop = FALSE], 
                        outputColumnName = "kME")
  }
  
  # 3. 筛选核心基因
  kME_col_name <- colnames(datKME)[1]
  top_genes_df <- data.frame(
    Gene = rownames(datKME),
    kME = datKME[, kME_col_name],
    Module = target_module
  ) %>%
    arrange(desc(abs(kME))) %>%
    mutate(Rank = 1:nrow(.),
           Hub_Gene = ifelse(Rank <= n_top_genes, TRUE, FALSE))
  
  # 保存核心基因列表
  write.csv(top_genes_df, 
            paste0("06_Advanced_Results/Tables/Hub_Genes_", target_module, ".csv"))
  
  top_genes <- top_genes_df$Gene[1:min(n_top_genes, nrow(top_genes_df))]
  
  # 4. 现代化热图
  top_expr <- datExpr[, top_genes, drop = FALSE]
  top_expr_scaled <- t(scale(top_expr))
  
  # 创建注释信息
  annotation_col <- data.frame(
    Module = rep(target_module, ncol(top_expr_scaled))
  )
  rownames(annotation_col) <- colnames(top_expr_scaled)
  
  # 关键修复：为当前模块动态创建颜色映射
  annotation_colors <- list(
    Module = setNames(target_module, target_module)
  )
  
  # 确保因子水平匹配
  annotation_col$Module <- factor(annotation_col$Module, 
                                  levels = names(annotation_colors$Module)
                                  )
  
  p_heatmap <- pheatmap(top_expr_scaled,
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        show_colnames = ncol(top_expr_scaled) <= 30,
                        show_rownames = TRUE,
                        annotation_col = annotation_col,
                        annotation_colors = annotation_colors,
                        color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                        fontsize_row = 9,
                        fontsize_col = 8,
                        border_color = NA,
                        main = paste("Expression Pattern of Top", length(top_genes), 
                                     "Hub Genes in", target_module),
                        filename = paste0("06_Advanced_Results/Figures/Heatmap_", 
                                          target_module, ".pdf"),
                        width = 12,
                        height = 10)
  
  # 5. 现代化网络图
  if(length(top_genes) >= 5) {
    # 计算 Top 基因之间的相关性作为边权重
    cor_matrix <- cor(top_expr, use = "pairwise.complete.obs", method = "pearson")
    
    # 构建 Graph 对象
    adj_matrix <- cor_matrix
    # 关键修改1：将相关性取绝对值，使权重为正
    adj_matrix <- abs(adj_matrix)
    # 关键修改2：过滤弱连接，只保留强相关 (例如保留相关性绝对值>0.5的边)
    adj_matrix[adj_matrix < 0.5] <- 0
    diag(adj_matrix) <- 0
    
    # 构建 Graph 对象
    g <- graph_from_adjacency_matrix(adj_matrix,
                                     mode = "undirected",
                                     weighted = TRUE,
                                     diag = FALSE)
    
    # 如果图中没有边，则跳过绘图
    if (ecount(g) == 0) {
      cat(paste0("    提示: 模块 ", target_module, " 的核心基因间无强连接（相关性绝对值均<0.5），跳过网络图绘制。\n"))
      return(top_genes_df) # 提前返回，不执行后续绘图代码
    }
    
    # 设置节点属性（保持不变）
    V(g)$size <- 15 + 10 * scale(abs(top_genes_df$kME[match(V(g)$name, top_genes_df$Gene)]))[,1]
    # 关键修改3：根据原始相关性的正负为边上色，保留生物学意义
    edge_weights <- E(g)$weight # 此时weight是相关性的绝对值
    edge_signs <- sapply(E(g), function(e) {
      sign(cor_matrix[ends(g, e)[1], ends(g, e)[2]]) # 获取原始相关性的正负号
    })
    V(g)$color <- ifelse(V(g)$size > median(V(g)$size),
                         alpha("#E64B35FF", 0.8),
                         alpha("#4DBBD5FF", 0.6))
    V(g)$label.cex <- 0.8 + V(g)$size/50
    V(g)$label.color <- "black"
    V(g)$frame.color <- target_module  # 使用模块颜色作为边框色，或者使用 "gray50" 等固定颜色
    
    # 设置边属性
    E(g)$width <- edge_weights * 5 # 边的粗细基于相关性绝对值
    E(g)$color <- ifelse(edge_signs > 0,
                         alpha("#00A087FF", 0.7), # 正相关用绿色
                         alpha("#F39B7FFF", 0.7)) # 负相关用橙色
    
    # 布局算法 - 现在权重全为正，不会报错
    layout <- layout_with_fr(g)
    
    # 绘制网络图
    pdf(paste0("06_Advanced_Results/Figures/Network_", target_module, ".pdf"), 
        width = 12, height = 10)
    
    par(mar = c(1, 1, 3, 1))
    plot(g,
         layout = layout,
         main = paste("Hub Gene Interaction Network:", target_module),
         vertex.label = ifelse(V(g)$size > median(V(g)$size), 
                               V(g)$name, ""),  # 只标注大节点
         vertex.label.dist = 0.5,
         vertex.label.font = 2,
         edge.curved = 0.2)
    
    # 添加图例
    legend("bottomright",
           legend = c("High kME", "Low kME", "Positive cor", "Negative cor"),
           col = c(alpha("#E64B35FF", 0.8), alpha("#4DBBD5FF", 0.6),
                   alpha("#00A087FF", 0.5), alpha("#F39B7FFF", 0.5)),
           pch = c(19, 19, NA, NA),
           lty = c(NA, NA, 1, 1),
           lwd = c(NA, NA, 3, 3),
           bty = "n",
           cex = 0.9,
           title = "Legend")
    
    dev.off()
    
    cat(paste0("    ✅ 网络图已保存 (", vcount(g), " 节点, ", ecount(g), " 边)\n"))
  }
  
  # 6. 核心基因kME分布图
  p_dist <- ggplot(top_genes_df, aes(x = Rank, y = abs(kME), color = Hub_Gene)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line(data = top_genes_df[top_genes_df$Hub_Gene, ], 
              aes(x = Rank, y = abs(kME)), 
              color = "#E64B35FF", linewidth = 0.5) +
    scale_color_manual(values = c("FALSE" = "#8491B4FF", "TRUE" = "#E64B35FF")) +
    labs(title = paste("Module Membership (kME) Distribution:", target_module),
         x = "Gene Rank",
         y = "|Module Membership (kME)|",
         color = "Hub Gene") +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = mean(abs(top_genes_df$kME)), 
               linetype = "dashed", color = "grey50")
  
  ggsave(paste0("06_Advanced_Results/Figures/kME_Distribution_", target_module, ".pdf"),
         plot = p_dist, width = 10, height = 6)
  
  return(top_genes_df)
}

# 自动分析所有显著模块
for(module_color in as.character(module_summary$module_colors[1:5])) {  # 分析前5大模块
  analyze_module_hub_modern(module_color, datExpr, MEs, n_top_genes = 20)
}

# ==============================================================================
# 7. 模块相关性热图（新增）
# ==============================================================================
cat("📈 绘制模块相关性热图...\n")

if(exists("MEs") && !is.null(MEs)) {
  # 计算模块特征基因相关性
  MEs_clean <- MEs[, !grepl("^MEgrey$", colnames(MEs))]  # 去除grey模块
  module_cor <- cor(MEs_clean, use = "pairwise.complete.obs")
  module_pheatmap <- pheatmap(module_cor,
                              color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
                              border_color = NA,
                              display_numbers = TRUE,
                              number_format = "%.2f",
                              number_color = "black",
                              fontsize_number = 8,
                              main = "Module-Trait Relationships",
                              filename = "06_Advanced_Results/Figures/Module_Correlation_Heatmap.pdf",
                              width = 12,
                              height = 10)
}

# ==============================================================================
# 8. 生成分析报告摘要
# ==============================================================================
cat("📄 生成分析报告摘要...\n")

generate_summary_report <- function(module_colors, datExpr) {
  summary_stats <- data.frame(
    Module = names(table(module_colors)),
    Gene_Count = as.numeric(table(module_colors)),
    Proportion = round(as.numeric(table(module_colors))/length(module_colors)*100, 2)
  ) %>%
    arrange(desc(Gene_Count))
  
  write.csv(summary_stats, "06_Advanced_Results/Tables/Module_Summary_Statistics.csv")
  
  # 创建可视化摘要
  p_summary <- ggplot(summary_stats, aes(x = reorder(Module, Gene_Count), y = Gene_Count, fill = Module)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = Gene_Count), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = setNames(summary_stats$Module, summary_stats$Module)) +
    coord_flip() +
    labs(title = "WGCNA Module Size Distribution",
         x = "Module",
         y = "Number of Genes") +
    theme(legend.position = "none",
          axis.text.y = element_text(color = summary_stats$Module, face = "bold"))
  
  ggsave("06_Advanced_Results/Figures/Module_Size_Distribution.pdf",
         plot = p_summary, width = 10, height = 8)
  
  return(summary_stats)
}

module_summary <- generate_summary_report(module_colors, datExpr)

# ==============================================================================
# 9. 完成信息
# ==============================================================================
cat(paste("\n", strrep("=", 60), "\n", sep=""))
cat("🎉 WGCNA进阶分析完成！\n")
cat(strrep("-", 60) + "\n")
cat("📁 结果目录结构:\n")
cat("  06_Advanced_Results/\n")
cat("  ├── Figures/          # 所有图形文件 (PDF & PNG)\n")
cat("  ├── Tables/           # 数据表格 (CSV格式)\n")
cat("  └── Enrichment/       # 富集分析结果\n")
cat(paste(strrep("-", 60), "\n", sep=""))
cat("📊 分析统计:\n")
cat(paste("  • 总基因数:", length(module_colors), "\n"))
cat(paste("  • 模块数量:", length(unique(module_colors)), "\n"))
cat(paste("  • 最大模块:", module_summary$Module[1], 
          "(", module_summary$Gene_Count[1], "genes,", 
          module_summary$Proportion[1], "%)\n"))
cat(paste("  • 最小模块:", tail(module_summary$Module, 1), 
          "(", tail(module_summary$Gene_Count, 1), "genes)\n"))
cat(paste(strrep("=", 60), "\n", sep=""))

# 保存会话信息
writeLines(capture.output(sessionInfo()), 
           "06_Advanced_Results/Analysis_Session_Info.txt")