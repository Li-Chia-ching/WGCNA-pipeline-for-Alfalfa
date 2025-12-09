# ==============================================================================
# Medicago sativa WGCNA 完整分析流程 v7.0-Final (修复版)
# 硬件适配: AMD Ryzen 7 PRO/ 15GB RAM
# 功能: 全面可视化 + 特定模块分析 + KEGG引导
# ==============================================================================

#### 1. 环境与配置 ####
# ==============================================================================
cat("\n[Step 1] Environment Initialization...\n")
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 加载必要包
pkgs <- c("WGCNA", "dplyr", "stringr", "ggplot2", "pheatmap", "igraph", 
          "ComplexHeatmap", "circlize", "RColorBrewer", "viridis")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# 并行设置
n_cores <- max(1, parallel::detectCores() - 2)
enableWGCNAThreads(nThreads = n_cores)

# 创建完整目录结构
dirs <- c("01_Data_Clean", "02_Network_Topology", "03_Modules", 
          "04_Module_Trait_Correlations", "05_Hub_Genes", 
          "05_Hub_Genes/magenta", "05_Hub_Genes/darkgreen", # 预留特定目录
          "06_Enrichment_Materials", "08_Summary_Reports")
for(d in dirs) if(!dir.exists(d)) dir.create(d, recursive=TRUE)

# 现代化配色
modern_pal <- list(
  net_node = "#2C3E50", net_edge = "grey80",
  heat_map = colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
)

#### 2. 数据清洗 (针对 15GB RAM 优化 - 增强版) ####
# ==============================================================================
cat("\n[Step 2] Data Processing (Robust Filtering & Mapping)...\n")

# 2.1 读取数据
if (!exists("Expression_with_annotation")) stop("Missing: Expression_with_annotation")
if (!exists("Transcriptome_grouping_summary")) stop("Missing: Transcriptome_grouping_summary")

# 确保数据是 data.frame 格式（防止 tibble 导致行名设置失败）
raw_expr <- as.data.frame(Expression_with_annotation)
meta_data <- as.data.frame(Transcriptome_grouping_summary)

#### 2.2 格式化表达矩阵 (针对混合列名的修复版) ####
# ==============================================================================
cat("\n[Step 2.2] Formatting Expression Matrix (Extracting FPKM/TPM)...\n")

# 1. 提取基因ID (第一列)
gene_ids <- as.character(raw_expr[, 1]) 

# 2. 识别 FPKM 或 TPM 列
# 根据您的 Debug Info，列名包含 ":fpkm" 或 ":read count"
all_cols <- colnames(raw_expr)
fpkm_cols_idx <- grep(":fpkm$", all_cols, ignore.case = TRUE)

# 如果没找到 FPKM，尝试找 TPM
if(length(fpkm_cols_idx) == 0) {
  fpkm_cols_idx <- grep(":tpm$", all_cols, ignore.case = TRUE)
}

# 3. 检查是否找到
if(length(fpkm_cols_idx) == 0) {
  # 如果都没有，打印前几个列名帮助诊断
  cat("错误：无法识别表达量列。当前列名示例：\n")
  print(head(all_cols, 10))
  stop("无法找到以 ':fpkm' 或 ':tpm' 结尾的列！请检查 Expression_with_annotation 数据格式。")
}

cat(sprintf("  -> 识别到 %d 个 FPKM/TPM 数据列，正在提取...\n", length(fpkm_cols_idx)))

# 4. 提取数据构建矩阵
datExpr0 <- raw_expr[, fpkm_cols_idx, drop = FALSE]

# 5. 清洗列名：去除 ":fpkm" 或 ":tpm" 后缀
# 这样 "BPGS_1:fpkm" 就会变成 "BPGS_1"，就能和分组表匹配了
current_names <- colnames(datExpr0)
clean_names <- gsub(":(fpkm|tpm)$", "", current_names, ignore.case = TRUE)
colnames(datExpr0) <- clean_names

# 6. 设置行名 (基因ID)
valid_gene_names <- make.names(gene_ids, unique = TRUE)

# 保存映射表 (Important)
gene_name_mapping <- data.frame(
  Original_ID = gene_ids,
  Valid_Rowname = valid_gene_names,
  stringsAsFactors = FALSE
)
write.csv(gene_name_mapping, "01_Data_Clean/Gene_Name_Mapping_Full.csv", row.names = FALSE)

rownames(datExpr0) <- valid_gene_names
cat(sprintf("  -> 矩阵构建完成: %d 基因 x %d 样本 (已去除后缀)\n", nrow(datExpr0), ncol(datExpr0)))


#### 2.3 样本对齐 (基于清洗后的列名) ####
# ==============================================================================
cat("\n[Step 2.3] Sample Alignment...\n")

# 1. 获取分组表中的 Sample 列
# 根据 Debug Info，您的列名确实叫 "Sample"，且第一行是 "BPGS_1"
sample_col_idx <- grep("^Sample$", colnames(meta_data), ignore.case = TRUE)[1]
if(is.na(sample_col_idx)) stop("Error: Grouping summary missing 'Sample' column.")

meta_samples <- as.character(meta_data[[sample_col_idx]])
meta_samples <- trimws(gsub('^"|"$', '', meta_samples))

# 2. 获取表达矩阵列名 (现在已经是清洗过的 "BPGS_1" 等)
expr_samples <- colnames(datExpr0)

# 3. 执行匹配
common_samples <- intersect(expr_samples, meta_samples)

# 4. 诊断与报错
if(length(common_samples) == 0) {
  cat("\n=======================================================\n")
  cat("【CRITICAL ERROR DEBUG INFO】\n")
  cat("Cleaned Expr Columns: ", paste(head(expr_samples, 5), collapse=", "), "\n")
  cat("Meta Sample IDs:      ", paste(head(meta_samples, 5), collapse=", "), "\n")
  cat("=======================================================\n")
  stop("无法匹配样本！请对比上方'Cleaned Expr Columns'和'Meta Sample IDs'是否一致。")
}

# 5. 对齐数据
datExpr0 <- datExpr0[, common_samples, drop = FALSE]
meta_data$SampleID_Clean <- meta_samples # 确保有一列标准ID
meta_data <- meta_data[match(common_samples, meta_data$SampleID_Clean), ]

cat(sprintf("  -> 成功对齐样本数: %d\n", ncol(datExpr0)))

# 2.4 关键过滤 (适配 15GB RAM)
# ------------------------------------------------------------------
cat("  -> Filtering genes to fit 15GB RAM limit...\n")

# A. 丰度过滤 (保留在 >20% 样本中 TPM/FPKM >= 1 的基因)
keep_abd <- rowSums(datExpr0 >= 1) >= (0.2 * ncol(datExpr0))
datExpr_abd <- datExpr0[keep_abd, , drop = FALSE]
cat(sprintf("     Abundance filter: %d genes removed.\n", nrow(datExpr0) - nrow(datExpr_abd)))

# B. 方差过滤 (保留 Top 25,000 高变异基因)
target_n <- 25000
if(nrow(datExpr_abd) > target_n){
  mads <- apply(datExpr_abd, 1, mad)
  # 获取保留基因的索引
  keep_var_idx <- order(mads, decreasing = TRUE)[1:target_n]
  datExpr_final <- datExpr_abd[keep_var_idx, , drop = FALSE]
  cat(sprintf("     Variance filter: Kept top %d genes (out of %d).\n", target_n, nrow(datExpr_abd)))
} else {
  datExpr_final <- datExpr_abd
  cat("     Variance filter: Gene count below 25,000. Kept all.\n")
}

# [增强] 保存最终分析用的基因映射表
# 这一步非常重要：后续找 Hub 基因时，你手里只有 Valid_Rowname，
# 这个表能帮你快速查回 Original_ID
final_genes <- rownames(datExpr_final)
final_mapping <- gene_name_mapping[gene_name_mapping$Valid_Rowname %in% final_genes, ]
write.csv(final_mapping, "01_Data_Clean/Gene_Name_Mapping_Filtered.csv", row.names = FALSE)
cat("  -> Filtered gene mapping saved: 01_Data_Clean/Gene_Name_Mapping_Filtered.csv\n")

#### 2.5 转置并取 log2 (标准化) ####
# ==============================================================================
# 接上文 Step 2.4 ...
datExpr <- as.data.frame(t(log2(datExpr_final + 1)))
cat(sprintf("  -> Final Matrix for WGCNA: %d Samples x %d Genes\n", nrow(datExpr), ncol(datExpr)))

#### 2.6 构建性状矩阵 (修复: 必须在绘图前完成) ####
# ==============================================================================
cat("\n[Step 2.6] Constructing Trait Matrix...\n")
datTraits <- data.frame(row.names = rownames(datExpr))
groups <- unique(meta_data$Group_Short_Name)

for(g in groups) {
  # 将分组名转换为合法的变量名 (例如 "High-Temp" -> "High_Temp")
  safe_g <- gsub("-", "_", g)
  datTraits[[safe_g]] <- ifelse(meta_data$Group_Short_Name == g, 1, 0)
}
cat("  -> Trait matrix constructed.\n")

#### 2.7 样本聚类、离群值检测与增强可视化 (修复依赖后) ####
# ==============================================================================
cat("\n[Step 2.7] Sample Dendrogram with Dynamic Marking & Validation...\n")

# 1. 计算样本间距离与 Z-score 检测
# ------------------------------------------------------------------
sample_dist_matrix <- 1 - cor(t(datExpr), use = "p")
sampleTree <- hclust(as.dist(sample_dist_matrix), method = "average")

# Z-score 离群值判定 (阈值 2.5)
mean_dist_per_sample <- colMeans(sample_dist_matrix)
z_scores <- scale(mean_dist_per_sample)
outlier_threshold <- 2.5 
outliers <- rownames(z_scores)[which(z_scores > outlier_threshold)]

# 2. 准备绘图颜色 (现在 datTraits 已经存在了，不会报错)
traitColors <- numbers2colors(datTraits, signed = FALSE)

# 3. 输出主图：聚类树 + 性状热图 + 动态标记
# ------------------------------------------------------------------
pdf("01_Data_Clean/Sample_Dendrogram_Trait_Heatmap_QC.pdf", width = 14, height = 10)
par(mar = c(1, 4, 3, 1))

# A. 绘制基础树图
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample Dendrogram and Trait Heatmap (With QC)",
                    cex.dendroLabels = 0.6,
                    marAll = c(1, 4, 3, 1))

# B. 动态添加离群点标记 (防止重叠)
if(length(outliers) > 0) {
  # 获取样本在图中的 X 轴坐标
  labels_in_plot_order <- sampleTree$labels[sampleTree$order]
  outlier_x_coords <- match(outliers, labels_in_plot_order)
  
  # 创建临时数据框用于排序处理
  marker_data <- data.frame(x = outlier_x_coords, sample = outliers)
  marker_data <- marker_data[order(marker_data$x), ] # 按位置从左到右排序
  
  # 获取树的高度范围，设定基础标记高度
  y_max <- max(sampleTree$height)
  base_y <- y_max * 0.02 # 基础高度设为树高的2%
  
  # 初始化 Y 坐标
  marker_data$y <- rep(base_y, nrow(marker_data))
  
  # 动态错位逻辑
  if(nrow(marker_data) > 1) {
    for(i in 2:nrow(marker_data)) {
      if((marker_data$x[i] - marker_data$x[i-1]) < 2) {
        if(marker_data$y[i-1] > base_y) {
          marker_data$y[i] <- base_y
        } else {
          marker_data$y[i] <- base_y + (y_max * 0.04) # 抬高
        }
      }
    }
  }
  
  # 绘制红点
  points(x = marker_data$x, y = marker_data$y, 
         pch = 19, col = "red", cex = 1.5)
  
  # 绘制星号
  text(x = marker_data$x, y = marker_data$y, 
       labels = "*", pos = 3, col = "red", cex = 2, font = 2)
  
  # 添加图例
  legend("topright", legend = c("Potential Outlier", "Normal Sample"),
         pch = c(19, NA), col = c("red", NA), 
         text.col = c("red", "black"), 
         bty = "n", cex = 0.9, title = "QC Status")
  
  mtext(paste("Warning:", length(outliers), "outlier(s) detected"), 
        side = 3, line = -2, col = "red", font = 2, adj = 0.98)
}

dev.off()
cat("  -> Saved: 01_Data_Clean/Sample_Dendrogram_Trait_Heatmap_QC.pdf\n")

# 4. 输出验证图：Z-score 分布与散点图
# ------------------------------------------------------------------
if(length(outliers) > 0) {
  pdf("01_Data_Clean/Outlier_Validation_Plot.pdf", width = 12, height = 6)
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
  
  # 左图
  hist(z_scores, breaks = 30, col = "skyblue", border = "white",
       main = "Distribution of Sample Distance Z-scores", 
       xlab = "Z-score")
  abline(v = outlier_threshold, col = "red", lty = 2, lwd = 2)
  
  # 右图
  plot_colors <- ifelse(z_scores > outlier_threshold, "red", "gray50")
  plot_pch <- ifelse(z_scores > outlier_threshold, 19, 1)
  
  plot(mean_dist_per_sample, z_scores, 
       pch = plot_pch, col = plot_colors, cex = 1.2,
       main = "Outlier Validation: Distance vs Z-score",
       xlab = "Mean Distance", ylab = "Z-score")
  abline(h = outlier_threshold, col = "red", lty = 2)
  
  # 标记名称
  text_y_offset <- (max(z_scores) - min(z_scores)) * 0.02
  if(length(outliers) > 0) {
    text(mean_dist_per_sample[outliers], z_scores[outliers] + text_y_offset, 
         labels = outliers, col = "red", cex = 0.8, font = 2)
  }
  
  dev.off()
  cat("  -> Saved Validation Plot: 01_Data_Clean/Outlier_Validation_Plot.pdf\n")
}

# 5. 保存离群点列表
if(length(outliers) > 0) {
  cat("\n  ! WARNING: Potential outlier samples detected (Z-score > 2.5):\n")
  cat(paste("     -", outliers), sep = "\n")
  outlier_info <- data.frame(Sample = outliers, Z_Score = z_scores[outliers,], Mean_Dist = mean_dist_per_sample[outliers])
  write.csv(outlier_info, "01_Data_Clean/Potential_Outliers.csv", row.names = FALSE)
} else {
  cat("  -> QC Passed: No significant outliers detected.\n")
}

#### 2.8 保存清洗后的数据 (Checkpoint) ####
# ==============================================================================
save(datExpr, datTraits, meta_data, gene_name_mapping, 
     file = "01_Data_Clean/Input_Data.RData")

cat("  -> [Checkpoint] Processed data saved to: 01_Data_Clean/Input_Data.RData\n")

# 4. 输出验证图：Z-score 分布与散点图 (新增)
# ------------------------------------------------------------------
if(length(outliers) > 0) {
  pdf("01_Data_Clean/Outlier_Validation_Plot.pdf", width = 12, height = 6)
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
  
  # 左图：Z-score 分布直方图
  hist(z_scores, breaks = 30, col = "skyblue", border = "white",
       main = "Distribution of Sample Distance Z-scores", 
       xlab = "Z-score (Mean Distance to others)")
  abline(v = outlier_threshold, col = "red", lty = 2, lwd = 2)
  legend("topright", legend = paste("Threshold =", outlier_threshold), 
         col = "red", lty = 2, bty = "n")
  
  # 右图：距离 vs Z-score 散点图
  plot_colors <- ifelse(z_scores > outlier_threshold, "red", "gray50")
  plot_pch <- ifelse(z_scores > outlier_threshold, 19, 1)
  
  plot(mean_dist_per_sample, z_scores, 
       pch = plot_pch, col = plot_colors, cex = 1.2,
       main = "Outlier Validation: Distance vs Z-score",
       xlab = "Mean Distance to Other Samples", ylab = "Z-score")
  abline(h = outlier_threshold, col = "red", lty = 2)
  
  # 为离群点添加文字标签 (避免重叠)
  if(requireNamespace("ggrepel", quietly = TRUE)) {
    # 基础绘图如果不使用 ggplot2，用 text 简单标记
    text(mean_dist_per_sample[outliers], z_scores[outliers], 
         labels = outliers, pos = 4, col = "red", cex = 0.8, font = 2)
  } else {
    text(mean_dist_per_sample[outliers], z_scores[outliers], 
         labels = outliers, pos = 4, col = "red", cex = 0.8)
  }
  
  dev.off()
  cat("  -> Saved Validation Plot: 01_Data_Clean/Outlier_Validation_Plot.pdf\n")
}

# 5. 保存离群点数据
if(length(outliers) > 0) {
  cat("\n  ! WARNING: Potential outlier samples detected (Z-score > 2.5):\n")
  cat(paste("     -", outliers), sep = "\n")
  
  outlier_info <- data.frame(
    Sample = outliers,
    Z_Score = z_scores[outliers, ],
    Mean_Distance = mean_dist_per_sample[outliers]
  )
  write.csv(outlier_info, "01_Data_Clean/Potential_Outliers.csv", row.names = FALSE)
  cat("  -> List saved to: 01_Data_Clean/Potential_Outliers.csv\n")
} else {
  cat("  -> QC Passed: No significant outliers detected.\n")
}

#### 3. 网络拓扑分析 (找回软阈值筛选图) ####
# ==============================================================================
cat("\n[Step 3] Soft Thresholding & Plotting...\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
sft_power <- sft$powerEstimate
if(is.na(sft_power)) sft_power <- 9 # 稳健默认值

# === 输出位置 1: 软阈值筛选图 ===
pdf("02_Network_Topology/Soft_Threshold_Plots.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
# SFT index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
abline(h=0.85,col="red")
# Mean Connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()
cat("  -> Saved: 02_Network_Topology/Soft_Threshold_Plots.pdf\n")

#### 4. 模块识别 (综合优化与增强统计版) ####
# ==============================================================================
cat("\n[Step 4] Module Construction (Robust & Optimized)...\n")
start_time <- Sys.time()

# 4.1 智能降维 (适配 15GB RAM)
# ------------------------------------------------------------------
# 目标：聚焦于变异度最高的 20,000 个基因以保证计算速度和内存安全
target_gene_size <- 20000 
if (ncol(datExpr) > target_gene_size) {
  cat(sprintf("  -> Optimizing input: Focusing on top %d most variable genes...\n", target_gene_size))
  gene_mads <- apply(datExpr, 2, mad)
  keep_idx <- order(gene_mads, decreasing = TRUE)[1:target_gene_size]
  
  # [关键] 更新全局 datExpr，确保后续分析基于同一基因集
  datExpr <- datExpr[, keep_idx] 
  
  # 同步更新映射表（如果存在）
  if(exists("gene_name_mapping")) {
    gene_name_mapping <- gene_name_mapping[gene_name_mapping$Valid_Rowname %in% colnames(datExpr), ]
  }
}
cat(sprintf("  -> Analysis Universe: %d Genes x %d Samples\n", ncol(datExpr), nrow(datExpr)))

# 4.2 构建网络 (初始参数: 倾向于大模块)
# ------------------------------------------------------------------
# minModuleSize=120: 避免产生过多细碎模块
# mergeCutHeight=0.25: 初始合并阈值 (相似度 0.75)
net <- blockwiseModules(
  datExpr,
  power = sft_power,
  maxBlockSize = 25000,     # 一次性计算所有基因，避免分块造成的边缘效应
  TOMType = "unsigned",
  minModuleSize = 120,      
  mergeCutHeight = 0.25,    
  deepSplit = 2,            
  numericLabels = TRUE,
  saveTOMs = FALSE,         # 内存计算更快，减少磁盘IO
  verbose = 2,
  nThreads = n_cores
)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
n_modules_init <- length(unique(moduleColors)) - 1 # 排除grey

# 4.3 模块数量检查与动态合并 (处理 n > 15 的情况)
# ------------------------------------------------------------------
cat(sprintf("\n  -> Initial Module Count: %d (excluding grey)\n", n_modules_init))

if (n_modules_init > 15) {
  cat("  -> Module count > 15. Attempting to merge closely related modules...\n")
  
  # 计算模块特征基因的聚类树
  MEs_temp <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
  MEDiss <- 1 - cor(MEs_temp, use = "p")
  METree <- hclust(as.dist(MEDiss), method = "average")
  
  # 设定稍微激进的阈值 0.30 (相似度 0.70)
  # 警告：这可能会合并生物学功能略有差异的模块，但在追求概览图时是可接受的
  merge_cut_force <- 0.30 
  cat(sprintf("  -> Applying merge threshold: %.2f (Correlation > %.2f)\n", 
              merge_cut_force, 1-merge_cut_force))
  
  merge_res <- mergeCloseModules(datExpr, moduleColors, cutHeight = merge_cut_force, verbose = 0)
  
  # 更新结果
  moduleColors <- merge_res$colors
  MEs <- merge_res$newMEs
  n_modules_final <- length(unique(moduleColors)) - 1
  
  # 更新 net 对象
  net$MEs <- MEs
  # 注意：这里暂时不更新 net$colors 为数字，因为 mergeCloseModules 返回的是颜色
  # 我们将在 4.4 统一生成新的数字标签
  
  cat(sprintf("  -> Merge complete. Final Module Count: %d\n", n_modules_final))
} else {
  cat("  -> Module count is within acceptable range. No forced merging needed.\n")
}

# 4.4 生成稳健的数字标签 (修复: Grey映射问题)
# ------------------------------------------------------------------
# 使用您建议的稳健映射法，确保 Grey = 0
# standardColors(200) 提供了足够多的颜色池
all_possible_colors <- c("grey", standardColors(200))
# 匹配颜色到索引，然后减1 (使得 grey 对应的 1 变为 0)
moduleLabels <- match(moduleColors, all_possible_colors) - 1

# 检查是否有匹配失败的情况 (NA)
if(any(is.na(moduleLabels))) {
  cat("  ! Warning: Some colors could not be mapped to labels. Using WGCNA default re-labeling.\n")
  # 备用方案：标准 WGCNA 转换
  moduleLabels <- as.numeric(factor(moduleColors, levels = c("grey", unique(moduleColors[moduleColors!="grey"])))) - 1
}

# 更新 net 对象中的 colors
net$colors <- moduleLabels

# 4.5 统计与质量控制 (新增建议功能)
# ------------------------------------------------------------------
# A. 模块大小分布
module_sizes <- table(moduleColors)
cat("\n  -> Module Size Distribution:\n")
# 排序输出：先大后小 (排除grey)
sorted_mods <- names(sort(module_sizes[names(module_sizes) != "grey"], decreasing = TRUE))
for(mod in sorted_mods) {
  cat(sprintf("     - %-12s: %d genes\n", mod, module_sizes[mod]))
}

# B. 灰色模块分析
grey_size <- sum(moduleColors == "grey")
grey_prop <- (grey_size / length(moduleColors)) * 100
cat(sprintf("  -> Grey (Unassigned): %d genes (%.1f%%)\n", grey_size, grey_prop))

if(grey_prop > 30) {
  cat("     ! Note: High proportion of unassigned genes (>30%). \n")
  cat("             Consider lowering 'minModuleSize' or 'deepSplit' in next run if this is unexpected.\n")
}

# C. 检查空模块 (虽然 mergeCloseModules 也会处理，但多加一道防线)
empty_modules <- names(which(module_sizes == 0))
if(length(empty_modules) > 0) {
  cat(sprintf("  ! Warning: %d empty modules detected. (This is rare after merging)\n", length(empty_modules)))
}

# #### 4.6 导出结果 (增强版: 含 kME 值) ####
# ==============================================================================
# 计算 kME (基因与模块特征向量的相关性)
cat("  -> Calculating kME (Module Membership)...\n")
datME <- MEs
datKME <- signedKME(datExpr, datME, outputColumnName = "kME_")

# 构建主表
module_assignment <- data.frame(
  Gene_ID = colnames(datExpr),
  Module_Color = moduleColors,
  Module_Label = moduleLabels,
  stringsAsFactors = FALSE
)

# 合并 kME 数据 (只保留该基因所属模块的kME值，或者保留全部)
# 这里我们将所有模块的kME值都合并进去，方便筛选
module_assignment <- cbind(module_assignment, datKME)

# 如果有原始基因名映射，合并进去
if(exists("gene_name_mapping")) {
  merged_map <- gene_name_mapping[match(module_assignment$Gene_ID, gene_name_mapping$Valid_Rowname), ]
  module_assignment$Original_ID <- merged_map$Original_ID
  
  # 重新排列列顺序：把 Original_ID 放最前
  # 获取除了 Original_ID 之外的所有列名
  other_cols <- setdiff(colnames(module_assignment), "Original_ID")
  module_assignment <- module_assignment[, c("Original_ID", other_cols)]
}

# 保存
write.csv(module_assignment, "03_Modules/Module_Assignment_With_kME.csv", row.names = FALSE)
cat("  -> Module assignment (with kME) saved: 03_Modules/Module_Assignment_With_kME.csv\n")

# 保存 RData
geneTree <- net$dendrograms[[1]]
MEs <- orderMEs(MEs)
save(net, moduleLabels, moduleColors, MEs, geneTree, module_assignment, 
     file = "03_Modules/Network_Results.RData")

# 4.7 可视化
# ------------------------------------------------------------------
pdf("03_Modules/Gene_Dendrogram_Final.pdf", width = 12, height = 6)
plotDendroAndColors(geneTree, moduleColors, "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste("Gene Dendrogram (Final:", length(unique(moduleColors))-1, "modules)"))
dev.off()

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\n[Step 4] Completed in %.1f minutes.\n", total_time))

#### 5. 模块-性状关联分析 (可视化) ####
# ==============================================================================
cat("\n[Step 5] Module-Trait Correlation Heatmap...\n")

MEs0 <- orderMEs(MEs)
moduleTraitCor <- cor(MEs0, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 准备 ComplexHeatmap 注释
mod_names_num <- gsub("ME", "", rownames(moduleTraitCor))
mod_real_cols <- labels2colors(as.numeric(mod_names_num)) # 数字转回颜色
mod_sizes <- table(moduleColors)
size_vec <- as.numeric(mod_sizes[mod_real_cols])

ha_row <- rowAnnotation(
  Module = mod_real_cols,
  Size = anno_barplot(size_vec, width = unit(1.5, "cm"), gp = gpar(fill = "grey80")),
  col = list(Module = setNames(mod_real_cols, mod_real_cols)),
  show_legend = FALSE
)

# === 输出位置 3: 模块-性状关联热图 ===
pdf("04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf", width = 10, height = 8)
ht <- Heatmap(moduleTraitCor,
              name = "Cor",
              col = modern_pal$heat_map,
              right_annotation = ha_row,
              cluster_rows = TRUE, cluster_columns = TRUE,
              column_title = "Group / Sample Relationships",
              row_names_side = "left",
              cell_fun = function(j, i, x, y, width, height, fill) {
                p <- moduleTraitPvalue[i, j]
                if(!is.na(p) && p < 0.05) {
                  grid.text(ifelse(p<0.001, "***", "*"), x, y)
                }
              })
draw(ht)
dev.off()
cat("  -> Saved: 04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf\n")

#### 6. 核心基因深度分析 (特定模块: Magenta/Darkgreen) ####
# ==============================================================================
cat("\n[Step 6] Hub Gene Analysis (Network & Heatmap)...\n")

# 定义绘图函数
plot_module_analysis <- function(target_color, output_dir) {
  # 检查模块是否存在
  if (!(target_color %in% moduleColors)) {
    cat(sprintf("  ! Warning: Module '%s' not found (likely merged or filtered).\n", target_color))
    return(NULL)
  }
  
  # 创建目录
  mod_out <- file.path(output_dir, target_color)
  if(!dir.exists(mod_out)) dir.create(mod_out, recursive=TRUE)
  
  # 1. 提取基因
  mod_genes <- colnames(datExpr)[moduleColors == target_color]
  # 取前30个hub基因 (基于kME)
  datME <- MEs0[, paste0("ME", match(target_color, labels2colors(0:(ncol(MEs0)-1))))]
  kME <- abs(cor(datExpr[, mod_genes], datME, use="p"))
  top_genes <- rownames(kME)[order(kME[,1], decreasing=TRUE)[1:min(30, length(mod_genes))]]
  
  # 2. === 输出位置 4: 表达热图 ===
  expr_mat <- t(scale(datExpr[, top_genes]))
  # 简单注释
  anno_col <- data.frame(Group = meta_data$Group_Short_Name)
  rownames(anno_col) <- rownames(datExpr)
  
  pdf(file.path(mod_out, paste0("Heatmap_", target_color, ".pdf")), width=10, height=8)
  pheatmap(expr_mat, 
           annotation_col = anno_col, 
           show_colnames = FALSE,
           main = paste("Top 30 Hub Genes Expression - Module", target_color))
  dev.off()
  
  # 3. === 输出位置 5: 互作网络图 ===
  # 计算简单的相关性网络
  adj_mat <- abs(cor(datExpr[, top_genes]))
  adj_mat[adj_mat < 0.3] <- 0 # 阈值过滤
  diag(adj_mat) <- 0
  
  g <- graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=TRUE)
  
  if(vcount(g) > 0) {
    pdf(file.path(mod_out, paste0("Network_", target_color, ".pdf")), width=8, height=8)
    plot(g, layout = layout_with_fr(g),
         vertex.label = V(g)$name, vertex.label.cex = 0.7,
         vertex.color = target_color, vertex.frame.color = NA,
         vertex.size = 15, edge.width = E(g)$weight * 5,
         main = paste("Interaction Network - Module", target_color))
    dev.off()
  }
  
  cat(sprintf("  -> Analyzed Module: %s\n", target_color))
}

# 6.1 分析指定模块 (magenta, darkgreen)
# 注意：由于重新计算了，颜色分配可能发生变化。如果找不到，脚本会自动跳过。
plot_module_analysis("magenta", "05_Hub_Genes")
plot_module_analysis("darkgreen", "05_Hub_Genes")

# 6.2 补充分析：如果指定模块不存在，分析关联性最强的前2个模块
top_cor_mods <- rownames(moduleTraitCor)[order(apply(moduleTraitCor, 1, max), decreasing=TRUE)[1:2]]
top_colors <- labels2colors(as.numeric(gsub("ME", "", top_cor_mods)))

for(col in top_colors) {
  if(!col %in% c("magenta", "darkgreen")) { # 避免重复
    plot_module_analysis(col, "05_Hub_Genes")
  }
}

#### 7. KEGG 外部引导与 Summary 输出 ####
# ==============================================================================
cat("\n[Step 7] Generating Materials & Reports...\n")

# 7.1 === 输出位置 6: KEGG 引导文件 ===
# 导出每个模块的基因列表
for(mod in unique(moduleColors)) {
  genes <- colnames(datExpr)[moduleColors == mod]
  write.table(genes, file.path("06_Enrichment_Materials", paste0("Genes_", mod, ".txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

guide_text <- c(
  "【KEGG Enrichment Guide】",
  "Since R packages often lack optimal support for Medicago sativa (Alfalfa),",
  "please use the gene lists in this folder with the following tools:",
  "",
  "1. KOBAS (http://kobas.cbi.pku.edu.cn/)",
  "   - Upload 'Genes_[color].txt'",
  "   - Species: Medicago truncatula (closely related model)",
  "",
  "2. ShinyGO (http://bioinformatics.sdstate.edu/go/)",
  "   - Paste the gene IDs",
  "   - Select Medicago truncatula"
)
writeLines(guide_text, "06_Enrichment_Materials/READ_ME_KEGG.txt")

# 7.2 === 输出位置 7: 综合 Summary 报告 ===
summary_txt <- c(
  "===========================================================",
  "             WGCNA Analysis Summary Report (v7.0)",
  "===========================================================",
  paste("Date:", Sys.time()),
  "",
  "[1. Data Overview]",
  paste("   Total Genes Processed:", ncol(datExpr)),
  paste("   Samples Analyzed:     ", nrow(datExpr)),
  paste("   Filtering Strategy:    TPM>1 & Top 25,000 Variance (Optimized for 15GB RAM)"),
  "",
  "[2. Network Construction]",
  paste("   Soft Power:", sft_power),
  paste("   Modules Detected:", length(unique(moduleColors))),
  "   (See 03_Modules/ for dendrograms)",
  "",
  "[3. Key Visualizations Generated]",
  "   - Topology:    02_Network_Topology/Soft_Threshold_Plots.pdf",
  "   - Correlation: 04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf",
  "   - Hub Genes:   05_Hub_Genes/ (Check subfolders for specific modules)",
  "",
  "[4. Next Steps]",
  "   - Use gene lists in '06_Enrichment_Materials' for KEGG/GO analysis.",
  "   - Correlate module eigengenes with specific phenotype data.",
  "==========================================================="
)
writeLines(summary_txt, "08_Summary_Reports/Analysis_Summary.txt")

cat("\n======================================================\n")
cat("   Analysis Complete! All requested files restored.\n")
cat("======================================================\n")
