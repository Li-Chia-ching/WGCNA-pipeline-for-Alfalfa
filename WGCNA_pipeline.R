# ==============================================================================
# Medicago sativa WGCNA 完整分析流程 v7.1-Cleaned
# 硬件适配: AMD Ryzen 7 PRO / 15GB RAM
# 修复内容: 去除冗余绘图代码，优化内存逻辑，增强 Z-score 离群点检测流程
# ==============================================================================

#### 1. 环境与配置 ####
# ==============================================================================
cat("\n[Step 1] Environment Initialization...\n")
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 1.1 加载必要包
pkgs <- c("WGCNA", "dplyr", "stringr", "ggplot2", "pheatmap", "igraph", 
          "ComplexHeatmap", "circlize", "RColorBrewer", "viridis")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# 1.2 并行设置与随机种子
n_cores <- max(1, parallel::detectCores() - 2)
enableWGCNAThreads(nThreads = n_cores)
set.seed(12345) # 保证结果可复现

# 1.3 创建完整目录结构
dirs <- c("01_Data_Clean", "02_Network_Topology", "03_Modules", 
          "04_Module_Trait_Correlations", "05_Hub_Genes", 
          "05_Hub_Genes/magenta", "05_Hub_Genes/darkgreen", 
          "06_Enrichment_Materials", "08_Summary_Reports")
for(d in dirs) if(!dir.exists(d)) dir.create(d, recursive=TRUE)

# 1.4 现代化配色
modern_pal <- list(
  net_node = "#2C3E50", net_edge = "grey80",
  heat_map = colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
)

#### 2. 数据清洗 (针对 15GB RAM 优化 - 增强版) ####
# ==============================================================================
cat("\n[Step 2] Data Processing (Robust Filtering & Mapping)...\n")

# 2.1 读取数据检查
if (!exists("Expression_with_annotation")) stop("Missing: Expression_with_annotation")
if (!exists("Transcriptome_grouping_summary")) stop("Missing: Transcriptome_grouping_summary")

raw_expr <- as.data.frame(Expression_with_annotation)
meta_data <- as.data.frame(Transcriptome_grouping_summary)

#### 2.2 格式化表达矩阵 (混合列名自动识别) ####
# ------------------------------------------------------------------------------
cat("\n[Step 2.2] Formatting Expression Matrix (Extracting FPKM/TPM)...\n")

# 1. 提取基因ID
gene_ids <- as.character(raw_expr[, 1]) 

# 2. 识别 FPKM 或 TPM 列
all_cols <- colnames(raw_expr)
fpkm_cols_idx <- grep(":(fpkm|tpm)$", all_cols, ignore.case = TRUE)

if(length(fpkm_cols_idx) == 0) {
  cat("错误：无法识别表达量列。当前列名示例：\n")
  print(head(all_cols, 10))
  stop("无法找到以 ':fpkm' 或 ':tpm' 结尾的列！")
}

cat(sprintf("  -> 识别到 %d 个 FPKM/TPM 数据列，正在提取...\n", length(fpkm_cols_idx)))

# 3. 提取并清洗列名
datExpr0 <- raw_expr[, fpkm_cols_idx, drop = FALSE]
colnames(datExpr0) <- gsub(":(fpkm|tpm)$", "", colnames(datExpr0), ignore.case = TRUE)

# 4. 设置行名 (基因ID去重)
valid_gene_names <- make.names(gene_ids, unique = TRUE)
rownames(datExpr0) <- valid_gene_names

# 5. 保存基因映射表 (重要：用于后续找回原始ID)
gene_name_mapping <- data.frame(
  Original_ID = gene_ids,
  Valid_Rowname = valid_gene_names,
  stringsAsFactors = FALSE
)
write.csv(gene_name_mapping, "01_Data_Clean/Gene_Name_Mapping_Full.csv", row.names = FALSE)

#### 2.3 样本对齐 ####
# ------------------------------------------------------------------------------
cat("\n[Step 2.3] Sample Alignment...\n")

# 1. 获取分组表中的 Sample 列
sample_col_idx <- grep("^Sample$", colnames(meta_data), ignore.case = TRUE)[1]
if(is.na(sample_col_idx)) stop("Error: Grouping summary missing 'Sample' column.")

meta_samples <- trimws(gsub('^"|"$', '', as.character(meta_data[[sample_col_idx]])))
expr_samples <- colnames(datExpr0)

# 2. 匹配与对齐
common_samples <- intersect(expr_samples, meta_samples)

if(length(common_samples) == 0) {
  stop("无法匹配样本！请检查表达矩阵列名与分组表Sample列是否一致。")
}

datExpr0 <- datExpr0[, common_samples, drop = FALSE]
meta_data$SampleID_Clean <- meta_samples 
meta_data <- meta_data[match(common_samples, meta_data$SampleID_Clean), ]

cat(sprintf("  -> 成功对齐样本数: %d\n", ncol(datExpr0)))

#### 2.4 关键过滤 (适配 15GB RAM) ####
# ------------------------------------------------------------------------------
cat("  -> Filtering genes to fit 15GB RAM limit...\n")

# A. 丰度过滤 (>20% 样本中 >= 1)
keep_abd <- rowSums(datExpr0 >= 1) >= (0.2 * ncol(datExpr0))
datExpr_abd <- datExpr0[keep_abd, , drop = FALSE]
cat(sprintf("     Abundance filter: %d genes removed.\n", nrow(datExpr0) - nrow(datExpr_abd)))

# B. 方差过滤 (保留 Top 25,000)
target_n <- 25000
if(nrow(datExpr_abd) > target_n){
  mads <- apply(datExpr_abd, 1, mad)
  keep_var_idx <- order(mads, decreasing = TRUE)[1:target_n]
  datExpr_final <- datExpr_abd[keep_var_idx, , drop = FALSE]
  cat(sprintf("     Variance filter: Kept top %d genes (out of %d).\n", target_n, nrow(datExpr_abd)))
} else {
  datExpr_final <- datExpr_abd
  cat("     Variance filter: Gene count below 25,000. Kept all.\n")
}

# 保存过滤后的映射表
final_genes <- rownames(datExpr_final)
final_mapping <- gene_name_mapping[gene_name_mapping$Valid_Rowname %in% final_genes, ]
write.csv(final_mapping, "01_Data_Clean/Gene_Name_Mapping_Filtered.csv", row.names = FALSE)

#### 2.5 转置并标准化 ####
# ------------------------------------------------------------------------------
datExpr <- as.data.frame(t(log2(datExpr_final + 1)))

#### 2.6 构建性状矩阵 ####
# ------------------------------------------------------------------------------
datTraits <- data.frame(row.names = rownames(datExpr))
groups <- unique(meta_data$Group_Short_Name)
for(g in groups) {
  safe_g <- gsub("-", "_", g)
  datTraits[[safe_g]] <- ifelse(meta_data$Group_Short_Name == g, 1, 0)
}

#### 2.7 样本聚类与离群值检测 (整合版) ####
# ------------------------------------------------------------------------------
cat("\n[Step 2.7] Sample QC: Dendrogram & Z-score Analysis...\n")

# 1. 计算距离与 Z-score
sample_dist_matrix <- 1 - cor(t(datExpr), use = "p")
sampleTree <- hclust(as.dist(sample_dist_matrix), method = "average")

mean_dist_per_sample <- colMeans(sample_dist_matrix)
z_scores <- scale(mean_dist_per_sample)
outlier_threshold <- 2.5 
outliers <- rownames(z_scores)[which(z_scores > outlier_threshold)]
traitColors <- numbers2colors(datTraits, signed = FALSE)

# 2. 主图：聚类树 + 性状 + 标记
pdf("01_Data_Clean/Sample_Dendrogram_Trait_Heatmap_QC.pdf", width = 14, height = 10)
par(mar = c(1, 4, 3, 1))
plotDendroAndColors(sampleTree, traitColors, groupLabels = colnames(datTraits),
                    main = "Sample Dendrogram and Trait Heatmap (With QC)",
                    cex.dendroLabels = 0.6)

# 若有离群点，添加星号标记
if(length(outliers) > 0) {
  labels_in_plot_order <- sampleTree$labels[sampleTree$order]
  outlier_x_coords <- match(outliers, labels_in_plot_order)
  
  # 简单添加红色星号
  y_mark <- max(sampleTree$height) * 0.05
  text(x = outlier_x_coords, y = rep(y_mark, length(outliers)), 
       labels = "*", col = "red", cex = 2, font = 2)
  mtext(paste("Outliers detected:", length(outliers)), side=3, line=-1, col="red")
}
dev.off()

# 3. 辅助图：Z-score 分布与散点图 (如果存在离群点)
if(length(outliers) > 0) {
  cat("\n  ! WARNING: Potential outliers detected (Z > 2.5):\n")
  cat(paste("     -", outliers), sep = "\n")
  
  pdf("01_Data_Clean/Outlier_Validation_Plot.pdf", width = 12, height = 6)
  par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
  
  # 直方图
  hist(z_scores, breaks = 30, col = "skyblue", border = "white",
       main = "Z-score Distribution", xlab = "Z-score")
  abline(v = outlier_threshold, col = "red", lty = 2)
  
  # 散点图
  plot_colors <- ifelse(z_scores > outlier_threshold, "red", "gray50")
  plot(mean_dist_per_sample, z_scores, pch = 19, col = plot_colors,
       main = "Distance vs Z-score", xlab = "Mean Distance", ylab = "Z-score")
  abline(h = outlier_threshold, col = "red", lty = 2)
  text(mean_dist_per_sample[outliers], z_scores[outliers], 
       labels = outliers, pos = 4, col = "red", cex = 0.8)
  
  dev.off()
  
  # 保存 CSV
  outlier_info <- data.frame(Sample = outliers, Z_Score = z_scores[outliers,])
  write.csv(outlier_info, "01_Data_Clean/Potential_Outliers.csv", row.names = FALSE)
} else {
  cat("  -> QC Passed: No significant outliers detected.\n")
}

#### 2.8 保存 Checkpoint ####
# ------------------------------------------------------------------------------
save(datExpr, datTraits, meta_data, gene_name_mapping, 
     file = "01_Data_Clean/Input_Data.RData")
cat("  -> [Checkpoint] Processed data saved to: 01_Data_Clean/Input_Data.RData\n")

#### 3. 网络拓扑分析 (软阈值) ####
# ==============================================================================
cat("\n[Step 3] Soft Thresholding...\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
sft_power <- sft$powerEstimate
if(is.na(sft_power)) sft_power <- 9 # 默认值

pdf("02_Network_Topology/Soft_Threshold_Plots.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

#### 4. 模块识别 (15GB RAM 安全版) ####
# ==============================================================================
cat("\n[Step 4] Module Construction...\n")
start_time <- Sys.time()

# 4.1 再次确认基因数量 (防止爆内存)
# 这一步是双重保险，如果 2.4 步已经做了，这里不会有动作
target_gene_size <- 20000 
if (ncol(datExpr) > target_gene_size) {
  cat(sprintf("  -> Optimizing input: Focusing on top %d most variable genes...\n", target_gene_size))
  gene_mads <- apply(datExpr, 2, mad)
  keep_idx <- order(gene_mads, decreasing = TRUE)[1:target_gene_size]
  datExpr <- datExpr[, keep_idx]
}

# 4.2 构建网络
net <- blockwiseModules(
  datExpr,
  power = sft_power,
  maxBlockSize = 25000,     # 一次性计算
  TOMType = "unsigned",
  minModuleSize = 120,      
  mergeCutHeight = 0.25,    
  deepSplit = 2,            
  numericLabels = TRUE,
  saveTOMs = FALSE,         # 禁用 TOM 保存以节省时间和磁盘
  verbose = 2,
  nThreads = n_cores
)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# 4.3 模块合并 (如果模块过多)
n_modules <- length(unique(moduleColors)) - 1
if (n_modules > 15) {
  cat(sprintf("  -> Merging close modules (Current: %d)...\n", n_modules))
  merge_res <- mergeCloseModules(datExpr, moduleColors, cutHeight = 0.30, verbose = 0)
  moduleColors <- merge_res$colors
  MEs <- merge_res$newMEs
  moduleLabels <- match(moduleColors, c("grey", standardColors(200))) - 1
  if(any(is.na(moduleLabels))) moduleLabels <- as.numeric(factor(moduleColors)) - 1
  net$colors <- moduleLabels
  net$MEs <- MEs
}

# 4.4 导出结果 (含 kME)
cat("  -> Calculating kME and saving results...\n")
datKME <- signedKME(datExpr, MEs, outputColumnName = "kME_")
module_assignment <- data.frame(
  Gene_ID = colnames(datExpr),
  Module_Color = moduleColors,
  Module_Label = moduleLabels,
  stringsAsFactors = FALSE
)
module_assignment <- cbind(module_assignment, datKME)

# 合并原始ID
if(exists("gene_name_mapping")) {
  merged_map <- gene_name_mapping[match(module_assignment$Gene_ID, gene_name_mapping$Valid_Rowname), ]
  module_assignment$Original_ID <- merged_map$Original_ID
  module_assignment <- module_assignment[, c("Original_ID", setdiff(colnames(module_assignment), "Original_ID"))]
}
write.csv(module_assignment, "03_Modules/Module_Assignment_With_kME.csv", row.names = FALSE)

# 保存对象
geneTree <- net$dendrograms[[1]]
save(net, moduleLabels, moduleColors, MEs, geneTree, module_assignment, 
     file = "03_Modules/Network_Results.RData")

# 4.5 可视化
pdf("03_Modules/Gene_Dendrogram_Final.pdf", width = 12, height = 6)
plotDendroAndColors(geneTree, moduleColors, "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram")
dev.off()

cat(sprintf("\n[Step 4] Completed in %.1f minutes.\n", difftime(Sys.time(), start_time, units = "mins")))

#### 5. 模块-性状关联分析 ####
# ==============================================================================
cat("\n[Step 5] Module-Trait Correlation Heatmap...\n")

MEs0 <- orderMEs(MEs)
moduleTraitCor <- cor(MEs0, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 注释条
mod_real_cols <- labels2colors(as.numeric(gsub("ME", "", rownames(moduleTraitCor))))
mod_sizes <- table(moduleColors)
ha_row <- rowAnnotation(
  Module = mod_real_cols,
  Size = anno_barplot(as.numeric(mod_sizes[mod_real_cols]), width = unit(1.5, "cm"), gp = gpar(fill = "grey80")),
  col = list(Module = setNames(mod_real_cols, mod_real_cols)),
  show_legend = FALSE
)

pdf("04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf", width = 10, height = 8)
ht <- Heatmap(moduleTraitCor, name = "Cor", col = modern_pal$heat_map,
              right_annotation = ha_row, cluster_rows = TRUE, cluster_columns = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                p <- moduleTraitPvalue[i, j]
                if(!is.na(p) && p < 0.05) grid.text(ifelse(p<0.001, "***", "*"), x, y)
              })
draw(ht)
dev.off()

#### 6. 核心基因深度分析 (通用函数) ####
# ==============================================================================
cat("\n[Step 6] Hub Gene Analysis...\n")

plot_module_analysis <- function(target_color, output_dir) {
  if (!(target_color %in% moduleColors)) return(NULL)
  
  mod_out <- file.path(output_dir, target_color)
  if(!dir.exists(mod_out)) dir.create(mod_out, recursive=TRUE)
  
  # 提取 Top 30 Hub 基因
  mod_genes <- colnames(datExpr)[moduleColors == target_color]
  datME_curr <- MEs0[, paste0("ME", match(target_color, labels2colors(0:(ncol(MEs0)-1))))]
  kME_curr <- abs(cor(datExpr[, mod_genes], datME_curr, use="p"))
  top_genes <- rownames(kME_curr)[order(kME_curr[,1], decreasing=TRUE)[1:min(30, length(mod_genes))]]
  
  # 1. 热图
  expr_mat <- t(scale(datExpr[, top_genes]))
  anno_col <- data.frame(Group = meta_data$Group_Short_Name)
  rownames(anno_col) <- rownames(datExpr)
  
  pdf(file.path(mod_out, paste0("Heatmap_", target_color, ".pdf")), width=10, height=8)
  pheatmap(expr_mat, annotation_col = anno_col, show_colnames = FALSE,
           main = paste("Top 30 Hub Genes - Module", target_color))
  dev.off()
  
  # 2. 网络图
  adj_mat <- abs(cor(datExpr[, top_genes]))
  adj_mat[adj_mat < 0.3] <- 0
  diag(adj_mat) <- 0
  g <- graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=TRUE)
  
  if(vcount(g) > 0) {
    pdf(file.path(mod_out, paste0("Network_", target_color, ".pdf")), width=8, height=8)
    plot(g, layout = layout_with_fr(g), vertex.label = V(g)$name, 
         vertex.color = target_color, vertex.size = 15, edge.width = E(g)$weight * 5,
         main = paste("Interaction Network - Module", target_color))
    dev.off()
  }
  cat(sprintf("  -> Analyzed Module: %s\n", target_color))
}

# 分析特定模块 + Top 2 关联模块
target_list <- c("magenta", "darkgreen")
top_cor_mods <- rownames(moduleTraitCor)[order(apply(moduleTraitCor, 1, max), decreasing=TRUE)[1:2]]
target_list <- unique(c(target_list, labels2colors(as.numeric(gsub("ME", "", top_cor_mods)))))

for(col in target_list) plot_module_analysis(col, "05_Hub_Genes")

#### 7. KEGG 引导与报告 ####
# ==============================================================================
cat("\n[Step 7] Generating Reports...\n")

# 导出基因列表
for(mod in unique(moduleColors)) {
  genes <- colnames(datExpr)[moduleColors == mod]
  write.table(genes, file.path("06_Enrichment_Materials", paste0("Genes_", mod, ".txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Summary
summary_txt <- c(
  "=== WGCNA Analysis Summary ===",
  paste("Date:", Sys.time()),
  paste("Genes Processed:", ncol(datExpr)),
  paste("Modules Detected:", length(unique(moduleColors))),
  "Soft Power Used: 9 (or calculated)",
  "Results saved in 03_Modules and 05_Hub_Genes."
)
writeLines(summary_txt, "08_Summary_Reports/Analysis_Summary.txt")

cat("\n======================================================\n")
cat("   Analysis Complete! (Duplicated code removed)\n")
cat("======================================================\n")
