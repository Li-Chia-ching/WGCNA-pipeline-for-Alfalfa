# ==============================================================================
# Medicago sativa WGCNA 完整分析流程 v7.2-Final-Data (可视化+数据双输出版)
# 硬件适配: AMD Ryzen 7 PRO / 15GB RAM
# 功能: 核心计算 + PDF可视化 + CSV源数据导出 + 增强型质控
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
          "06_Enrichment_Materials", "08_Summary_Reports")
for(d in dirs) if(!dir.exists(d)) dir.create(d, recursive=TRUE)

# 配色设置
modern_pal <- list(
  heat_map = colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
)

#### 2. 数据清洗 (完整逻辑+数据导出) ####
# ==============================================================================
cat("\n[Step 2] Data Processing & Export...\n")

# 2.1 读取数据 (请确保变量已在工作空间中)
if (!exists("Expression_with_annotation")) stop("Missing: Expression_with_annotation")
if (!exists("Transcriptome_grouping_summary")) stop("Missing: Transcriptome_grouping_summary")

raw_expr <- as.data.frame(Expression_with_annotation)
meta_data <- as.data.frame(Transcriptome_grouping_summary)

# 2.2 格式化表达矩阵
# ------------------------------------------------------------------
gene_ids <- as.character(raw_expr[, 1])
fpkm_cols_idx <- grep(":(fpkm|tpm)$", colnames(raw_expr), ignore.case = TRUE)
if(length(fpkm_cols_idx) == 0) stop("无法找到 :fpkm 或 :tpm 列")

datExpr0 <- raw_expr[, fpkm_cols_idx, drop = FALSE]
colnames(datExpr0) <- gsub(":(fpkm|tpm)$", "", colnames(datExpr0), ignore.case = TRUE)
valid_gene_names <- make.names(gene_ids, unique = TRUE)
rownames(datExpr0) <- valid_gene_names

# 保存基因映射表
gene_name_mapping <- data.frame(Original_ID = gene_ids, Valid_Rowname = valid_gene_names)
write.csv(gene_name_mapping, "01_Data_Clean/Gene_Name_Mapping_Full.csv", row.names = FALSE)

# 2.3 样本对齐
# ------------------------------------------------------------------
sample_col_idx <- grep("^Sample$", colnames(meta_data), ignore.case = TRUE)[1]
meta_samples <- trimws(gsub('^"|"$', '', as.character(meta_data[[sample_col_idx]])))
common_samples <- intersect(colnames(datExpr0), meta_samples)
if(length(common_samples) == 0) stop("样本无法匹配")

datExpr0 <- datExpr0[, common_samples, drop = FALSE]
meta_data$SampleID_Clean <- meta_samples
meta_data <- meta_data[match(common_samples, meta_data$SampleID_Clean), ]

# 2.4 关键过滤 (15GB RAM 优化)
# ------------------------------------------------------------------
cat("  -> Filtering genes...\n")
keep_abd <- rowSums(datExpr0 >= 1) >= (0.2 * ncol(datExpr0))
datExpr_abd <- datExpr0[keep_abd, , drop = FALSE]

target_n <- 25000
if(nrow(datExpr_abd) > target_n){
  mads <- apply(datExpr_abd, 1, mad)
  datExpr_final <- datExpr_abd[order(mads, decreasing = TRUE)[1:target_n], , drop = FALSE]
} else {
  datExpr_final <- datExpr_abd
}
datExpr <- as.data.frame(t(log2(datExpr_final + 1)))

# 保存过滤后的映射表 (用于Step 4导出)
final_mapping <- gene_name_mapping[gene_name_mapping$Valid_Rowname %in% rownames(datExpr_final), ]
write.csv(final_mapping, "01_Data_Clean/Gene_Name_Mapping_Filtered.csv", row.names = FALSE)

# 2.6 构建性状矩阵
# ------------------------------------------------------------------
datTraits <- data.frame(row.names = rownames(datExpr))
for(g in unique(meta_data$Group_Short_Name)) {
  datTraits[[gsub("-", "_", g)]] <- ifelse(meta_data$Group_Short_Name == g, 1, 0)
}

# 2.7 样本聚类、QC与数据导出 (核心修改)
# ------------------------------------------------------------------
cat("\n[Step 2.7] Sample QC: Plotting & Data Export...\n")

# A. 计算距离并保存
sample_dist_matrix <- 1 - cor(t(datExpr), use = "p")
write.csv(sample_dist_matrix, "01_Data_Clean/Sample_Distance_Matrix.csv", row.names = TRUE)

# B. Z-score 计算与保存
mean_dist <- colMeans(sample_dist_matrix)
z_scores <- scale(mean_dist)
outlier_thresh <- 2.5
outliers <- rownames(z_scores)[which(z_scores > outlier_thresh)]

sample_stats <- data.frame(
  Sample = rownames(sample_dist_matrix),
  Mean_Distance = mean_dist,
  Z_Score = z_scores[,1],
  Is_Outlier = ifelse(z_scores[,1] > outlier_thresh, "YES", "NO")
)
write.csv(sample_stats, "01_Data_Clean/Sample_QC_Statistics.csv", row.names = FALSE)
cat("  -> Saved QC Data: Distance Matrix & Statistics CSVs\n")

# C. 绘图 (保留可视化)
sampleTree <- hclust(as.dist(sample_dist_matrix), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)

pdf("01_Data_Clean/Sample_Dendrogram_Trait_Heatmap_QC.pdf", width = 14, height = 10)
par(mar = c(1, 4, 3, 1))
plotDendroAndColors(sampleTree, traitColors, groupLabels = colnames(datTraits),
                    main = "Sample Dendrogram and Trait Heatmap")
if(length(outliers) > 0) {
  # 简单标记离群点
  x_coords <- match(outliers, sampleTree$labels[sampleTree$order])
  text(x = x_coords, y = rep(max(sampleTree$height)*0.05, length(outliers)), 
       labels = "*", col = "red", cex = 2)
}
dev.off()

# 保存 Checkpoint
save(datExpr, datTraits, meta_data, final_mapping, file = "01_Data_Clean/Input_Data.RData")

#### 3. 网络拓扑分析 (双输出) ####
# ==============================================================================
cat("\n[Step 3] Soft Thresholding: Plotting & Data Export...\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
sft_power <- sft$powerEstimate
if(is.na(sft_power)) sft_power <- 9

# A. 保存数据
write.csv(sft$fitIndices, "02_Network_Topology/Soft_Threshold_Fit_Indices.csv", row.names = FALSE)

# B. 绘图
pdf("02_Network_Topology/Soft_Threshold_Plots.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.85, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

#### 4. 模块识别 (完整导出) ####
# ==============================================================================
cat("\n[Step 4] Module Construction & Export...\n")

# 4.1 再次内存安全检查
target_gene_size <- 20000 
if (ncol(datExpr) > target_gene_size) {
  cat("  -> Optimizing input for module construction...\n")
  gene_mads <- apply(datExpr, 2, mad)
  datExpr <- datExpr[, order(gene_mads, decreasing = TRUE)[1:target_gene_size]]
  # 更新映射表
  final_mapping <- final_mapping[final_mapping$Valid_Rowname %in% colnames(datExpr), ]
}

# 4.2 构建网络
net <- blockwiseModules(datExpr, power = sft_power, maxBlockSize = 25000,
                        TOMType = "unsigned", minModuleSize = 120,
                        mergeCutHeight = 0.25, deepSplit = 2, numericLabels = TRUE,
                        saveTOMs = FALSE, verbose = 1, nThreads = n_cores)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# 4.3 模块合并处理
if ((length(unique(moduleColors))-1) > 15) {
  merge_res <- mergeCloseModules(datExpr, moduleColors, cutHeight = 0.30, verbose = 0)
  moduleColors <- merge_res$colors
  MEs <- merge_res$newMEs
}
# 重新映射颜色
moduleLabels <- match(moduleColors, c("grey", standardColors(200))) - 1
if(any(is.na(moduleLabels))) moduleLabels <- as.numeric(factor(moduleColors)) - 1
net$colors <- moduleLabels; net$MEs <- MEs

# 4.4 可视化
pdf("03_Modules/Gene_Dendrogram_Final.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors, "Module Colors",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, main = "Gene Dendrogram")
dev.off()

# 4.5 导出全量数据
cat("  -> Exporting module data...\n")
datKME <- signedKME(datExpr, MEs, outputColumnName = "kME_")

# 主表
module_assignment <- data.frame(
  Gene_ID = colnames(datExpr),
  Module_Color = moduleColors,
  Module_Label = moduleLabels,
  stringsAsFactors = FALSE
)
module_assignment <- cbind(module_assignment, datKME)

# 合并原始ID
if(exists("final_mapping")) {
  module_assignment$Original_ID <- final_mapping$Original_ID[match(module_assignment$Gene_ID, final_mapping$Valid_Rowname)]
  module_assignment <- module_assignment[, c("Original_ID", setdiff(colnames(module_assignment), "Original_ID"))]
}

write.csv(module_assignment, "03_Modules/Module_Assignment_Full.csv", row.names = FALSE)
write.csv(MEs, "03_Modules/Module_Eigengenes.csv", row.names = TRUE)

module_size_stats <- as.data.frame(table(Module_Color = moduleColors))
write.csv(module_size_stats, "03_Modules/Module_Size_Statistics.csv", row.names = FALSE)

#### 5. 模块-性状关联分析 (双输出) ####
# ==============================================================================
cat("\n[Step 5] Module-Trait Analysis: Heatmap & Data Export...\n")

MEs0 <- orderMEs(MEs)
moduleTraitCor <- cor(MEs0, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# A. 保存数据
write.csv(moduleTraitCor, "04_Module_Trait_Correlations/Module_Trait_Correlation.csv", row.names = TRUE)
write.csv(moduleTraitPvalue, "04_Module_Trait_Correlations/Module_Trait_Pvalue.csv", row.names = TRUE)

# 创建综合表格
combined_res <- data.frame(Module = rownames(moduleTraitCor))
for(trait in colnames(moduleTraitCor)) {
  combined_res[[paste0(trait, "_Cor")]] <- moduleTraitCor[, trait]
  combined_res[[paste0(trait, "_Pval")]] <- moduleTraitPvalue[, trait]
}
write.csv(combined_res, "04_Module_Trait_Correlations/Module_Trait_Combined_Results.csv", row.names = FALSE)

# B. 绘图 (ComplexHeatmap)
mod_sizes <- table(moduleColors)
row_ha <- rowAnnotation(
  Size = anno_barplot(as.numeric(mod_sizes[labels2colors(as.numeric(gsub("ME", "", rownames(moduleTraitCor))))]), 
                      width = unit(1.5, "cm"), gp = gpar(fill = "grey80")),
  show_legend = FALSE
)

pdf("04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf", width = 10, height = 8)
ht <- Heatmap(moduleTraitCor, name = "Cor", col = modern_pal$heat_map,
              right_annotation = row_ha, cluster_rows = TRUE, cluster_columns = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                p <- moduleTraitPvalue[i, j]
                if(!is.na(p) && p < 0.05) grid.text(ifelse(p<0.001, "***", "*"), x, y)
              })
draw(ht)
dev.off()

#### 6. 核心基因深度分析 (双输出) ####
# ==============================================================================
cat("\n[Step 6] Hub Gene Analysis: Networks, Heatmaps & CSVs...\n")

plot_module_analysis <- function(target_color, output_dir) {
  if (!(target_color %in% moduleColors)) return(NULL)
  
  mod_out <- file.path(output_dir, target_color)
  if(!dir.exists(mod_out)) dir.create(mod_out, recursive=TRUE)
  
  # 1. 数据准备
  mod_genes <- colnames(datExpr)[moduleColors == target_color]
  datME_curr <- MEs0[, paste0("ME", match(target_color, labels2colors(0:(ncol(MEs0)-1))))]
  kME_curr <- abs(cor(datExpr[, mod_genes], datME_curr, use="p"))
  
  # 2. 导出完整 Hub 基因排名
  hub_ranking <- data.frame(Gene_ID = rownames(kME_curr), kME = kME_curr[,1])
  hub_ranking <- hub_ranking[order(hub_ranking$kME, decreasing = TRUE), ]
  write.csv(hub_ranking, file.path(mod_out, paste0("Hub_Gene_Ranking_", target_color, ".csv")), row.names = FALSE)
  
  # 3. 提取 Top 30
  top_n <- min(30, length(mod_genes))
  top_genes <- hub_ranking$Gene_ID[1:top_n]
  
  # 4. 表达热图 (Plot + CSV)
  expr_mat <- t(scale(datExpr[, top_genes]))
  write.csv(expr_mat, file.path(mod_out, paste0("Top_Hub_Genes_Expression_", target_color, ".csv")), row.names = TRUE)
  
  pdf(file.path(mod_out, paste0("Heatmap_", target_color, ".pdf")), width=10, height=8)
  pheatmap(expr_mat, annotation_col = data.frame(Group = meta_data$Group_Short_Name, row.names=rownames(datExpr)),
           show_colnames = FALSE, main = paste("Top 30 Hub Genes - Module", target_color))
  dev.off()
  
  # 5. 网络图 (Plot + CSV Edge List)
  if(length(top_genes) > 1) {
    adj_mat <- abs(cor(datExpr[, top_genes]))
    adj_mat[adj_mat < 0.3] <- 0; diag(adj_mat) <- 0
    
    # 导出边列表
    g <- graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=TRUE)
    edge_df <- as_data_frame(g, what="edges")
    if(nrow(edge_df) > 0) {
      write.csv(edge_df, file.path(mod_out, paste0("Gene_Interaction_Edges_", target_color, ".csv")), row.names = FALSE)
      
      # 绘图
      pdf(file.path(mod_out, paste0("Network_", target_color, ".pdf")), width=8, height=8)
      plot(g, layout = layout_with_fr(g), vertex.label = V(g)$name, 
           vertex.color = target_color, vertex.size = 15, edge.width = E(g)$weight * 5,
           main = paste("Interaction Network - Module", target_color))
      dev.off()
    }
  }
  cat(sprintf("  -> Analyzed %s (Data & Plots saved)\n", target_color))
}

# 执行分析
target_list <- unique(c("magenta", "darkgreen", 
                        labels2colors(as.numeric(gsub("ME", "", rownames(moduleTraitCor)[order(apply(moduleTraitCor, 1, max), decreasing=TRUE)[1:2]])))))
for(col in target_list) plot_module_analysis(col, "05_Hub_Genes")

#### 7. 导出验证与总结 ####
# ==============================================================================
cat("\n[Step 7] Generating Summary & Verification...\n")

# KEGG列表导出
for(mod in unique(moduleColors)) {
  write.table(colnames(datExpr)[moduleColors == mod], 
              file.path("06_Enrichment_Materials", paste0("Genes_", mod, ".txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# 验证文件存在性
check_files <- c("Sample_Distance_Matrix.csv", "Soft_Threshold_Fit_Indices.csv", 
                 "Module_Assignment_Full.csv", "Module_Trait_Correlation.csv")
missing <- c()
for(f in check_files) {
  # 简单的递归搜索检查
  found <- list.files(pattern = f, recursive = TRUE)
  if(length(found) == 0) missing <- c(missing, f)
}

if(length(missing) == 0) {
  cat("\nSUCCESS: All critical data files and plots have been generated.\n")
  cat("Output directories: 01_Data_Clean to 06_Enrichment_Materials\n")
} else {
  cat("\nWARNING: Some files seem missing:\n")
  print(missing)
}

cat("\n======================================================\n")
cat("   Analysis Complete - Data & Visualization Ready\n")
cat("======================================================\n")
