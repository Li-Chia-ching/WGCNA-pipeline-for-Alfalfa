# ==============================================================================
# Medicago sativa WGCNA v9.3-Visual-Pro (可视化增强版)
# 硬件适配: 16GB RAM / 146k Genes (Filtered) / 54 Samples
# 核心升级: 
#   1. [Step 5] 模块-性状热图增加实体颜色注释条 (不再只有文字)
#   2. [Step 6] 优化配色方案，提升图表美观度
#   3. [全局] 保持 v9.2 的逻辑修复和性能优化
# ==============================================================================

#### 1. 环境与配置 ####
# ==============================================================================
cat("\n[Step 1] Environment Initialization...\n")
options(stringsAsFactors = FALSE)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 加载包
pkgs <- c("WGCNA", "dplyr", "stringr", "ggplot2", "pheatmap", "igraph", 
          "ComplexHeatmap", "circlize", "RColorBrewer", "viridis")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# 并行与内存清理
n_cores <- min(parallel::detectCores() - 1, 4)
enableWGCNAThreads(nThreads = n_cores)
clean_mem <- function() { gc(verbose = FALSE) }

# 目录结构
dirs <- c("01_Data_Clean", "02_Network_Topology", "03_Modules", 
          "04_Module_Trait_Correlations", "05_Hub_Genes", 
          "06_Enrichment_Materials", "08_Summary_Reports")
for(d in dirs) if(!dir.exists(d)) dir.create(d, recursive=TRUE)

# 定义出版级配色方案
pub_pal <- list(
  heat = colorRamp2(c(-1, 0, 1), c("#2166AC", "#FFFFFF", "#B2182B")), # 经典的蓝-白-红
  group = brewer.pal(8, "Set2") # 柔和的分组配色
)

#### 2. 数据清洗 ####
# ==============================================================================
cat("\n[Step 2] Data Processing & Export...\n")

if (!exists("Expression_with_annotation")) stop("Missing: Expression_with_annotation")
if (!exists("Transcriptome_grouping_summary")) stop("Missing: Transcriptome_grouping_summary")

raw_expr <- as.data.frame(Expression_with_annotation)
meta_data <- as.data.frame(Transcriptome_grouping_summary)

# 2.2 格式化
gene_ids <- as.character(raw_expr[, 1])
fpkm_cols_idx <- grep(":(fpkm|tpm)$", colnames(raw_expr), ignore.case = TRUE)
if(length(fpkm_cols_idx) == 0) stop("列名错误: 未找到 :fpkm 或 :tpm")

datExpr0 <- raw_expr[, fpkm_cols_idx, drop = FALSE]
colnames(datExpr0) <- gsub(":(fpkm|tpm)$", "", colnames(datExpr0), ignore.case = TRUE)
rownames(datExpr0) <- make.names(gene_ids, unique = TRUE)

gene_name_mapping <- data.frame(Original_ID = gene_ids, Valid_Rowname = rownames(datExpr0))
write.csv(gene_name_mapping, "01_Data_Clean/Gene_Name_Mapping_Full.csv", row.names = FALSE)
rm(raw_expr); clean_mem()

# 2.3 样本对齐
sample_col_idx <- grep("^Sample$", colnames(meta_data), ignore.case = TRUE)[1]
meta_samples <- trimws(gsub('^"|"$', '', as.character(meta_data[[sample_col_idx]])))
common_samples <- intersect(colnames(datExpr0), meta_samples)

if(length(common_samples) == 0) stop("错误：样本ID无法匹配")

meta_data$SampleID_Clean <- meta_samples 
datExpr0 <- datExpr0[, common_samples, drop = FALSE]
meta_data <- meta_data[match(common_samples, meta_data$SampleID_Clean), ]
cat(sprintf("  -> Samples aligned: %d samples.\n", ncol(datExpr0)))

# 2.4 智能过滤 (Top 25000)
cat("\n  [Filtering Strategy]...\n")
keep_abd <- rowSums(datExpr0 >= 1) >= (0.1 * ncol(datExpr0))
datExpr_abd <- datExpr0[keep_abd, , drop = FALSE]
rm(datExpr0); clean_mem()

target_n <- 25000 
if(nrow(datExpr_abd) > target_n){
  cat(sprintf("  -> Reducing gene count to %d (Top Variance)...\n", target_n))
  mads <- apply(datExpr_abd, 1, mad)
  keep_idx <- order(mads, decreasing = TRUE)[1:target_n]
  datExpr_final <- datExpr_abd[keep_idx, , drop = FALSE]
} else {
  datExpr_final <- datExpr_abd
}
rm(datExpr_abd); clean_mem()

datExpr <- as.data.frame(t(log2(datExpr_final + 1)))
final_mapping <- gene_name_mapping[gene_name_mapping$Valid_Rowname %in% colnames(datExpr), ]
write.csv(final_mapping, "01_Data_Clean/Gene_Name_Mapping_Filtered.csv", row.names = FALSE)

# 2.5 构建性状矩阵
datTraits <- data.frame(row.names = rownames(datExpr))
groups <- unique(meta_data$Group_Short_Name)
for(g in groups) {
  safe_g <- gsub("[- ]", "_", g) 
  datTraits[[safe_g]] <- ifelse(meta_data$Group_Short_Name == g, 1, 0)
}

# 2.6 样本QC (优化图表)
sample_dist_matrix <- 1 - cor(t(datExpr), use = "p")
write.csv(sample_dist_matrix, "01_Data_Clean/Sample_Distance_Matrix.csv", row.names = TRUE)

sampleTree <- hclust(as.dist(sample_dist_matrix), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)

pdf("01_Data_Clean/Sample_Dendrogram_Trait_Heatmap_QC.pdf", width = 14, height = 10)
plotDendroAndColors(sampleTree, traitColors, groupLabels = colnames(datTraits),
                    main = "Sample Dendrogram & Trait Heatmap",
                    cex.dendroLabels = 0.6)
dev.off()

qc_stats <- data.frame(Sample = rownames(sample_dist_matrix), 
                       Z_Score = scale(colMeans(sample_dist_matrix)), 
                       Mean_Dist = colMeans(sample_dist_matrix))
write.csv(qc_stats, "01_Data_Clean/Sample_QC_Stats.csv", row.names = FALSE)
save(datExpr, datTraits, meta_data, final_mapping, file = "01_Data_Clean/Input_Data.RData")
clean_mem()

#### 3. 网络拓扑分析 ####
# ==============================================================================
cat("\n[Step 3] Soft Thresholding...\n")
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
sft_power <- sft$powerEstimate
if(is.na(sft_power)) sft_power <- 9

write.csv(sft$fitIndices, "02_Network_Topology/Soft_Threshold_Fit_Indices.csv", row.names = FALSE)

pdf("02_Network_Topology/Soft_Threshold_Plots.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold", ylab="Scale Free Fit R^2", type="n", main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.85, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold", ylab="Mean Connectivity", type="n", main="Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()
clean_mem()

#### 4. 模块识别 ####
# ==============================================================================
cat("\n[Step 4] Module Construction (Optimized)...\n")

net <- blockwiseModules(
  datExpr,
  power = sft_power,
  maxBlockSize = 9000,
  TOMType = "unsigned",
  minModuleSize = 100,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3,
  nThreads = n_cores
)
clean_mem()

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# 4.3 模块合并
n_modules_init <- length(unique(moduleColors)) - 1
cat(sprintf("\n  -> Initial module count: %d (excluding grey)\n", n_modules_init))

if (n_modules_init > 15) {
  cat("  -> Merging closely related modules (threshold = 0.30)...\n")
  merge_res <- mergeCloseModules(datExpr, moduleColors, cutHeight = 0.30, verbose = 0)
  moduleColors <- merge_res$colors
  MEs <- merge_res$newMEs
  
  moduleLabels <- match(moduleColors, c("grey", standardColors(200))) - 1
  if(any(is.na(moduleLabels))) moduleLabels <- as.numeric(factor(moduleColors)) - 1
  net$colors <- moduleLabels
  net$MEs <- MEs
  cat(sprintf("  -> After merging: %d modules\n", length(unique(moduleColors))-1))
}

datKME <- signedKME(datExpr, MEs, outputColumnName = "kME_")

# 4.4 导出
mod_df <- data.frame(Gene_ID = colnames(datExpr), Module = moduleColors)
mod_df <- cbind(mod_df, datKME)
mod_df$Original_ID <- final_mapping$Original_ID[match(mod_df$Gene_ID, final_mapping$Valid_Rowname)]
mod_df <- mod_df[, c("Original_ID", setdiff(colnames(mod_df), "Original_ID"))]

write.csv(mod_df, "03_Modules/Module_Assignment_Full.csv", row.names = FALSE)
write.csv(MEs, "03_Modules/Module_Eigengenes.csv", row.names = TRUE)
write.csv(as.data.frame(table(moduleColors)), "03_Modules/Module_Size_Statistics.csv", row.names = FALSE)

# 4.5 导出树状图
pdf("03_Modules/Gene_Dendrogram_Final_All_Blocks.pdf", width = 12, height = 6)
for(b in 1:length(net$dendrograms)){
  plotDendroAndColors(net$dendrograms[[b]], moduleColors[net$blockGenes[[b]]],
                      "Module Colors", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste("Gene Dendrogram - Block", b))
}
dev.off()

#### 5. 模块-性状关联 (视觉增强 - 最终修复版) ####
# ==============================================================================
cat("\n[Step 5] Module-Trait Analysis & Visualization...\n")
MEs0 <- orderMEs(MEs)
moduleTraitCor <- cor(MEs0, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 5.1 导出CSV数据（保持不变）
write.csv(moduleTraitCor, "04_Module_Trait_Correlations/Module_Trait_Correlation.csv", row.names = TRUE)
write.csv(moduleTraitPvalue, "04_Module_Trait_Correlations/Module_Trait_Pvalue.csv", row.names = TRUE)

combined_res <- data.frame(Module = rownames(moduleTraitCor))
for(trait in colnames(moduleTraitCor)) {
  combined_res[[paste0(trait, "_Cor")]] <- moduleTraitCor[, trait]
  combined_res[[paste0(trait, "_Pval")]] <- moduleTraitPvalue[, trait]
}
write.csv(combined_res, "04_Module_Trait_Correlations/Module_Trait_Combined_Results.csv", row.names = FALSE)

# 5.2 准备模块信息（关键：确保正确提取）
mod_names_long <- rownames(moduleTraitCor)  # "ME1", "ME2"等
mod_numbers <- as.numeric(gsub("ME", "", mod_names_long))
mod_real_colors <- labels2colors(mod_numbers)

# 5.3 创建现代化热图（修复注释条问题）
pdf("04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf", width = 14, height = 10)

# 颜色映射函数
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

# ==========================================================================
# 【关键修复】正确创建注释对象
# ==========================================================================

# A. 创建右侧注释：模块大小条形图
mod_sizes <- table(moduleColors)
# 确保顺序与热图行一致
mod_size_vec <- numeric(length(mod_numbers))
for(i in 1:length(mod_numbers)) {
  mod_color <- mod_real_colors[i]
  mod_size_vec[i] <- ifelse(mod_color %in% names(mod_sizes), 
                            mod_sizes[mod_color], 0)
}

ha_row <- rowAnnotation(
  Size = anno_barplot(mod_size_vec, 
                      width = unit(2, "cm"),
                      gp = gpar(fill = "#2C3E50", col = NA),
                      border = FALSE),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

# B. 创建顶部注释：性状分组（如果需要）
# 如果性状有分组信息，可以创建顶部注释
# 这里假设每个性状是独立的，所以不添加顶部注释

# ==========================================================================
# 创建热图对象
# ==========================================================================

ht <- Heatmap(moduleTraitCor,
              name = "Correlation\n(r)", # 图例名称更清晰
              col = col_fun,
              
              # --- 矩阵样式 ---
              rect_gp = gpar(col = "white", lwd = 0.5),
              na_col = "grey90",
              
              # --- 行设置：添加模块颜色标签 ---
              row_labels = paste0(mod_names_long, " (", mod_real_colors, ")"),
              row_names_gp = gpar(fontsize = 11, 
                                  col = ifelse(mod_real_colors == "white", 
                                               "black", "black")),
              
              # --- 行聚类设置 ---
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "complete",
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              
              # --- 列设置 ---
              column_names_gp = gpar(fontsize = 11),
              column_names_rot = 45,
              column_names_centered = FALSE,
              
              # --- 标题 ---
              column_title = "Module-Trait Relationships",
              column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "#2C3E50"),
              column_title_side = "top",
              
              # --- 【关键】添加注释对象 ---
              right_annotation = ha_row,
              
              # --- 图例设置 ---
              heatmap_legend_param = list(
                title = "Pearson r",
                title_gp = gpar(fontsize = 11, fontface = "bold"),
                labels_gp = gpar(fontsize = 10),
                legend_height = unit(3, "cm"),
                grid_width = unit(0.5, "cm"),
                at = c(-1, -0.5, 0, 0.5, 1),
                border = "black"
              ),
              
              # --- 单元格内添加显著性星号（优化版本）---
              cell_fun = function(j, i, x, y, width, height, fill) {
                pval <- moduleTraitPvalue[i, j]
                if(!is.na(pval) && pval < 0.05) {
                  star <- ifelse(pval < 0.001, "***", 
                                 ifelse(pval < 0.01, "**", "*"))
                  
                  # 智能选择星号颜色
                  bg_color <- col2rgb(fill)
                  brightness <- (bg_color[1]*0.299 + bg_color[2]*0.587 + bg_color[3]*0.114) / 255
                  star_col <- ifelse(brightness > 0.5, "black", "white")
                  
                  # 将星号放在单元格上半部分
                  grid.text(star, x, y - height*0.25, 
                            gp = gpar(fontsize = 9, fontface = "bold", col = star_col))
                }
              },
              
              # --- 其他优化 ---
              show_heatmap_legend = TRUE,
              use_raster = FALSE,
              border = TRUE,
              border_gp = gpar(col = "grey80", lwd = 0.5)
)

# 5.4 绘制并保存
draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 4, 2), "mm") # 上右下左边距
)

dev.off()
cat("  -> Saved enhanced heatmap: 04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf\n")

# 5.5 可选：生成PNG版本（用于快速查看）
png("04_Module_Trait_Correlations/Module_Trait_Heatmap.png", 
    width = 1400, height = 1000, res = 150)
draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 4, 2), "mm"))
dev.off()
cat("  -> Saved PNG version: 04_Module_Trait_Correlations/Module_Trait_Heatmap.png\n")

#### 6. Hub 基因分析 (视觉增强版 - 修复ME匹配) ####
# ==============================================================================
cat("\n[Step 6] Hub Gene Analysis (Visual Enhanced)...\n")

plot_module_analysis <- function(target_color, output_dir) {
  if (!(target_color %in% moduleColors)) {
    cat(sprintf("  Skipping %s: Module not found.\n", target_color))
    return(NULL)
  }
  
  mod_out <- file.path(output_dir, target_color)
  if(!dir.exists(mod_out)) dir.create(mod_out, recursive=TRUE)
  
  # 数据准备 - 修复ME列名匹配
  mod_genes <- colnames(datExpr)[moduleColors == target_color]
  
  # 关键修复：找到对应模块的ME列名
  # 方法1：遍历MEs0的列名，提取颜色部分进行匹配
  ME_col_name <- NULL
  for(col in colnames(MEs0)) {
    if(grepl("^ME", col)) {
      mod_num <- as.numeric(gsub("ME", "", col))
      mod_col_from_num <- labels2colors(mod_num)
      if(mod_col_from_num == target_color) {
        ME_col_name <- col
        break
      }
    }
  }
  
  if(is.null(ME_col_name)) {
    cat(sprintf("  Warning: Cannot find ME column for module %s\n", target_color))
    return(NULL)
  }
  
  datME_curr <- MEs0[, ME_col_name, drop=FALSE]
  kME_curr <- abs(cor(datExpr[, mod_genes], datME_curr, use="p"))
  
  # 导出 Ranking
  hub_ranking <- data.frame(Gene_ID = rownames(kME_curr), kME = kME_curr[,1])
  hub_ranking <- hub_ranking[order(hub_ranking$kME, decreasing = TRUE), ]
  hub_ranking$Original_ID <- final_mapping$Original_ID[match(hub_ranking$Gene_ID, final_mapping$Valid_Rowname)]
  hub_ranking <- hub_ranking[, c("Original_ID", "Gene_ID", "kME")] 
  write.csv(hub_ranking, file.path(mod_out, paste0("Hub_Gene_Ranking_", target_color, ".csv")), row.names = FALSE)
  
  # Top 30 表达
  top_genes <- hub_ranking$Gene_ID[1:min(30, length(mod_genes))]
  expr_mat <- t(scale(datExpr[, top_genes]))
  write.csv(expr_mat, file.path(mod_out, paste0("Top_Hub_Genes_Expression_", target_color, ".csv")), row.names = TRUE)
  
  # [Visual Upgrade] 热图美化
  # 构建分组注释颜色
  groups_unique <- unique(meta_data$Group_Short_Name)
  group_colors <- setNames(pub_pal$group[1:length(groups_unique)], groups_unique)
  anno_colors <- list(Group = group_colors)
  
  pdf(file.path(mod_out, paste0("Heatmap_", target_color, ".pdf")), width=10, height=8)
  pheatmap(expr_mat, 
           annotation_col = data.frame(Group = meta_data$Group_Short_Name, row.names=rownames(datExpr)),
           annotation_colors = anno_colors, # 自定义分组颜色
           color = colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(50), # 统一红蓝配色
           show_colnames = FALSE, 
           border_color = NA, # 去除网格线，更现代
           main = paste("Top 30 Hub Genes - Module", target_color))
  dev.off()
  
  # 网络图
  if(length(top_genes) > 1) {
    adj_mat <- abs(cor(datExpr[, top_genes]))
    adj_mat[adj_mat < 0.3] <- 0; diag(adj_mat) <- 0
    g <- graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=TRUE)
    
    edge_df <- as_data_frame(g, what="edges")
    if(nrow(edge_df) > 0) {
      write.csv(edge_df, file.path(mod_out, paste0("Gene_Interaction_Edges_", target_color, ".csv")), row.names = FALSE)
      
      pdf(file.path(mod_out, paste0("Network_", target_color, ".pdf")), width=8, height=8)
      plot(g, layout = layout_with_fr(g), 
           vertex.label = V(g)$name, vertex.label.cex = 0.7, 
           vertex.color = target_color, # 节点颜色 = 模块颜色
           vertex.frame.color = "white", # 节点白边
           vertex.size = 15, 
           edge.width = E(g)$weight * 5,
           edge.color = "grey80",
           main = paste("Interaction Network - Module", target_color))
      dev.off()
    }
  }
  cat(sprintf("  -> Analyzed %s (Visuals Enhanced)\n", target_color))
}

# 自动选择模块
top_2_ME_names <- rownames(moduleTraitCor)[order(apply(moduleTraitCor, 1, max), decreasing=TRUE)[1:2]]
top_2_colors <- gsub("ME", "", top_2_ME_names)
target_list <- unique(c("magenta", "darkgreen", top_2_colors))

for(col in target_list) plot_module_analysis(col, "05_Hub_Genes")

#### 7. KEGG 引导文件 ####
# ==============================================================================
cat("\n[Step 7] Generating Enrichment Lists...\n")
for(mod in unique(moduleColors)) {
  genes <- colnames(datExpr)[moduleColors == mod]
  write.table(genes, file.path("06_Enrichment_Materials", paste0("Genes_", mod, ".txt")),
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

#### 8. 输出完整性验证 ####
# ==============================================================================
cat("\n[Step 8] Output Integrity Verification...\n")
required_files <- list(
  "01_Data_Clean/Sample_Distance_Matrix.csv",
  "01_Data_Clean/Sample_QC_Stats.csv",
  "01_Data_Clean/Gene_Name_Mapping_Full.csv",
  "01_Data_Clean/Gene_Name_Mapping_Filtered.csv",
  "02_Network_Topology/Soft_Threshold_Fit_Indices.csv",
  "03_Modules/Module_Assignment_Full.csv",
  "03_Modules/Module_Eigengenes.csv",
  "03_Modules/Module_Size_Statistics.csv",
  "04_Module_Trait_Correlations/Module_Trait_Correlation.csv",
  "04_Module_Trait_Correlations/Module_Trait_Pvalue.csv",
  "04_Module_Trait_Correlations/Module_Trait_Combined_Results.csv"
)

missing_files <- c()
for (f in required_files) {
  if (!file.exists(f)) {
    cat(sprintf("  ✗ MISSING: %s\n", f))
    missing_files <- c(missing_files, f)
  } else {
    cat(sprintf("  ✓ %s\n", f))
  }
}

required_pdfs <- c(
  "01_Data_Clean/Sample_Dendrogram_Trait_Heatmap_QC.pdf",
  "02_Network_Topology/Soft_Threshold_Plots.pdf",
  "03_Modules/Gene_Dendrogram_Final_All_Blocks.pdf",
  "04_Module_Trait_Correlations/Module_Trait_Heatmap.pdf"
)

for (pdf_file in required_pdfs) {
  if (!file.exists(pdf_file)) {
    cat(sprintf("  ✗ MISSING: %s\n", pdf_file))
    missing_files <- c(missing_files, pdf_file)
  } else {
    cat(sprintf("  ✓ %s\n", pdf_file))
  }
}

if (length(missing_files) == 0) {
  cat("\n✅ SUCCESS: All output files generated successfully!\n")
} else {
  cat(sprintf("\n⚠️  WARNING: %d files missing:\n", length(missing_files)))
  for (f in missing_files) cat(sprintf("    - %s\n", f))
}

cat("\n======================================================\n")
cat("   Analysis Complete!\n")
cat("======================================================\n")
