# ==============================================================================
#  Medicago sativa WGCNA 定制分析流程 (适配人工筛选基因集 + 表型关联)
#  输入对象: 
#    1. Summary_of_manually_filtered_genes (转录组)
#    2. Transcriptome_grouping_summary (分组信息)
# ==============================================================================

# ------------------------------------------------------------------------------
# [Step 1] 环境初始化与包加载
# ------------------------------------------------------------------------------
cat("\n[Step 1] 初始化环境...\n")
options(stringsAsFactors = FALSE)

# 自动检测并加载必要的包
required_packages <- c("WGCNA", "dplyr", "stringr", "ggplot2", "reshape2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "WGCNA" || pkg == "GO.db" || pkg == "impute") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# 开启多线程 (根据您的电脑配置，建议使用总核心数的 3/4)
enableWGCNAThreads(nThreads = 6) 

# 创建输出目录结构
dirs <- c("01_Data_Clean", "02_SoftThreshold", "03_Network", "04_Module_Trait", "05_Cytoscape")
for (d in dirs) if (!dir.exists(d)) dir.create(d)

# ------------------------------------------------------------------------------
# [Step 2] 数据读取与深度清洗 (核心定制部分)
# ------------------------------------------------------------------------------
cat("\n[Step 2] 加载并清洗数据...\n")

# 1. 检查对象是否存在
if (!exists("Summary_of_manually_filtered_genes")) stop("错误: 未找到转录组数据 'Summary_of_manually_filtered_genes'")
if (!exists("Transcriptome_grouping_summary")) stop("错误: 未找到分组数据 'Transcriptome_grouping_summary'")

# 2. 处理转录组数据
dat_raw <- Summary_of_manually_filtered_genes

# 识别 Sample:fpkm 列 (兼容 :fpkm 或 .fpkm)
fpkm_cols <- grep("(:|\\.)fpkm$", colnames(dat_raw), value = TRUE, ignore.case = TRUE)
if(length(fpkm_cols) == 0) stop("错误: 未能在数据中找到以 ':fpkm' 结尾的列。")

cat(sprintf("   -> 识别到 %d 个样本列 (自动忽略 NR 列)\n", length(fpkm_cols)))

# 提取表达矩阵
datExpr0 <- as.data.frame(t(dat_raw[, fpkm_cols]))
colnames(datExpr0) <- dat_raw[[1]] # 第一列为 Gene_ID

# 清洗样本名称 (去除 :fpkm 后缀)
clean_names <- str_replace_all(rownames(datExpr0), "(:|\\.)fpkm$", "")
rownames(datExpr0) <- clean_names
cat("   -> 样本名清洗示例:", paste(head(clean_names, 3), collapse=", "), "\n")

# 3. 处理表型/分组数据并对齐
trait_data <- Transcriptome_grouping_summary

# 假设第一列是 Sample 名
trait_sample_col <- colnames(trait_data)[1] 
# 确保表型数据的样本列与清洗后的表达谱样本名一致
common_samples <- intersect(rownames(datExpr0), trait_data[[trait_sample_col]])

if(length(common_samples) < nrow(datExpr0)) {
  cat("   ⚠️ 注意: 表达谱中有部分样本在分组表中未找到，将只分析共有样本。\n")
}

# 重新取子集并排序
datExpr <- datExpr0[common_samples, ]
trait_data <- trait_data[match(common_samples, trait_data[[trait_sample_col]]), ]

# 4. 预处理：Log2 转换 (FPKM 通常需要 log 转换以符合正态分布假设)
# 使用 log2(x + 1) 防止 log(0)
datExpr <- log2(datExpr + 1)
cat("   -> 数据已完成 Log2 转化\n")

# ------------------------------------------------------------------------------
# [Step 3] 构建表型数值矩阵 (用于关联分析)
# ------------------------------------------------------------------------------
cat("\n[Step 3] 构建表型特征矩阵 (0/1编码)...\n")

# 创建一个空的数值矩阵
datTraits <- data.frame(row.names = rownames(datExpr))

# --- A. 自动根据 Group_Short_Name 生成分组特征 ---
# 例如: Group_LLGS, Group_SLRS 等
groups <- unique(trait_data$Group_Short_Name)
for(g in groups){
  datTraits[[paste0("Group_", g)]] <- ifelse(trait_data$Group_Short_Name == g, 1, 0)
}

# --- B. 解析生物学特征 (品系与光周期) ---
# 根据样本名或分组名中的关键字进行提取
# 逻辑：
#   GS = Green Stem (绿茎)
#   RS = Red Stem (红茎)
#   LL = Long Light (长光照)
#   SL = Short Light (短光照)

# 1. 品系特征 (Stem Color)
datTraits$Stem_Green <- ifelse(grepl("GS", trait_data$Group_Short_Name), 1, 0)
datTraits$Stem_Red   <- ifelse(grepl("RS", trait_data$Group_Short_Name), 1, 0)

# 2. 光周期特征 (Photoperiod)
datTraits$Light_Long  <- ifelse(grepl("LL", trait_data$Group_Short_Name), 1, 0)
datTraits$Light_Short <- ifelse(grepl("SL", trait_data$Group_Short_Name), 1, 0)

# 3. 移除方差为0的列 (全0或全1的列，无法计算相关性)
datTraits <- datTraits[, apply(datTraits, 2, var) > 0]

cat("   -> 已生成以下表型特征用于关联分析:\n")
print(colnames(datTraits))

# 保存清洗后的数据备份
save(datExpr, datTraits, file = "01_Data_Clean/Cleaned_Input_Data.RData")

# ------------------------------------------------------------------------------
# [Step 4] 软阈值筛选 (Soft Thresholding)
# ------------------------------------------------------------------------------
cat("\n[Step 4] 计算软阈值...\n")

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned")

# 自动判断最佳 Power
# 优先选择 R^2 > 0.85 的最小 power
estimate_power <- sft$powerEstimate
if (is.na(estimate_power)) {
  # 如果自动判定失败，手动查找第一个 > 0.8 的
  idx <- which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] > 0.80)
  if(length(idx) > 0) estimate_power <- powers[idx[1]] else estimate_power <- 8 # 默认兜底
}
cat(sprintf("   -> 选定的软阈值 Power = %d\n", estimate_power))

# 绘图
pdf("02_SoftThreshold/Soft_Threshold_Plots.pdf", width=9, height=5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, R^2",
     type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

# ------------------------------------------------------------------------------
# [Step 5] 构建共表达网络 (One-step Method)
# ------------------------------------------------------------------------------
cat("\n[Step 5] 构建网络与模块识别...\n")

# 因为只有 ~5486 个基因，我们可以设置较大的 maxBlockSize 一次性算完
net <- blockwiseModules(datExpr,
                        power = estimate_power,
                        maxBlockSize = 8000,      # 设大一点以覆盖所有基因
                        TOMType = "unsigned", 
                        minModuleSize = 30,       # 最小模块基因数
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25,    # 合并相似度 > 0.75 的模块
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "03_Network/Medicago_TOM",
                        verbose = 3)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# 保存结果
save(net, moduleLabels, moduleColors, MEs, file = "03_Network/Network_Result.RData")

# 绘制层级聚类图
pdf("03_Network/Module_Dendrogram.pdf", width=10, height=6)
plotDendroAndColors(geneTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

cat("\n   -> 模块统计:\n")
print(table(moduleColors))

# ------------------------------------------------------------------------------
# [Step 6] 模块与表型关联分析 (关键步骤)
# ------------------------------------------------------------------------------
cat("\n[Step 6] 分析模块与表型的相关性...\n")

# 1. 计算相关性
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0) # 重新排序特征基因

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# 2. 绘制热图
pdf("04_Module_Trait/Module_Trait_Heatmap.pdf", width=10, height=min(12, nrow(moduleTraitCor)*0.5 + 4))

# 准备热图文本 (相关系数 + P值)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 9, 3, 3)) # 调整下、左、上、右边距
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), # 蓝-白-红 渐变
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))
dev.off()

# 3. 导出相关性表格
trait_results <- data.frame(Module = rownames(moduleTraitCor))
for (col in colnames(datTraits)) {
  trait_results[[paste0(col, "_Cor")]] <- moduleTraitCor[, col]
  trait_results[[paste0(col, "_Pval")]] <- moduleTraitPvalue[, col]
}
write.csv(trait_results, "04_Module_Trait/Module_Trait_Correlation_Table.csv", row.names = FALSE)

cat("   -> 关联热图已保存至 04_Module_Trait 文件夹\n")

# ------------------------------------------------------------------------------
# [Step 7] 导出结果与 Cytoscape 文件
# ------------------------------------------------------------------------------
cat("\n[Step 7] 导出基因列表与网络文件...\n")

# 1. 导出基因-模块对应表
gene_info <- data.frame(GeneID = colnames(datExpr),
                        ModuleLabel = moduleLabels,
                        ModuleColor = moduleColors)
write.csv(gene_info, "03_Network/Gene_Module_List.csv", row.names = FALSE)

# 2. 导出 Cytoscape 网络文件 (筛选 Hub 基因以防卡死)
# 逻辑：只导出那些在热图中显著的模块，或者全部导出
# 这里演示导出全部非灰色模块，每个模块取 kME 前 150 个点

modules_to_export <- unique(moduleColors)
modules_to_export <- modules_to_export[modules_to_export != "grey"]

for (mod in modules_to_export) {
  # 获取该模块基因
  inModule <- is.finite(match(moduleColors, mod))
  modGenes <- colnames(datExpr)[inModule]
  
  # 如果基因太多，只取 Hub 基因 (基于 kME)
  if (length(modGenes) > 150) {
    # 计算该模块内的 kME
    ME_name <- paste0("ME", mod)
    datME <- MEs[, ME_name]
    gene_kME <- abs(cor(datExpr[, modGenes], datME, use="p"))
    # 排序取前 150
    ranked_genes <- rownames(gene_kME)[order(-gene_kME, decreasing = TRUE)]
    modGenes <- ranked_genes[1:150]
  }
  
  # 重新加载 TOM (如果有多个 block 这里的逻辑比较复杂，但因为我们设了 maxBlockSize > geneNum，只有一个 TOM)
  # 为简单起见，这里重新计算小子集的 TOM，速度很快
  modExpr <- datExpr[, modGenes]
  TOM <- TOMsimilarityFromExpr(modExpr, power = estimate_power, networkType = "unsigned")
  
  # 设定显示的边阈值 (Top 20% 强边)
  threshold <- quantile(TOM[lower.tri(TOM)], probs = 0.80)
  
  # 导出
  cytoscape_file_edge <- paste0("05_Cytoscape/Cytoscape_Edge_", mod, ".txt")
  cytoscape_file_node <- paste0("05_Cytoscape/Cytoscape_Node_", mod, ".txt")
  
  exportNetworkToCytoscape(TOM,
                           edgeFile = cytoscape_file_edge,
                           nodeFile = cytoscape_file_node,
                           weighted = TRUE,
                           threshold = threshold,
                           nodeNames = modGenes,
                           nodeAttr = moduleColors[inModule][1:length(modGenes)])
}

cat("\n[Done] 分析全部完成！请查看工作目录下的文件夹。\n")
