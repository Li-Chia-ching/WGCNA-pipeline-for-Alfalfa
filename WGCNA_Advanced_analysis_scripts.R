# ==============================================================================
# WGCNA 进阶分析脚本 - v2.0 (修复版)
# 功能：可视化、富集分析、核心基因挖掘
# 修复：
#   1. 修正了"模块-性状"热图的命名错误，改为"模块-模块"相关性。
#   2. 增加了 RStudio 代码折叠目录 (####)。
#   3. 清晰区分了网络图、热图和聚类图的输出。
# ==============================================================================

#### 1. 初始化与包加载 ####
# ==============================================================================
cat("╔══════════════════════════════════════╗\n")
cat("║      WGCNA 进阶分析 v2.0 启动        ║\n")
cat("╚══════════════════════════════════════╝\n")

options(stringsAsFactors = FALSE)

# 绘图与数据处理
library(WGCNA)
library(ggplot2)
library(pheatmap)   # 绘制热图
library(igraph)     # 绘制网络图
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(patchwork)

# 富集分析 (需确保已安装)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db) # 注意：如果是其他物种，请修改此处

# 设置多线程
enableWGCNAThreads(nThreads = ifelse(parallel::detectCores() > 4, 6, 4))

# 设置绘图主题
theme_set(theme_minimal(base_size = 11) + 
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", hjust = 0.5)))

#### 2. 数据加载与校验 ####
# ==============================================================================
cat("\n[Step 2] 加载并验证 Step1 的分析结果...\n")

# 定义加载函数
load_and_validate_data <- function() {
  # 1. 加载 RData
  if(!file.exists("03_Network/WGCNA_Network_Object.RData")) {
    stop("❌ 错误：找不到 03_Network/WGCNA_Network_Object.RData")
  }
  load("03_Network/WGCNA_Network_Object.RData")
  
  # 2. 加载表达矩阵
  if(!file.exists("01_InputData/Preprocessed_Expression_Matrix.csv")) {
    stop("❌ 错误：找不到 01_InputData/Preprocessed_Expression_Matrix.csv")
  }
  datExpr <- read.csv("01_InputData/Preprocessed_Expression_Matrix.csv", row.names = 1)
  datExpr <- as.data.frame(t(datExpr)) 
  
  # 3. 修复 module_colors 名字丢失问题
  if (is.null(names(module_colors))) {
    if (length(module_colors) == ncol(datExpr)) {
      names(module_colors) <- colnames(datExpr)
      cat("  ✅ 已修复 module_colors 的基因名称。\n")
    } else {
      stop("❌ module_colors 长度与基因数不匹配，无法继续。")
    }
  }
  
  # 4. 对齐数据
  common_genes <- intersect(colnames(datExpr), names(module_colors))
  datExpr <- datExpr[, common_genes, drop = FALSE]
  module_colors <- module_colors[common_genes]
  
  cat(paste("  ✅ 数据加载完成: ", nrow(datExpr), "样本 x", ncol(datExpr), "基因\n"))
  return(list(datExpr = datExpr, module_colors = module_colors, net = net, MEs = MEs, softPower = softPower))
}

# 执行加载
data_list <- load_and_validate_data()
list2env(data_list, envir = .GlobalEnv)

# 创建目录结构
dir.create("06_Advanced_Results", showWarnings = FALSE)
dir.create("06_Advanced_Results/Figures", showWarnings = FALSE)
dir.create("06_Advanced_Results/Tables", showWarnings = FALSE)
dir.create("06_Advanced_Results/Enrichment", showWarnings = FALSE)

#### 3. 可视化：基因层次聚类树 (Dendrogram) ####
# ==============================================================================
# 目的：展示所有基因是如何被切分成不同颜色模块的
cat("\n[Step 3] 绘制基因聚类与模块切割图...\n")

pdf("06_Advanced_Results/Figures/1_Gene_Dendrogram_and_Modules.pdf", width = 12, height = 7)
plotDendroAndColors(net$dendrograms[[1]], 
                    module_colors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Assignment")
dev.off()
cat("  -> 输出: Figures/1_Gene_Dendrogram_and_Modules.pdf\n")

#### 4. 可视化：软阈值拓扑分析 (Topology) ####
# ==============================================================================
# 目的：验证之前选择的 Power 值是否符合无尺度网络分布
cat("\n[Step 4] 绘制网络拓扑分析图...\n")

# 为了绘图需快速重算一次统计量
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 0)
sft_df <- data.frame(Power = sft$fitIndices$Power,
                     R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                     MeanK = sft$fitIndices$mean.k.)

# 绘图
p1 <- ggplot(sft_df, aes(x=Power, y=R2)) + 
  geom_point(color="red") + geom_line(color="red") +
  geom_hline(yintercept = 0.85, linetype="dashed") +
  geom_vline(xintercept = softPower, linetype="dashed", color="blue") +
  labs(title="Scale Independence", y="Scale Free Topology Model Fit (R^2)")

p2 <- ggplot(sft_df, aes(x=Power, y=MeanK)) + 
  geom_point(color="blue") + geom_line(color="blue") +
  scale_y_log10() +
  labs(title="Mean Connectivity", y="Mean Connectivity")

combined_plot <- p1 + p2
ggsave("06_Advanced_Results/Figures/2_Soft_Threshold_Topology.pdf", combined_plot, width=12, height=6)
cat("  -> 输出: Figures/2_Soft_Threshold_Topology.pdf\n")

#### 5. 功能富集分析 (Enrichment) ####
# ==============================================================================
# 目的：对主要模块进行 KEGG/GO 功能注释
cat("\n[Step 5] 进行功能富集分析...\n")

# 定义富集函数
run_enrichment <- function(module, genes, db, code) {
  # 转换ID
  ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=db)
  if(nrow(ids) < 5) return(NULL)
  
  # KEGG
  kk <- enrichKEGG(ids$ENTREZID, organism=code, pvalueCutoff=0.05)
  if(!is.null(kk) && nrow(kk) > 0) {
    write.csv(kk@result, paste0("06_Advanced_Results/Enrichment/KEGG_", module, ".csv"))
    p <- dotplot(kk, showCategory=15) + ggtitle(paste("KEGG:", module))
    ggsave(paste0("06_Advanced_Results/Figures/Enrich_KEGG_", module, ".pdf"), p, width=8, height=7)
  }
}

# 选取前5大模块进行分析
mod_stats <- sort(table(module_colors), decreasing=TRUE)
target_mods <- names(mod_stats)[names(mod_stats) != "grey"][1:min(5, length(mod_stats))]

for(mod in target_mods) {
  cat(paste("  正在分析模块:", mod, "...\n"))
  genes <- names(module_colors)[module_colors == mod]
  # 注意：这里默认是人类 (hsa, org.Hs.eg.db)。如需其他物种请手动修改
  tryCatch({
    run_enrichment(mod, genes, org.Hs.eg.db, "hsa")
  }, error=function(e) { cat("   警告: 富集分析失败 (可能是网络或物种包问题)\n") })
}

#### 6. 可视化：核心基因互作网络与热图 (Hub Genes) ####
# ==============================================================================
# 目的：展示模块内部最核心的基因及其相互作用
cat("\n[Step 6] 绘制核心基因互作网络与热图...\n")

analyze_hub_genes <- function(target_module, datExpr, MEs, n_top=20) {
  # 1. 准备数据
  genes <- names(module_colors)[module_colors == target_module]
  ME_vec <- MEs[, paste0("ME", target_module)]
  
  # 2. 计算 kME (模块身份度)
  kME <- signedKME(datExpr[, genes], MEs[, paste0("ME", target_module), drop=FALSE])
  colnames(kME) <- "kME"
  top_genes <- rownames(kME)[order(-abs(kME$kME))][1:min(n_top, length(genes))]
  
  # 3. 绘制热图 (Expression Heatmap)
  expr_mat <- t(scale(datExpr[, top_genes])) # Z-score 标准化
  pheatmap(expr_mat, 
           main = paste("Top", n_top, "Hub Genes Expression -", target_module),
           cluster_rows = TRUE, cluster_cols = TRUE,
           fontsize_row = 8,
           filename = paste0("06_Advanced_Results/Figures/3_Heatmap_Hub_", target_module, ".pdf"),
           width = 8, height = 8)
  
  # 4. 绘制网络图 (Interaction Network)
  # 计算这些 Top 基因间的相关性
  adj_mat <- abs(cor(datExpr[, top_genes]))
  adj_mat[adj_mat < 0.5] <- 0 # 只保留强相关
  diag(adj_mat) <- 0
  
  g <- graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted=TRUE)
  if(ecount(g) > 0) {
    pdf(paste0("06_Advanced_Results/Figures/3_Network_Hub_", target_module, ".pdf"), width=8, height=8)
    plot(g, vertex.label=V(g)$name, 
         vertex.size=15, 
         vertex.color=target_module,
         vertex.label.cex=0.7,
         edge.width=E(g)$weight*3,
         layout=layout_with_fr(g),
         main=paste("Hub Gene Network -", target_module))
    dev.off()
  }
}

# 对前3大模块进行绘图
for(mod in target_mods[1:3]) {
  cat(paste("  绘制核心图:", mod, "\n"))
  analyze_hub_genes(mod, datExpr, MEs)
}

#### 7. 可视化：模块间相关性热图 (Module-Module) ####
# ==============================================================================
# 目的：展示不同模块之间的表达模式是否相似 (用于判断是否需要合并)
# 注意：这*不是*表型关联图，因为没有表型数据。
cat("\n[Step 7] 绘制模块-模块特征基因相关性热图...\n")

if(exists("MEs") && ncol(MEs) > 1) {
  # 去除 grey 模块
  MEs_clean <- MEs[, !grepl("grey", colnames(MEs))]
  
  # 计算相关性
  mod_cor <- cor(MEs_clean)
  
  # 绘图
  pheatmap(mod_cor,
           main = "Module-Module Eigengene Correlation (Similarity)", # 修正后的标题
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 6,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           filename = "06_Advanced_Results/Figures/4_Module_Module_Correlation.pdf",
           width = 10, height = 10)
  
  cat("  -> 输出: Figures/4_Module_Module_Correlation.pdf (已修正标题)\n")
} else {
  cat("  警告: 模块数量不足，跳过模块相关性热图。\n")
}

#### 8. 完成总结 ####
# ==============================================================================
cat("\n╔══════════════════════════════════════╗\n")
cat("║       分析完成！结果已分类保存       ║\n")
cat("╚══════════════════════════════════════╝\n")
cat("输出文件说明:\n")
cat("1. 基因聚类图 (Figures/1_...): 展示基因如何被分配到模块。\n")
cat("2. 拓扑分析图 (Figures/2_...): 验证软阈值选择是否合理。\n")
cat("3. 核心基因图 (Figures/3_...): 包含特定模块的【热图】和【互作网络图】。\n")
cat("4. 模块相关图 (Figures/4_...): 展示模块间的相似性 (非性状关联)。\n")
