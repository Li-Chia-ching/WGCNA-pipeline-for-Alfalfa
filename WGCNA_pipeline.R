# ============================================
# WGCNA 全流程分析脚本 (适配 R 4.5.2 + 清华镜像)
# 作者：李嘉庆
# 日期：2025年12月5日(v1.1.0)
# ============================================

cat("
╔══════════════════════════════════════╗
║        WGCNA 全流程分析启动          ║
║        环境：R 4.5.2 | Win10         ║
╚══════════════════════════════════════╝\n")

# ---------- 第一部分：初始化与包管理 (稳健化) ----------
cat("\n[1/7] 初始化环境与安装依赖包...\n")
options(stringsAsFactors = FALSE)
options(timeout = 600) # 超时设置为10分钟

# 设置清华镜像（已配置成功，此处显式声明确保一致性）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioc/")
cat("   -> 镜像已设置为清华源\n")

# 创建有序的输出目录
output_dirs <- c("01_InputData", "02_QC", "03_Network", "04_Modules", "05_Results")
for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = FALSE)
    cat("   -> 创建目录:", dir, "\n")
  }
}

# 定义包安装函数（支持重试）
install_with_retry <- function(pkg, source = "CRAN", retry = 2) {
  installed <- FALSE
  attempt <- 1
  
  while (!installed && attempt <= retry) {
    cat("   尝试安装 (第", attempt, "次): ", pkg, " [", source, "]\n", sep="")
    tryCatch({
      if (source == "CRAN") {
        install.packages(pkg, quiet = TRUE, dependencies = TRUE)
      } else if (source == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) 
          install.packages("BiocManager", quiet = TRUE)
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
      }
      installed <- TRUE
      cat("     √ 成功安装\n")
    }, error = function(e) {
      cat("     × 尝试失败: ", e$message, "\n")
      attempt <- attempt + 1
      if (attempt > retry) warning("包 ‘", pkg, "’ 安装最终失败，可能影响部分功能。")
    })
  }
  return(installed)
}

# 安装并加载核心包 (修复版)
cat("\n  处理核心分析包...\n")
required_cran <- c("WGCNA", "dplyr", "ggplot2", "reshape2", "stringr", "RColorBrewer")
required_bioc <- c("GO.db", "impute")

# 修复1：安装CRAN包 (使用显式循环，避免“pkg”未定义)
cat("  安装CRAN包:\n")
for (current_pkg in required_cran) {
  if (!requireNamespace(current_pkg, quietly = TRUE)) {
    cat("    正在安装:", current_pkg, "... ")
    # 调用之前定义好的安装函数，并传入具体的包名
    success <- install_with_retry(current_pkg, source = "CRAN", retry = 2)
    if(success) {
      cat("完成\n")
    } else {
      cat("失败\n")
    }
  } else {
    cat("    ", current_pkg, "已安装，跳过\n", sep="")
  }
}

# 修复2：安装Bioconductor包
cat("\n  处理Bioconductor包:\n")
for (current_pkg in required_bioc) {
  if (!requireNamespace(current_pkg, quietly = TRUE)) {
    cat("    正在安装:", current_pkg, "... ")
    success <- install_with_retry(current_pkg, source = "Bioconductor", retry = 2)
    if(success) {
      cat("完成\n")
    } else {
      cat("失败\n")
    }
  } else {
    cat("    ", current_pkg, "已安装，跳过\n", sep="")
  }
}

# 安全加载包
cat("\n  加载包至会话...\n")
suppressPackageStartupMessages({
  library(WGCNA)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})
cat("   -> 核心包加载完成\n")

# 启用多线程（为8核CPU优化）
enableWGCNAThreads(nThreads = 6)
cat("   -> 并行计算已启用 (使用6个线程)\n")

# ---------- 第二部分：数据加载与验证 (修复版) ----------
cat("\n[2/7] 加载并验证表达数据...\n")

# 检查数据对象是否存在
if (!exists("Expression_with_annotation")) {
  stop("❌ 错误: 全局环境中未找到 ‘Expression_with_annotation’ 对象。\n",
       "请确保：\n",
       "1. 已使用 readr/read.csv 正确导入CSV文件\n",
       "2. 对象名称拼写完全一致\n",
       "3. 运行命令类似: Expression_with_annotation <- read.csv('your_file.csv', sep='\\t')")
}
dat0 <- Expression_with_annotation
cat("   -> 数据框加载成功，维度:", dim(dat0)[1], "行 x", dim(dat0)[2], "列\n")

# 自动识别基因ID列和FPKM列（关键步骤）
cat("   -> 自动识别列结构...\n")
gene_id_column <- NULL
possible_gene_ids <- c("Gene_ID", "gene_id", "GeneID", "Geneid", "geneId", colnames(dat0)[1])
for (id in possible_gene_ids) {
  if (id %in% colnames(dat0)) {
    gene_id_column <- id
    break
  }
}
if (is.null(gene_id_column)) stop("未找到基因ID列，请检查数据")
cat("       基因ID列: ‘", gene_id_column, "’\n", sep="")

# 识别所有FPKM列（兼容冒号或点号格式）- 修复正则表达式
fpkm_patterns <- c(":fpkm$", "\\\\.fpkm$", "_fpkm$", ":FPKM$", "\\\\.FPKM$")
fpkm_columns <- unique(unlist(sapply(fpkm_patterns, 
                                     function(p) grep(p, colnames(dat0), value=TRUE, ignore.case=TRUE))))

if (length(fpkm_columns) < 3) {
  cat("⚠️  警告: 仅找到", length(fpkm_columns), "个FPKM列。前10个列名为:\n")
  print(head(colnames(dat0), 10))
  stop("FPKM列数不足，请确认数据格式是否为‘样本:fpkm’或‘样本.fpkm’")
}
cat("       找到", length(fpkm_columns), "个FPKM表达量列\n")

# 提取表达矩阵
expr_raw <- as.matrix(dat0[, fpkm_columns, drop=FALSE])
rownames(expr_raw) <- as.character(dat0[[gene_id_column]])

# 清理样本名（修复正则表达式错误）
clean_sample_names <- function(names) {
  # 使用str_remove_all，更安全
  names <- stringr::str_remove_all(names, "(?i)(:fpkm$|\\\\.fpkm$|_fpkm$)")
  # 进一步清理：移除可能的其他后缀
  names <- stringr::str_remove_all(names, "(?i)(:read count$|\\\\.read count$)")
  return(names)
}

sample_names <- clean_sample_names(fpkm_columns)
colnames(expr_raw) <- sample_names
cat("   -> 表达矩阵创建完成，样本名示例:", paste(head(sample_names, 3), collapse=", "), "\n")

# ---------- 第三部分：数据预处理与严格过滤 ----------
cat("\n[3/7] 数据预处理与质控过滤...\n")

# 转换并处理缺失值
expr_numeric <- apply(expr_raw, 2, as.numeric)
rownames(expr_numeric) <- rownames(expr_raw)
expr_numeric[is.na(expr_numeric) | is.infinite(expr_numeric)] <- 0
cat("   -> 缺失值/无限值已处理\n")

# 过滤步骤1：去除在超过50%样本中表达量为0的基因
min_samples <- ceiling(ncol(expr_numeric) * 0.5)
keep <- rowSums(expr_numeric > 1) >= min_samples
expr_f1 <- expr_numeric[keep, ]
cat("   -> 低表达过滤: ", sum(keep), "/", nrow(expr_numeric), " 个基因保留\n")

# 过滤步骤2：去除表达方差最小的后25%基因（大幅提升速度）
gene_var <- apply(expr_f1, 1, var)
var_cut <- quantile(gene_var, probs=0.25)
keep_var <- gene_var > var_cut
expr_final <- expr_f1[keep_var, ]
cat("   -> 低变异过滤: ", sum(keep_var), "/", nrow(expr_f1), " 个基因保留\n")
cat("   -> 最终用于分析的基因数: ", nrow(expr_final), "\n")

# 保存过滤统计
qc_stats <- data.frame(
  Stage = c("原始数据", "低表达过滤后", "低变异过滤后"),
  Genes = c(nrow(expr_numeric), nrow(expr_f1), nrow(expr_final)),
  Samples = rep(ncol(expr_numeric), 3)
)
write.csv(qc_stats, file="01_InputData/QC_Filtering_Stats.csv", row.names=FALSE)

# Log2转换
datExpr_log2 <- log2(expr_final + 1) # 伪计数+1

# 保存预处理后的矩阵
write.csv(datExpr_log2, 
          file="01_InputData/Preprocessed_Expression_Matrix.csv",
          quote=FALSE)
cat("   -> 预处理矩阵已保存\n")

# 样本聚类检查离群值
cat("   -> 进行样本层级聚类...\n")
sample_tree <- hclust(dist(t(datExpr_log2)), method="average")
pdf("02_QC/Sample_Clustering_Outliers.pdf", width=12, height=6)
plot(sample_tree, main="Sample Clustering to Detect Outliers", 
     xlab="", sub="", cex=0.7)
abline(h=mean(sample_tree$height) + 2*sd(sample_tree$height), col="red", lty=2)
dev.off()
cat("   -> 样本聚类图已生成\n")

# ---------- 第四部分：软阈值选择 ----------
cat("\n[4/7] 选择网络构建软阈值...\n")
datExpr_forWGCNA <- as.data.frame(t(datExpr_log2))

# 测试一系列软阈值
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr_forWGCNA, 
                         powerVector=powers, 
                         networkType="unsigned",
                         verbose=5, 
                         corFnc="cor", # 使用标准相关性，更快
                         corOptions=list(use='p', method='pearson'))

# 自动选择或推荐软阈值
if (!is.na(sft$powerEstimate)) {
  softPower <- sft$powerEstimate
  cat("   -> 自动选择的软阈值 power =", softPower, "\n")
} else {
  # 找出R²首次超过0.85的阈值
  R_sq <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  idx <- which(R_sq >= 0.85)
  if (length(idx) > 0) {
    softPower <- sft$fitIndices[min(idx), 1]
    cat("   -> 基于R²≥0.85选择软阈值 power =", softPower, "\n")
  } else {
    softPower <- 6
    cat("   -> 未找到合适阈值，使用默认值 power =", softPower, "\n")
  }
}

# 绘制软阈值选择图
png("02_QC/Soft_Threshold_Selection.png", width=1000, height=500)
par(mfrow=c(1,2))
# 图1：无标度拓扑拟合
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, R²",
     type="n", main="Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.85, col="red", lty=2)
# 图2：平均连接度
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity",
     type="n", main="Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()
cat("   -> 软阈值选择图已保存\n")

# ---------- 第五部分：构建共表达网络与模块识别 ----------
cat("\n[5/7] 构建基因共表达网络 (此步骤最耗时)...\n")
cat("      开始时间: ", format(Sys.time(), "%H:%M:%S"), "\n")

# 根据数据规模调整参数（性能优化关键）
n_genes <- ncol(datExpr_forWGCNA)
if (n_genes > 8000) {
  maxBlockSize <- 6000; minModuleSize <- 50
} else if (n_genes > 3000) {
  maxBlockSize <- 4000; minModuleSize <- 40
} else {
  maxBlockSize <- 2000; minModuleSize <- 30
}
cat("   -> 参数设置: 基因数=", n_genes, 
    ", 分块大小=", maxBlockSize, 
    ", 最小模块大小=", minModuleSize, "\n", sep="")

# 核心网络构建函数
net <- blockwiseModules(datExpr_forWGCNA,
                        power = softPower,
                        maxBlockSize = maxBlockSize,     # 分块计算防内存溢出
                        TOMType = "unsigned",
                        minModuleSize = minModuleSize,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,            # 模块用数字标签
                        pamRespectsDendro = FALSE,       # 设为FALSE可加快计算
                        saveTOMs = TRUE,                 # 保存TOM矩阵
                        saveTOMFileBase = "03_Network/WGCNA_TOM",
                        verbose = 3,
                        nThreads = 6)                    # 使用6个线程

cat("      结束时间: ", format(Sys.time(), "%H:%M:%S"), "\n")
cat("   -> 网络构建完成！\n")

# 提取结果
module_labels <- net$colors
module_colors <- labels2colors(module_labels)
MEs <- net$MEs

# 模块统计
module_stats <- as.data.frame(table(module_colors))
colnames(module_stats) <- c("ModuleColor", "GeneCount")
module_stats <- module_stats[order(-module_stats$GeneCount), ]
cat("\n  模块大小统计:\n")
print(module_stats)

# 保存网络对象
save(net, module_colors, MEs, softPower, 
     file="03_Network/WGCNA_Network_Object.RData")
cat("   -> 网络对象已保存\n")

# ---------- 第六部分：结果可视化 ----------
# ====== 修复方案：生成模块特征基因热图 ======
cat("\n[修复运行] 重新生成模块特征基因相关热图...\n")

# 1. 确保WGCNA包已加载
if (!"WGCNA" %in% .packages()) {
  cat("正在加载WGCNA包...\n")
  library(WGCNA)
}

# 2. 检查所需对象是否存在
if (exists("MEs")) {
  cat("MEs对象存在，维度:", dim(MEs), "\n")
} else {
  # 尝试从保存的文件中加载
  if (file.exists("03_Network/WGCNA_Network_Object.RData")) {
    cat("从文件加载网络对象...\n")
    load("03_Network/WGCNA_Network_Object.RData")
  } else {
    stop("错误: 未找到MEs对象，请确保已成功运行网络构建步骤。")
  }
}

# 3. 使用完整的命名空间调用函数（核心修复）
if (ncol(MEs) > 2) {
  cat("正在对模块特征基因进行排序...\n")
  
  # 方法1: 使用 WGCNA:: 显式调用
  MEs_sorted <- WGCNA::orderMEs(MEs)
  
  # 如果仍然失败，尝试方法2: 直接使用函数定义
  if (!exists("MEs_sorted")) {
    cat("尝试备用方法...\n")
    # orderMEs的简化实现
    MEs_sorted <- MEs
    if (ncol(MEs) > 1) {
      corME <- cor(MEs_sorted, use = "p")
      disME <- as.dist(1 - corME)
      clustME <- hclust(disME, method = "average")
      MEs_sorted <- MEs_sorted[, clustME$order]
    }
  }
  
  # 计算相关性矩阵
  module_cor <- cor(MEs_sorted, use="p")
  
  # 生成热图
  cat("正在生成热图...\n")
  png("04_Modules/Module_Eigengene_Correlation.png", width=800, height=800)
  
  # 使用更可靠的heatmap调用方式
  tryCatch({
    heatmap(module_cor, 
            symm = TRUE,
            main = "Correlation between Module Eigengenes",
            margins = c(12, 12),  # 底部和左侧边距
            cexRow = 0.8, 
            cexCol = 0.8,
            col = colorRampPalette(c("blue", "white", "red"))(50))
  }, error = function(e) {
    cat("热图生成警告:", e$message, "\n")
    # 备用方案: 使用基础plot
    plot(1, type="n", xlab="", ylab="", main="Module Correlation Heatmap")
    text(1, 1, "热图生成失败，请检查数据")
  })
  
  dev.off()
  cat("   -> 模块特征基因相关热图已保存: 04_Modules/Module_Eigengene_Correlation.png\n")
  
  # 保存相关性矩阵
  write.csv(module_cor, 
            file = "04_Modules/Module_Eigengene_Correlation_Matrix.csv",
            quote = FALSE)
  cat("   -> 模块特征基因相关性矩阵已保存\n")
  
} else {
  cat("   -> 模块特征基因数量不足 (只有", ncol(MEs), "个)，跳过热图生成。\n")
}

# 4. 验证修复结果
cat("\n[验证] 检查函数可用性:\n")
cat("orderMEs函数存在:", exists("orderMEs"), "\n")
cat("WGCNA::orderMEs存在:", exists("orderMEs", where = asNamespace("WGCNA")), "\n")

# ---------- 第七部分：结果导出与整合 ----------
cat("\n[7/7] 导出结果文件...\n")

# 1. 基因-模块对应关系（包含原始注释）
gene_module_assign <- data.frame(
  GeneID = colnames(datExpr_forWGCNA),
  ModuleLabel = module_labels,
  ModuleColor = module_colors,
  stringsAsFactors=FALSE
)

# 合并原始注释信息（如果存在）
annotation_cols <- setdiff(colnames(dat0), c(gene_id_column, fpkm_columns))
if (length(annotation_cols) > 0) {
  # 提取注释
  gene_annot <- dat0[, c(gene_id_column, annotation_cols[1:min(3, length(annotation_cols))]), drop=FALSE]
  colnames(gene_annot)[1] <- "GeneID"
  # 合并
  gene_module_assign <- merge(gene_module_assign, gene_annot, by="GeneID", all.x=TRUE)
  cat("   -> 已合并", length(annotation_cols), "列注释信息\n")
}

write.csv(gene_module_assign, 
          file="05_Results/Gene_Module_Assignments_Full.csv",
          row.names=FALSE, quote=FALSE)
cat("   -> 基因-模块分配表已保存 (包含", nrow(gene_module_assign), "个基因)\n")

# 2. 按模块输出基因列表
for (color in unique(module_colors)) {
  if (color == "grey") next
  genes_in_module <- gene_module_assign$GeneID[gene_module_assign$ModuleColor == color]
  if (length(genes_in_module) > 5) {
    write.table(genes_in_module,
                file=paste0("05_Results/Module_", color, "_GeneList.txt"),
                row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
}
cat("   -> 各模块基因列表已保存\n")

# 3. 模块特征基因表达量
MEs_with_names <- orderMEs(MEs)
rownames(MEs_with_names) <- rownames(datExpr_forWGCNA)
write.csv(MEs_with_names, 
          file="05_Results/Module_Eigengenes_Expression.csv",
          quote=FALSE)
cat("   -> 模块特征基因表达量表已保存\n")

# 4. 分析总结报告
sink("05_Results/WGCNA_Analysis_Summary.txt")
cat("WGCNA 共表达网络分析总结报告\n")
cat("========================================\n\n")
cat("分析日期: ", date(), "\n")
cat("R 版本: ", R.version.string, "\n")
cat("WGCNA包版本: ", as.character(packageVersion("WGCNA")), "\n\n") ##packageVersion("WGCNA") 返回的是一个 package_version 类对象
cat("------- 数据概览 -------\n")
cat("原始基因数: ", nrow(expr_raw), "\n")
cat("过滤后基因数: ", nrow(expr_final), "\n")
cat("样本数: ", ncol(expr_final), "\n\n")
cat("------- 网络参数 -------\n")
cat("软阈值 (Power): ", softPower, "\n")
cat("网络类型: 无符号 (unsigned)\n")
cat("最小模块大小: ", minModuleSize, "\n")
cat("合并截断高度: 0.25\n\n")
cat("------- 模块检测结果 -------\n")
cat("检测到的模块总数 (含grey): ", length(unique(module_colors)), "\n\n")
for (i in 1:nrow(module_stats)) {
  cat(module_stats$ModuleColor[i], ": ", module_stats$GeneCount[i], " 个基因\n")
}
cat("\n------- 输出文件清单 -------\n")
all_outputs <- list.files(recursive=TRUE, pattern="\\.(csv|txt|png|pdf|RData)$")
for (f in all_outputs) {
  fs <- file.size(f)
  cat(f, sprintf(" (%.1f KB)\n", fs/1024))
}
sink()

cat("\n")
cat("╔══════════════════════════════════════╗\n")
cat("║       分析流程全部完成！             ║\n")
cat("║       请查看 05_Results/ 目录下的文件║\n")
cat("╚══════════════════════════════════════╝\n")
cat("\n下一步建议:\n")
cat("1. 检查 ‘Gene_Module_Assignments_Full.csv’ 查看基因所属模块\n")
cat("2. 使用 ‘Module_Eigengenes_Expression.csv’ 进行模块-性状关联分析\n")
cat("3. 对感兴趣模块的基因进行功能富集分析\n")
cat("\n输出目录结构:\n")
cat("01_InputData/    # 预处理数据\n")
cat("02_QC/           # 质控图表\n")
cat("03_Network/      # 网络对象与TOM矩阵\n")
cat("04_Modules/      # 模块可视化\n")
cat("05_Results/      # 所有结果表格与总结\n")

# ---------- 第八部分：导出Cytoscape网络文件 （修正与优化版）----------
cat("\n[7.4] 正在为 Cytoscape 导出网络文件（内存优化 + kME筛选模式）...\n")

# 创建 Cytoscape 输出目录
cytoscape_dir <- "05_Results/Cytoscape_Networks"
if (!dir.exists(cytoscape_dir)) {
  dir.create(cytoscape_dir, recursive = TRUE)
  cat("   创建目录:", cytoscape_dir, "\n")
}

# 获取当前内存使用情况
get_memory_usage <- function() {
  if (.Platform$OS.type == "windows") {
    mem <- system("wmic OS get FreePhysicalMemory /Value", intern = TRUE)
    mem <- as.numeric(gsub("[^0-9]", "", mem[grepl("FreePhysicalMemory", mem)]))
    return(mem / 1024 / 1024) # 转换为GB
  }
  return(NA)
}

cat("   当前可用内存（估计）:", round(get_memory_usage(), 2), "GB\n")

# 设置网络类型（与网络构建时一致）
networkType <- "unsigned"
cat("   网络类型:", networkType, "\n")

# v1.2.0 修正：使用 module_colors (与第五部分定义一致)
if (!exists("module_colors")) {
  stop("错误：找不到 module_colors 对象。请确保已运行第五部分（网络构建）。")
}

# 获取模块信息（跳过灰色模块）
modules <- unique(module_colors)
modules <- modules[modules != "grey"]

cat("   将处理", length(modules), "个模块的 Cytoscape 导出...\n")
cat("   [内存保护] 每个模块独立计算TOM\n")
cat("   [核心筛选] 对于大模块，优先导出 kME 最高的 Hub 基因\n")

# 定义每个模块导出网络的最大基因数
max_genes_per_module <- 300 

# 逐模块处理
successful_exports <- 0
for (i in seq_along(modules)) {
  mod <- modules[i]
  cat(sprintf("   [%d/%d] 处理模块: %s", i, length(modules), mod))
  
  # 修正：使用 module_colors 获取基因
  modGenes <- colnames(datExpr_forWGCNA)[module_colors == mod]
  n_genes <- length(modGenes)
  cat(sprintf(" (%d 个基因)\n", n_genes))
  
  # 跳过基因数过少的模块
  if (n_genes < 10) {
    cat("       → 跳过（基因数<10，不适合网络可视化）\n")
    next
  }
  
  # --- 核心优化：基于 kME 筛选 Top 基因 ---
  if (n_genes > max_genes_per_module) {
    cat(sprintf("       → 模块较大，基于 kME (Hub基因) 筛选前 %d 个关键基因...\n", max_genes_per_module))
    
    # 1. 匹配当前模块的特征向量名
    ME_name <- paste0("ME", mod)
    
    # 2. 检查 MEs 对象是否存在
    if (exists("MEs") && ME_name %in% colnames(MEs)) {
      # 3. 提取特征向量
      curr_ME <- MEs[, ME_name]
      # 4. 提取表达矩阵
      curr_datExpr <- datExpr_forWGCNA[, modGenes]
      # 5. 计算 kME (绝对值)
      gene_kME <- abs(cor(curr_datExpr, curr_ME, use = "p"))
      # 6. 排序
      kME_ranking <- data.frame(GeneID = rownames(gene_kME), kME_Value = as.vector(gene_kME))
      kME_ranking <- kME_ranking[order(-kME_ranking$kME_Value), ]
      # 7. 截取
      modGenes <- kME_ranking$GeneID[1:max_genes_per_module]
      n_genes <- length(modGenes)
      
      cat(sprintf("       √ 已筛选 kME 最高的 %d 个基因 (kME范围: %.2f - %.2f)\n", 
                  n_genes, max(kME_ranking$kME_Value[1:max_genes_per_module]), 
                  min(kME_ranking$kME_Value[1:max_genes_per_module])))
    } else {
      cat("       ⚠️ 警告: 未找到模块特征向量，回退到随机抽样\n")
      set.seed(123)
      modGenes <- sample(modGenes, max_genes_per_module)
      n_genes <- length(modGenes)
    }
  }
  # ---------------------------------------
  
  tryCatch({
    # 提取该模块的表达数据子集
    modExpr <- datExpr_forWGCNA[, modGenes, drop = FALSE]
    
    # 计算该模块的TOM
    cat("       - 计算模块内TOM相似性...")
    TOM_mod <- TOMsimilarityFromExpr(
      modExpr,
      power = softPower,
      networkType = networkType,
      verbose = 0,
      corType = "pearson", 
      maxPOutliers = 0.05
    )
    cat("完成\n")
    
    # 设置连接阈值 (保留前10%强连接)
    threshold <- quantile(TOM_mod[lower.tri(TOM_mod)], probs = 0.90)
    cat(sprintf("       - 连接阈值: %.3f (Top 10%%)\n", threshold))
    
    # 生成输出文件名
    edge_file <- file.path(cytoscape_dir, sprintf("Cytoscape_%s_edges.txt", mod))
    node_file <- file.path(cytoscape_dir, sprintf("Cytoscape_%s_nodes.txt", mod))
    
    # 导出为Cytoscape格式
    cat("       - 导出文件...")
    exportNetworkToCytoscape(
      TOM_mod,
      edgeFile = edge_file,
      nodeFile = node_file,
      weighted = TRUE,
      threshold = threshold,
      nodeNames = modGenes,
      altNodeNames = modGenes,
      nodeAttr = rep(mod, n_genes)
    )
    cat("完成\n")
    
    successful_exports <- successful_exports + 1
    
    # 内存清理
    rm(TOM_mod, modExpr)
    gc(full = TRUE, verbose = FALSE)
    
  }, error = function(e) {
    cat(sprintf("       → 模块 %s 导出失败: %s\n", mod, e$message))
  })
}

cat("\n   → Cytoscape 导出总结:\n")
cat(sprintf("       成功导出 %d/%d 个模块\n", successful_exports, length(modules)))
cat(sprintf("       输出目录: %s/\n", cytoscape_dir))

# 添加分隔线
cat(paste0("\n", strrep("=", 60), "\n"))
cat("CYTOSCAPE 导出完成！\n")
cat(strrep("=", 60), "\n")

# 显示文件列表
if (dir.exists(cytoscape_dir)) {
  files <- list.files(cytoscape_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) > 0) {
    cat(sprintf("\n共生成 %d 个文件 (显示前5个):\n", length(files)))
    for (f in head(files, 5)) {
      size_kb <- round(file.info(f)$size / 1024, 1)
      cat(sprintf("  %s (%.1f KB)\n", basename(f), size_kb))
    }
  }

}
