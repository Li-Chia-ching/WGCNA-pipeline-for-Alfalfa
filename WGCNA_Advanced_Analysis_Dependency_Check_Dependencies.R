# ==============================================================================
# WGCNA v2.4 环境诊断与修复工具
# 功能：一键检测、安装依赖、验证环境是否满足 16GB 优化版需求
# ==============================================================================

cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║      WGCNA v2.4 环境诊断与修复工具 (终极版)       ║\n")
cat("╚═══════════════════════════════════════════════════╝\n\n")

# 1. 基础环境检测
# ==============================================================================
cat("[Step 1] 系统环境检测...\n")
r_ver <- R.Version()
cat(paste("   R 版本:", r_ver$major, ".", r_ver$minor, "\n"))
cat(paste("   操作系统:", Sys.info()["sysname"], "\n"))

# 内存检测 (粗略估计)
if (.Platform$OS.type == "windows") {
  mem_free <- tryCatch({
    wmic_out <- system("wmic OS get FreePhysicalMemory /Value", intern = TRUE)
    mem_kb <- as.numeric(gsub("[^0-9]", "", wmic_out[grepl("FreePhysicalMemory", wmic_out)]))
    round(mem_kb / 1024, 2)
  }, error = function(e) NA)
  
  if (!is.na(mem_free)) {
    cat(paste("   当前可用内存:", mem_free, "MB\n"))
    if (mem_free < 4000) cat("   ⚠️ 警告: 可用内存低于 4GB，建议关闭其他软件。\n")
  }
}

# 2. 网络连接与镜像设置
# ==============================================================================
cat("\n[Step 2] 网络连接检测...\n")
# 设置超时
options(timeout = 600)

# 优先使用清华源，速度快且稳定
cat("   正在配置清华大学 CRAN/Bioc 镜像源...\n")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioc/")

# 简单测试连接
tryCatch({
  readLines("https://mirrors.tuna.tsinghua.edu.cn/CRAN/", n = 1)
  cat("   ✅ 网络连接正常 (清华源)\n")
}, error = function(e) {
  cat("   ⚠️ 连接清华源失败，尝试切换回官方源...\n")
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  options(BioC_mirror = "https://bioconductor.org")
})

# 3. 定义依赖包清单 (精准匹配 v2.4)
# ==============================================================================
# 基础绘图与数据处理
cran_pkgs <- c(
  "WGCNA", "ggplot2", "dplyr", "RColorBrewer", "viridis", 
  "pheatmap", "igraph", "patchwork", "ggrepel", "data.table" # data.table 用于加速读取
)

# 生物信息学专用包
bioc_pkgs <- c(
  "BiocManager", "clusterProfiler", "enrichplot", "org.Hs.eg.db", 
  "GOSemSim", "AnnotationDbi", "impute", "preprocessCore", "GO.db"
)

# 4. 智能安装函数
# ==============================================================================
install_if_missing <- function(pkg, type = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("   📦 正在安装: ", pkg, " (", type, ")...\n"))
    
    tryCatch({
      if (type == "CRAN") {
        install.packages(pkg, quiet = TRUE)
      } else {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
      }
      
      if (requireNamespace(pkg, quietly = TRUE)) {
        cat(paste0("     ✅ ", pkg, " 安装成功\n"))
        return(TRUE)
      } else {
        cat(paste0("     ❌ ", pkg, " 安装失败\n"))
        return(FALSE)
      }
    }, error = function(e) {
      cat(paste0("     ❌ 错误: ", e$message, "\n"))
      return(FALSE)
    })
  } else {
    cat(paste0("   ✅ ", pkg, " 已安装\n"))
    return(TRUE)
  }
}

# 5. 执行安装流程
# ==============================================================================
cat("\n[Step 3] 检查 CRAN 依赖包...\n")
cran_results <- sapply(cran_pkgs, install_if_missing, type = "CRAN")

cat("\n[Step 4] 检查 Bioconductor 依赖包...\n")
# 预先检查 BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 修复常见的 clusterProfiler 依赖问题
if (!requireNamespace("GOSemSim", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE)) {
  cat("   ℹ️ 预先安装 clusterProfiler 核心依赖...\n")
  BiocManager::install(c("GOSemSim", "GO.db", "AnnotationDbi"), update = FALSE, ask = FALSE)
}

bioc_results <- sapply(bioc_pkgs, install_if_missing, type = "Bioc")

# 6. WGCNA 专用修复 (解决 preprocessCore 问题)
# ==============================================================================
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  cat("\n[Step 5] 尝试修复 WGCNA...\n")
  BiocManager::install(c("impute", "preprocessCore", "GO.db"), update = FALSE, ask = FALSE)
  install.packages("WGCNA")
}

# 7. 最终诊断报告
# ==============================================================================
cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║               环境诊断报告                        ║\n")
cat("╚═══════════════════════════════════════════════════╝\n")

all_pkgs <- c(cran_pkgs, bioc_pkgs)
missing <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing) == 0) {
  cat("\n🎉 恭喜！所有依赖包状态完美。\n")
  cat("🚀 您可以直接运行 WGCNA v2.4 脚本了。\n\n")
  
  # 简单的加载测试
  cat("🔍 正在进行加载测试...\n")
  suppressPackageStartupMessages({
    library(WGCNA)
    library(clusterProfiler)
    library(pheatmap)
  })
  cat("✅ 核心包加载成功！\n")
  
} else {
  cat("\n⚠️  以下包仍然缺失或安装失败：\n")
  cat(paste("   ❌", missing, collapse = "\n"), "\n")
  
  cat("\n🔧 建议手动修复命令：\n")
  cat("--------------------------------------------------\n")
  cat('# CRAN 包手动安装:\n')
  cran_miss <- intersect(missing, cran_pkgs)
  if(length(cran_miss) > 0) 
    cat(sprintf('install.packages(c("%s"))\n', paste(cran_miss, collapse = '", "')))
  
  cat('\n# Bioconductor 包手动安装:\n')
  bioc_miss <- intersect(missing, bioc_pkgs)
  if(length(bioc_miss) > 0)
    cat(sprintf('BiocManager::install(c("%s"), force = TRUE)\n', paste(bioc_miss, collapse = '", "')))
  cat("--------------------------------------------------\n")
  
  if ("clusterProfiler" %in% missing) {
    cat("\n💡 提示: clusterProfiler 安装失败通常是因为网络原因中断。\n")
    cat("   请尝试重启 RStudio 后单独运行: BiocManager::install('clusterProfiler')\n")
  }
}
