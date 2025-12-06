# =======================================================
# 现代化依赖包检查与安装程序（终极版）
# 包含网络检测、自动重试、包修复和重启提示
# =======================================================

cat("🔧 WGCNA依赖包智能安装与修复程序\n")
cat("版本: 3.0 | 日期: 2024\n")
cat("==================================================\n\n")

# 0. 检查R版本和重启状态
cat("0. 检查系统环境...\n")
r_version <- R.Version()
cat(paste("   R版本:", r_version$major, ".", r_version$minor, "\n"))
cat(paste("   平台:", r_version$platform, "\n"))
cat(paste("   库路径: ", paste(.libPaths(), collapse = ", "), "\n\n"))

# 1. 首先测试网络连接
test_internet_connection <- function() {
  cat("1. 测试网络连接...\n")
  
  test_urls <- c(
    "https://cran.r-project.org",
    "https://bioconductor.org",
    "https://mirrors.tuna.tsinghua.edu.cn"
  )
  
  has_connection <- FALSE
  for(url in test_urls) {
    tryCatch({
      test <- suppressWarnings(
        readLines(url, n = 1, warn = FALSE)
      )
      cat(paste("   ✅ 可访问:", url, "\n"))
      has_connection <- TRUE
      break
    }, error = function(e) {
      cat(paste("   ❌ 无法访问:", url, "\n"))
    })
  }
  
  return(has_connection)
}

# 2. 智能选择镜像源
setup_mirrors <- function(has_connection) {
  cat("\n2. 配置镜像源...\n")
  
  if(has_connection) {
    # 尝试多个镜像，选择最快的一个
    mirrors <- list(
      tuna = list(
        CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
        Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioc/"
      ),
      ustc = list(
        CRAN = "https://mirrors.ustc.edu.cn/CRAN/",
        Bioc = "https://mirrors.ustc.edu.cn/bioc/"
      ),
      aliyun = list(
        CRAN = "https://mirrors.aliyun.com/CRAN/",
        Bioc = "https://mirrors.aliyun.com/bioc/"
      ),
      official = list(
        CRAN = "https://cloud.r-project.org",
        Bioc = "https://bioconductor.org"
      )
    )
    
    # 测试镜像速度
    test_mirror_speed <- function(mirror_url) {
      tryCatch({
        start <- Sys.time()
        test <- suppressWarnings(
          readLines(paste0(mirror_url, "web/checks/index.html"), n = 1, warn = FALSE)
        )
        end <- Sys.time()
        return(as.numeric(difftime(end, start, units = "secs")))
      }, error = function(e) {
        return(Inf)
      })
    }
    
    cat("   测试镜像源速度...\n")
    best_mirror <- "official"
    best_time <- Inf
    
    for(mirror_name in names(mirrors)) {
      time_taken <- test_mirror_speed(mirrors[[mirror_name]]$CRAN)
      if(time_taken < best_time) {
        best_time <- time_taken
        best_mirror <- mirror_name
      }
      cat(paste0("     ", mirror_name, ": ", 
                 ifelse(is.infinite(time_taken), "不可用", 
                        paste(round(time_taken, 2), "秒")), "\n"))
    }
    
    selected_mirror <- mirrors[[best_mirror]]
    cat(paste0("\n   ✅ 选择镜像源: ", best_mirror, "\n"))
    
  } else {
    cat("   ⚠️ 无网络连接，使用默认镜像\n")
    selected_mirror <- list(
      CRAN = "https://cloud.r-project.org",
      Bioc = "https://bioconductor.org"
    )
  }
  
  # 设置镜像
  options(repos = c(CRAN = selected_mirror$CRAN))
  options(BioC_mirror = selected_mirror$Bioc)
  
  return(selected_mirror)
}

# 3. 检查是否安装BiocManager
check_biocmanager <- function() {
  cat("\n3. 检查Bioconductor管理器...\n")
  
  if(!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("   正在安装BiocManager...\n")
    tryCatch({
      install.packages("BiocManager", quiet = TRUE)
      cat("   ✅ BiocManager安装成功\n")
    }, error = function(e) {
      cat("   ❌ BiocManager安装失败\n")
      cat("   请手动运行: install.packages('BiocManager')\n")
      return(FALSE)
    })
  } else {
    cat("   ✅ BiocManager已安装\n")
  }
  
  return(TRUE)
}

# 4. 增强型智能安装包（带重试机制和验证）
smart_install_package <- function(pkg, type = "CRAN", max_retries = 2) {
  cat(paste("   📦", pkg, paste0("(", type, "): ")))
  
  # 检查是否已安装且可加载
  if(check_package_status(pkg)) {
    cat("✅ 已安装且可加载\n")
    return(TRUE)
  }
  
  # 安装函数
  install_func <- if(type == "CRAN") {
    function() install.packages(pkg, quiet = TRUE, dependencies = TRUE)
  } else if(type == "Bioc") {
    function() BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
  }
  
  # 尝试安装（带重试）
  success <- FALSE
  for(attempt in 1:max_retries) {
    tryCatch({
      if(attempt > 1) {
        cat(paste0("    重试 (", attempt, "/", max_retries, ")... "))
      }
      
      suppressWarnings(install_func())
      
      # 验证安装
      if(check_package_status(pkg)) {
        success <- TRUE
        break
      }
    }, error = function(e) {
      if(attempt == max_retries) {
        # 最后一次尝试失败
        cat("❌ 安装失败\n")
      }
    })
  }
  
  if(success) {
    cat("✅ 安装成功\n")
  } else {
    cat("⚠️  安装失败，需要手动安装\n")
  }
  
  return(success)
}

# 5. 检查包状态（安装和加载）
check_package_status <- function(pkg) {
  # 检查是否已安装
  if(!requireNamespace(pkg, quietly = TRUE)) {
    return(FALSE)
  }
  
  # 尝试加载，捕获错误
  tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# 6. GOSemSim包专用修复函数
repair_gosemsim <- function() {
  cat("\n🔄 特别修复GOSemSim包...\n")
  cat("   ==================================\n")
  
  # 检查当前状态
  cat("   检查当前状态...\n")
  if(check_package_status("GOSemSim")) {
    cat("   ✅ GOSemSim已安装且可加载\n")
    return(TRUE)
  }
  
  # 方法1: 强制重新安装
  cat("   方法1: 强制重新安装GOSemSim...\n")
  tryCatch({
    # 先尝试卸载
    try(remove.packages("GOSemSim"), silent = TRUE)
    
    # 重新安装
    if(requireNamespace("BiocManager", quietly = TRUE)) {
      BiocManager::install("GOSemSim", update = FALSE, ask = FALSE, force = TRUE)
      
      # 验证
      if(check_package_status("GOSemSim")) {
        cat("   ✅ GOSemSim修复成功\n")
        return(TRUE)
      }
    }
  }, error = function(e) {
    cat(paste("   ❌ 方法1失败:", e$message, "\n"))
  })
  
  # 方法2: 从源码安装
  cat("   方法2: 从源码安装GOSemSim...\n")
  tryCatch({
    if(requireNamespace("BiocManager", quietly = TRUE)) {
      BiocManager::install("GOSemSim", type = "source", update = FALSE, ask = FALSE)
      
      if(check_package_status("GOSemSim")) {
        cat("   ✅ GOSemSim源码安装成功\n")
        return(TRUE)
      }
    }
  }, error = function(e) {
    cat(paste("   ❌ 方法2失败:", e$message, "\n"))
  })
  
  # 方法3: 检查并修复库路径
  cat("   方法3: 检查库路径...\n")
  lib_paths <- .libPaths()
  cat(paste("   当前库路径:\n"))
  for(i in seq_along(lib_paths)) {
    cat(paste("     ", i, ":", lib_paths[i], "\n"))
  }
  
  # 检查GOSemSim是否在正确的位置
  gosemsim_path <- find.package("GOSemSim", quiet = TRUE)
  if(length(gosemsim_path) > 0) {
    cat(paste("   GOSemSim安装位置:", gosemsim_path, "\n"))
  } else {
    cat("   GOSemSim未在任何库路径中找到\n")
  }
  
  return(FALSE)
}

# 7. clusterProfiler专用修复函数
repair_clusterprofiler <- function() {
  cat("\n🔄 特别修复clusterProfiler包...\n")
  cat("   ==================================\n")
  
  # 检查当前状态
  cat("   检查当前状态...\n")
  if(check_package_status("clusterProfiler")) {
    cat("   ✅ clusterProfiler已安装且可加载\n")
    return(TRUE)
  }
  
  # 先修复GOSemSim
  if(!check_package_status("GOSemSim")) {
    cat("   ℹ️ 需要先修复GOSemSim...\n")
    repair_gosemsim()
  }
  
  # 安装其他依赖
  cat("   安装clusterProfiler依赖...\n")
  dependencies <- c("AnnotationDbi", "IRanges", "BiocGenerics", "S4Vectors")
  for(pkg in dependencies) {
    if(!check_package_status(pkg)) {
      cat(paste("     安装", pkg, "...\n"))
      if(requireNamespace("BiocManager", quietly = TRUE)) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      }
    }
  }
  
  # 重新安装clusterProfiler
  cat("   重新安装clusterProfiler...\n")
  tryCatch({
    if(requireNamespace("BiocManager", quietly = TRUE)) {
      BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE, force = TRUE)
      
      if(check_package_status("clusterProfiler")) {
        cat("   ✅ clusterProfiler修复成功\n")
        return(TRUE)
      }
    }
  }, error = function(e) {
    cat(paste("   ❌ clusterProfiler安装失败:", e$message, "\n"))
  })
  
  return(FALSE)
}

# 8. 手动安装提示生成器
generate_manual_instructions <- function(failed_packages) {
  if(length(failed_packages) == 0) return()
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("⚠️  以下包安装失败，请手动安装：\n")
  cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
  
  for(i in seq_along(failed_packages)) {
    pkg <- failed_packages[[i]]
    if(pkg$type == "CRAN") {
      cat(paste0("📦 ", pkg$name, " (CRAN):\n"))
      cat(paste0("   安装命令: install.packages(\"", pkg$name, "\")\n"))
      cat(paste0("   主页: https://cran.r-project.org/package=", pkg$name, "\n"))
    } else if(pkg$type == "Bioc") {
      cat(paste0("🧬 ", pkg$name, " (Bioconductor):\n"))
      cat(paste0("   安装命令: BiocManager::install(\"", pkg$name, "\")\n"))
      cat(paste0("   主页: https://bioconductor.org/packages/", pkg$name, "\n"))
    }
    cat("\n")
  }
  
  cat("💡 建议的安装顺序：\n")
  cat("1. 确保R版本 >= 4.0.0\n")
  
  cran_failed <- sapply(failed_packages, function(x) x$type == "CRAN")
  bioc_failed <- sapply(failed_packages, function(x) x$type == "Bioc")
  
  if(any(cran_failed)) {
    cran_names <- sapply(failed_packages[cran_failed], function(x) x$name)
    cat(paste0("2. 先安装CRAN包: ", paste(cran_names, collapse = ", "), "\n"))
  }
  
  if(any(bioc_failed)) {
    bioc_names <- sapply(failed_packages[bioc_failed], function(x) x$name)
    cat(paste0("3. 再安装Bioconductor包: ", paste(bioc_names, collapse = ", "), "\n"))
  }
  
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
}

# 9. 主函数：执行完整的依赖检查流程
main_dependency_check <- function() {
  # 测试网络连接
  has_connection <- test_internet_connection()
  
  if(!has_connection) {
    cat("\n⚠️  警告: 未检测到网络连接！\n")
    cat("  将使用已安装的包进行分析，部分功能可能受限。\n")
    cat("  建议连接网络后重新运行此脚本。\n")
  }
  
  # 设置镜像
  current_mirror <- setup_mirrors(has_connection)
  
  # 检查BiocManager
  if(!check_biocmanager()) {
    cat("⚠️  BiocManager检查失败，跳过Bioconductor包安装\n")
    can_install_bioc <- FALSE
  } else {
    can_install_bioc <- TRUE
    # 设置BiocManager选项
    if(has_connection) {
      tryCatch({
        options(BioC_mirror = current_mirror$Bioc)
      }, error = function(e) NULL)
    }
  }
  
  # 包列表
  cran_packages <- c(
    "ggplot2", "pheatmap", "igraph", "dplyr", "tidyverse",
    "RColorBrewer", "viridis", "ggsci", "ggrepel", 
    "patchwork", "scales", "stringr", "tidyr", "readr",
    "reshape2", "gridExtra", "cowplot", "ggpubr",
    "WGCNA", "flashClust"
  )
  
  bioc_packages <- c(
    "clusterProfiler", "enrichplot", "org.Hs.eg.db",
    "DOSE", "ggplotify", "ggnewscale", "GOSemSim",
    "AnnotationDbi", "topGO", "pathview"
  )
  
  # 移除重复的包（如果有）
  cran_packages <- unique(cran_packages)
  bioc_packages <- unique(bioc_packages)
  
  # 检查已安装的包
  cat("\n4. 检查已安装的包...\n")
  cat(paste(rep("-", 40), collapse = ""), "\n", sep = "")
  
  # 安装CRAN包
  cat("\n5. 安装CRAN包...\n")
  cat(paste(rep("-", 40), collapse = ""), "\n", sep = "")
  
  cran_results <- list()
  for(pkg in cran_packages) {
    # 跳过已安装的包
    if(check_package_status(pkg)) {
      cat(paste("   📦", pkg, "(CRAN): ✅ 已安装且可加载\n"))
      cran_results[[pkg]] <- TRUE
    } else {
      success <- smart_install_package(pkg, type = "CRAN")
      cran_results[[pkg]] <- success
    }
  }
  
  # 安装Bioconductor包
  if(can_install_bioc) {
    cat("\n6. 安装Bioconductor包...\n")
    cat(paste(rep("-", 40), collapse = ""), "\n", sep = "")
    
    bioc_results <- list()
    for(pkg in bioc_packages) {
      # 跳过已安装的包
      if(check_package_status(pkg)) {
        cat(paste("   🧬", pkg, "(Bioc): ✅ 已安装且可加载\n"))
        bioc_results[[pkg]] <- TRUE
      } else {
        success <- smart_install_package(pkg, type = "Bioc")
        bioc_results[[pkg]] <- success
      }
    }
  } else {
    bioc_results <- list()
    for(pkg in bioc_packages) {
      bioc_results[[pkg]] <- FALSE
    }
    cat("\n⚠️  跳过Bioconductor包安装（需要先安装BiocManager）\n")
  }
  
  # 7. 特别修复关键包
  cat("\n7. 特别修复关键包...\n")
  cat(paste(rep("-", 40), collapse = ""), "\n", sep = "")
  
  # 检查并修复GOSemSim
  if("GOSemSim" %in% bioc_packages && !check_package_status("GOSemSim")) {
    repair_gosemsim()
  }
  
  # 检查并修复clusterProfiler
  if("clusterProfiler" %in% bioc_packages && !check_package_status("clusterProfiler")) {
    repair_clusterprofiler()
  }
  
  # 汇总结果
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("📊 安装结果汇总\n")
  cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
  
  # 创建结果列表
  all_packages <- list()
  
  for(pkg in cran_packages) {
    all_packages[[length(all_packages) + 1]] <- list(
      name = pkg, 
      type = "CRAN", 
      success = cran_results[[pkg]]
    )
  }
  
  for(pkg in bioc_packages) {
    all_packages[[length(all_packages) + 1]] <- list(
      name = pkg, 
      type = "Bioc", 
      success = bioc_results[[pkg]]
    )
  }
  
  # 统计
  total <- length(all_packages)
  installed <- sum(sapply(all_packages, function(x) x$success))
  failed <- total - installed
  
  cat(paste("📈 总计:", total, "个包\n"))
  cat(paste("✅ 成功:", installed, "个\n"))
  cat(paste("❌ 失败:", failed, "个\n"))
  
  # 特别检查关键包
  cat("\n🔍 关键包检查:\n")
  critical_packages <- c("WGCNA", "clusterProfiler", "ggplot2", "pheatmap", "igraph")
  for(pkg in critical_packages) {
    status <- ifelse(check_package_status(pkg), "✅ 已安装且可加载", "❌ 未安装或无法加载")
    cat(paste("   ", pkg, ":", status, "\n"))
  }
  
  # 检查clusterProfiler依赖
  if(check_package_status("clusterProfiler") && !check_package_status("GOSemSim")) {
    cat("\n⚠️  警告: clusterProfiler已安装但GOSemSim未安装，这可能导致clusterProfiler加载失败\n")
  }
  
  if(failed > 0) {
    # 获取失败的包
    failed_pkgs <- all_packages[sapply(all_packages, function(x) !x$success)]
    
    # 生成手动安装提示
    generate_manual_instructions(failed_pkgs)
    
    cat("\n💡 临时解决方案：\n")
    cat("   1. 可以直接运行主脚本，脚本会自动跳过缺失的包\n")
    cat("   2. 或手动安装失败包后重新运行此脚本\n")
    
    # 检查是否为clusterProfiler依赖问题
    if("clusterProfiler" %in% sapply(failed_pkgs, function(x) x$name)) {
      cat("\n⚠️  注意：clusterProfiler安装失败可能是由于缺少依赖包\n")
      cat("   尝试运行以下命令修复：\n")
      cat('   BiocManager::install(c("GOSemSim", "AnnotationDbi", "IRanges", "BiocGenerics"))\n')
      cat('   然后重新安装：BiocManager::install("clusterProfiler")\n')
    }
    
    # 询问用户是否继续
    if(interactive()) {
      cat("\n❓ 是否继续运行主分析脚本？(y/n): ")
      answer <- readline()
      
      if(tolower(answer) %in% c("y", "yes", "是")) {
        cat("\n🚀 继续运行主脚本...\n")
        return(TRUE)
      } else {
        cat("\n⏸️  请先手动安装缺失的包，然后重新运行\n")
        return(FALSE)
      }
    } else {
      cat("\n⏸️  检测到安装失败，建议先手动安装缺失的包\n")
      return(FALSE)
    }
  } else {
    cat("\n🎉 所有依赖包安装成功！\n")
    
    # 测试加载关键包
    cat("\n🧪 测试加载关键包...\n")
    test_packages <- c("WGCNA", "clusterProfiler", "ggplot2", "igraph", "pheatmap")
    for(pkg in test_packages) {
      if(check_package_status(pkg)) {
        version <- packageVersion(pkg)
        cat(paste("   ", pkg, "版本", version, ": ✅ 可正常加载\n"))
      } else {
        cat(paste("   ", pkg, ": ⚠️  安装但无法加载\n"))
      }
    }
    
    cat("\n现在可以运行主分析脚本了。\n")
    return(TRUE)
  }
}

# 10. 优雅的错误处理
cat("开始执行依赖包检查与安装...\n\n")

# 检查是否在交互式会话中
if(interactive()) {
  tryCatch({
    should_continue <- main_dependency_check()
    
    if(should_continue) {
      cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
      cat("📝 下一步操作：\n")
      cat("   1. 运行主分析脚本: source('WGCNA_Advanced_Analysis_Modern.R')\n")
      cat("   2. 如果仍有问题，请重启R会话后再次运行此脚本\n")
      cat("   3. 或手动运行以下测试命令验证安装:\n")
      cat('      library(clusterProfiler)\n')
      cat('      library(org.Hs.eg.db)\n')
      cat('      library(WGCNA)\n')
      cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
    }
    
  }, error = function(e) {
    cat("\n❌ 脚本执行过程中出现错误:\n")
    cat(paste("   错误信息:", e$message, "\n"))
    cat("\n💡 解决方案：\n")
    cat("   1. 检查网络连接\n")
    cat("   2. 确保R版本 >= 4.0.0\n")
    cat("   3. 尝试手动安装缺失的包\n")
    cat("   4. 或直接运行主脚本，跳过缺失的包\n")
  })
} else {
  # 非交互式环境，自动运行
  cat("⏱️  检测到非交互式环境，自动运行安装程序...\n")
  tryCatch({
    main_dependency_check()
  }, error = function(e) {
    cat("自动安装过程中出现错误，请手动检查。\n")
  })
}

# 11. 提供简化版安装命令（备用）
cat("\n📋 备用安装命令（如果上述脚本失败）：\n")
cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
cat("# 1. 基础CRAN包\n")
cat('install.packages(c("ggplot2", "pheatmap", "igraph", "dplyr", "tidyverse",\n')
cat('                  "RColorBrewer", "viridis", "ggsci", "ggrepel",\n')
cat('                  "patchwork", "scales", "stringr", "tidyr", "readr",\n')
cat('                  "reshape2", "gridExtra", "cowplot", "ggpubr",\n')
cat('                  "WGCNA", "flashClust"))\n\n')
cat("# 2. Bioconductor包（需要先安装BiocManager）\n")
cat('if (!require("BiocManager", quietly = TRUE))\n')
cat('    install.packages("BiocManager")\n')
cat('# 先安装依赖包\n')
cat('BiocManager::install(c("GOSemSim", "AnnotationDbi", "IRanges", "BiocGenerics"))\n')
cat('# 再安装主要包\n')
cat('BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db",\n')
cat('                      "DOSE", "ggplotify", "ggnewscale", "topGO", "pathview"))\n')
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

# 12. 提供快速修复命令和重启提示
cat("\n🚑 快速修复命令（针对clusterProfiler依赖问题）：\n")
cat(paste(rep("-", 60), collapse = ""), "\n", sep = "")
cat('# 如果clusterProfiler加载失败，运行以下命令：\n')
cat('if (!requireNamespace("BiocManager", quietly = TRUE))\n')
cat('    install.packages("BiocManager")\n')
cat('# 先强制删除有问题的包\n')
cat('try(remove.packages("GOSemSim"), silent = TRUE)\n')
cat('try(remove.packages("clusterProfiler"), silent = TRUE)\n')
cat('# 重新安装\n')
cat('BiocManager::install("GOSemSim", force = TRUE)\n')
cat('BiocManager::install("clusterProfiler", force = TRUE)\n')
cat('# 然后重启R会话\n')
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

# 13. 最后的重启提示
cat("\n💡 重要提示：\n")
cat("   如果安装成功但包仍然无法加载，请重启R会话（Session -> Restart R）\n")
cat("   然后重新运行此脚本或直接运行主分析脚本。\n")
cat("   重启R会话可以解决大部分包加载问题。\n")
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

# 14. 保存会话信息
cat("\n📝 保存安装日志...\n")
tryCatch({
  sink("dependency_install_log.txt")
  cat("依赖包安装日志\n")
  cat("安装时间:", Sys.time(), "\n")
  cat("R版本:", R.version.string, "\n")
  cat("平台:", R.version$platform, "\n\n")
  
  # 检查关键包状态
  cat("关键包状态:\n")
  key_packages <- c("WGCNA", "clusterProfiler", "ggplot2", "igraph", 
                    "pheatmap", "GOSemSim", "enrichplot")
  
  for(pkg in key_packages) {
    if(requireNamespace(pkg, quietly = TRUE)) {
      version <- packageVersion(pkg)
      cat(paste("  ", pkg, ":", version, "\n"))
    } else {
      cat(paste("  ", pkg, ": 未安装或无法加载\n"))
    }
  }
  
  sink()
  cat("✅ 日志已保存到: dependency_install_log.txt\n")
}, error = function(e) {
  cat("⚠️  无法保存日志文件\n")
})
