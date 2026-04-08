#' DoubleFinder-based Doublet Removal for Mouse Samples
#'
#' @title Doublet Removal with Preprocessing (Mouse Optimized)
#' @description 
#' 专门针对老鼠单细胞数据优化的双细胞去除流程。
#' 步骤：1. 创建输出目录；2. 数据预处理；3. 降维聚类；4. 寻找最佳 pK 值；
#' 5. 运行 DoubletFinder (高/低置信度)；6. 过滤并保存结果。
#'
#' @param seurat_obj 输入的老鼠 Seurat 对象。
#' @param max.dim 用于 UMAP 和聚类的 PCA 维度数 (建议 20-30)。
#' @param max.pcs DoubletFinder 模拟时使用的 PCA 维度数 (建议与 max.dim 一致)。
#' @param res 聚类分辨率 (默认 1.0)。
#' @param dbrate 预估的双细胞比率。10x 官方通常为 每 1000 个细胞 0.8% (8 * 1e-6)。
#' @param pN 模拟人工双细胞的比例 (默认 0.25)。
#' @param out_dir 输出路径。
#' @author WWY (Mouse Version Adaptation)
#' @export

sc_doublefinder_mouse <- function(seurat_obj = NULL,
                                  max.dim = 30, 
                                  max.pcs = 30,
                                  res = 1.0,
                                  dbrate = 8 * 1e-6,
                                  pN = 0.25,
                                  width = 8, 
                                  height = 6,
                                  out_dir = NULL) {
    ## 1. 加载库
    suppressMessages({
      library(DoubletFinder)
      library(Seurat)
      library(ggplot2)
      library(dplyr)
      library(qs)
    })
    set.seed(1234)

    # 备份原始对象用于最后提取（避免预处理污染）
    seurat_raw_backup <- seurat_obj
    
    ###-------------------------------- 0. 参数检查 -------------------------------------------###
    if (is.null(seurat_obj)) stop("错误: 必须提供 'seurat_obj'!")
    if (is.null(out_dir)) stop("错误: 必须提供 'out_dir'!")
    
    ###-------------------------------- 1. 创建输出目录 ---------------------------------------###
    # 尝试获取对象名
    obj_name <- deparse(substitute(seurat_obj))
    # 如果通过 pipe 或列表传入导致名字失效，给定默认名
    if(length(obj_name) > 1 || obj_name == "." ) obj_name <- "Mouse_Sample"
    
    dir_name <- sub("_[^_]*$", "", obj_name)
    output_dir <- file.path(out_dir, paste0("DoubletFinder_", dir_name))
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        message("已自动创建输出目录: ", output_dir)
    }

    ###-------------------------------- 2. 预处理 (Mouse Data) -------------------------------###
    message(">>> 正在进行预处理 (Normalization, PCA)...")
    # 老鼠数据的 scale.factor 建议使用常用值 10000 或中位数
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000) 
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)

    ###-------------------------------- 3. 降维与聚类 -----------------------------------------###
    message(">>> 正在运行 UMAP 降维与聚类...")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:max.dim, reduction.name = "umap_naive", verbose = FALSE)  
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:max.dim, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    
    ###-------------------------------- 4. DoubletFinder 核心计算 ----------------------------###
    # A. 寻找最佳 pK
    message(">>> 正在计算老鼠样本的最佳 pK 值...")
    sweep.res.list <- paramSweep(seurat_obj, PCs = 1:max.pcs, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    opt_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message("确定最佳 pK 为: ", opt_pK)
    
    # B. 估计同质双细胞比例 (Homotypic Doublets)
    homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)
    # 根据细胞总数计算预估双细胞数
    nExp_poi <- round(dbrate * (ncol(seurat_obj)^2 / 1000)) # 修正后的通用 10x 算法
    # 如果你的 dbrate 是直接比例 (如 0.08)，则用：nExp_poi <- round(dbrate * ncol(seurat_obj))
    
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    message("预估双细胞数: ", nExp_poi, " | 调整后 (去除同质): ", nExp_poi.adj)
  
    # C. 运行 DoubletFinder
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:max.pcs, pN = pN, pK = opt_pK, nExp = nExp_poi.adj, sct = FALSE)
    
    # 获取高置信度列名
    df_col_high <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
    pANN_col <- grep("pANN", colnames(seurat_obj@meta.data), value = TRUE)

    ###-------------------------------- 5. 可视化 ---------------------------------------------###
    p1 <- DimPlot(seurat_obj, reduction = "umap_naive", group.by = df_col_high, cols = c("Singlet"="skyblue", "Doublet"="red")) +
          ggtitle(paste0("Mouse Doublets (", dir_name, ")"))
    
    ggsave(filename = file.path(output_dir, paste0("Doublets_Plot_", dir_name, ".pdf")), 
           plot = p1, width = width, height = height)
    
    ###-------------------------------- 6. 提取 Singlets 并过滤 -------------------------------###
    # 提取 Singlet 细胞的 ID
    singlet_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[df_col_high]] == "Singlet"]
    removed_ids <- setdiff(colnames(seurat_obj), singlet_cells)
    
    # 保存被剔除的 ID
    write.table(data.frame(Doublet_IDs = removed_ids),
                file.path(output_dir, paste0(dir_name, "_removed_doublets.txt")),
                row.names = FALSE, quote = FALSE)

    # 从原始备份中提取干净的细胞（保持原始 Count 矩阵不被 PCA/ScaleData 修改过）
    seurat_clean <- subset(seurat_raw_backup, cells = singlet_cells)
    
    ###-------------------------------- 7. 保存并返回 -----------------------------------------###
    qsave(seurat_clean, file = file.path(output_dir, paste0(dir_name, "_clean_mouse.qs")))
    
    message("\n===== 完成! =====")
    message("原始细胞数: ", ncol(seurat_obj))
    message("保留细胞数 (Singlets): ", length(singlet_cells))
    message("剔除双细胞数: ", length(removed_ids))
    
    return(seurat_clean)
}
