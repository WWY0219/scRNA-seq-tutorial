
sc_qc_mouse <- function(seurat_obj = NULL,
                        out_dir = NULL,
                        width = 10,
                        height = 12,
                        dpi = 300) {
    # 加载必要的R包
    suppressMessages({
        library(Seurat)
        library(ggplot2)
        library(dplyr)
        library(cowplot)
        library(patchwork)
        library(RColorBrewer)
    })
    Sys.setenv(LANGUAGE = "en")
    options(stringsAsFactors = FALSE)
    
    # 获取Seurat对象的名称
    obj_name <- deparse(substitute(seurat_obj))
    
    # ---------------------------------------- 00. Pre-Check ---------------------------------------------------------
    if (is.null(seurat_obj)) stop("Please provide a Seurat object (seurat_obj)!")
    if (is.null(out_dir)) stop("Please specify the output directory (out_dir)!")
    
    # ---------------------------------------- 01. Create Output Directory -------------------------------------------
    output_dir <- file.path(out_dir, obj_name)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        message(paste("Output directory automatically created:", output_dir))
    }
        
    # ---------------------------------------- 02. Quality Control (Mouse Specific) ------------------------
    # 1. 计算线粒体基因比例 (老鼠通常以 mt- 开头)
    seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^mt-")

    # 2. 计算核糖体基因比例 (老鼠通常以 Rps 或 Rpl 开头)
    # 使用 ignore.case = FALSE 确保匹配老鼠特有的命名规则
    seurat_obj[["percent.ribo1"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^Rps")
    seurat_obj[["percent.ribo2"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^Rpl")

    # 3. 计算血红蛋白/红细胞比例 (老鼠通常以 Hb 后面接字母开头，如 Hba, Hbb)
    seurat_obj[["percent_RBC"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^Hb[ab]-")
    
    # 如果上面的Hb[ab]-没有匹配到，也可以尝试更通用的 ^Hb (注意不要匹配到 Hbp 等其他基因)
    if(all(seurat_obj[["percent_RBC"]] == 0)) {
        seurat_obj[["percent_RBC"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^Hb")
    }

    # ---------------------------------------- 03. Generate QC Plots ------------------------------------------------
    # (绘图部分逻辑与人类一致)
    
    # RNA count vs. feature count scatter plot
    plot1 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_nRNA.pdf")), 
           plot = plot1, width = width, height = height, dpi = dpi)
    
    # Mitochondrial percentage scatter plot
    plot2 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_MT.pdf")), 
           plot = plot2, width = width, height = height, dpi = dpi)

    # Ribosomal gene scatter plots
    plot3 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo1")
    plot4 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo2")
    qc1 <- plot3 + plot4 
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_RPLS.pdf")), 
           plot = qc1, width = width, height = height, dpi = dpi)

    # Hemoglobin (RBC) percentage scatter plot
    plot5 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent_RBC")
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_RBC.pdf")), 
           plot = plot5, width = width, height = height, dpi = dpi)
    
    # Violin plots for all QC metrics
    qc_meta_violin <- Seurat::VlnPlot(
        seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo1", "percent.ribo2", "percent_RBC"),
        ncol = 3, pt.size = 0)
    
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_metaviolin.pdf")), 
           plot = qc_meta_violin, width = width, height = height, dpi = dpi)

    # ---------------------------------------- 04. Metadata Preview -----------------------------------------------
    cat("\n=== First 5 rows of metadata (Mouse QC metrics) ===\n")
    print(head(seurat_obj@meta.data, 5))
    
    # ---------------------------------------- 05. Return Updated Seurat Object -----------------------------------
    return(seurat_obj)
}
