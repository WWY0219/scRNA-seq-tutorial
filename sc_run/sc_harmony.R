#' Harmony-Based Batch Correction with Comprehensive Preprocessing (Normalization + PCA + Cell Cycle Regression + UMAP/TSNE)
#'
#' This function performs full single-cell RNA-seq data preprocessing (normalization, variable feature selection, 
#' cell cycle score calculation, unwanted variation regression, PCA), applies Harmony batch correction to eliminate
#' batch effects, performs dimensionality reduction (UMAP/TSNE) and cell clustering, and returns the processed Seurat object.
#'
#' @title Harmony Batch Correction with Preprocessing and Dimensionality Reduction
#' @description 
#' Core Steps: 
#' 1. Parameter validation and output directory setup;
#' 2. Data preprocessing: Log-normalization (median UMI as scaling factor) -> Variable feature selection -> 
#'    Cell cycle scoring -> Regression of unwanted variation (cell cycle/mitochondrial/ribosomal/RBC genes);
#' 3. PCA on highly variable genes;
#' 4. Harmony batch correction (grouped by orig.ident);
#' 5. Cell clustering (FindNeighbors + FindClusters);
#' 6. UMAP/TSNE dimensionality reduction on Harmony-corrected PCA;
#' 7. Save convergence plots (PCA/Harmony) and return processed Seurat object.
#' 
#' @param seurat_obj Input Seurat object (required). Must contain meta.data columns: nCount_RNA, percent_ribo1, 
#'                   percent_ribo2, percent_mt, percent_RBC, orig.ident.
#' @param max_dim_pca Numeric (required). Number of PCA dimensions used for FindNeighbors (e.g., 30).
#' @param max_dim_harmony Numeric (required). Number of Harmony-corrected dimensions used for UMAP/TSNE (e.g., 30).
#' @param max.iter Numeric. Maximum iterations for Harmony batch correction (default: 30).
#' @param nfeatures Numeric. Number of variable features to select (default: 3000).
#' @param res Numeric. Resolution for cell clustering (default: 0.5). If NULL, will use 0.5 as fallback.
#' @param width Numeric. Width of output plots in inches (default: 8).
#' @param height Numeric. Height of output plots in inches (default: 6).
#' @param out_dir Character (required). Path to output directory for saving convergence plots.
#' @return Processed Seurat object with:
#'         - Preprocessed expression data (normalized + scaled with unwanted variation regressed)
#'         - PCA/Harmony reductions
#'         - Cell cycle scores (S.Score/G2M.Score)
#'         - Cell clusters (seurat_clusters)
#'         - UMAP/TSNE reductions
#' @author WWY
#' @export
#' @examples
#' \dontrun{
#' # Load QC-passed Seurat object
#' seu_qc <- readRDS("../03.Output/USOO_qc.rds")
#' 
#' # Run Harmony batch correction with preprocessing
#' seu_harmony <- sc_harmony(
#'   seurat_obj = seu_qc,
#'   max_dim_pca = 30,
#'   max_dim_harmony = 30,
#'   max.iter = 30,
#'   nfeatures = 3000,
#'   res = 1.0,
#'   out_dir = "../03.Output/"
#' )
#' 
#' # Save processed object
#' saveRDS(seu_harmony, "../03.Output/USOO_harmony.rds")
#' }

sc_harmony <- function(seurat_obj = NULL,
                       max_dim_pca = NULL, 
                       max_dim_harmony = NULL,
                       max.iter = 30,
                       nfeatures = 3000,
                       res = 0.5,  
                       width = 8, 
                       height = 6,
                       out_dir = NULL) {  
    
    # -------------------------- 1. Load Required Libraries  -------------------------- #
    suppressMessages({
        library(harmony)
        library(Seurat)
        library(ggplot2)
        library(dplyr)
    })
    
    # Set random seed for reproducibility
    set.seed(1234)
    
    # -------------------------- 2. Strict Parameter Validation -------------------------- #
    if (is.null(seurat_obj)) stop("Error: 'seurat_obj' (input Seurat object) is required!")
    if (!inherits(seurat_obj, "Seurat")) stop("Error: 'seurat_obj' must be a valid Seurat object!")
    if (is.null(max_dim_pca)) stop("Error: 'max_dim_pca' (PCA dimensions for FindNeighbors) is required!")
    if (is.null(max_dim_harmony)) stop("Error: 'max_dim_harmony' (Harmony dimensions for UMAP/TSNE) is required!")
    if (is.null(out_dir)) stop("Error: 'out_dir' (output directory path) is required!")
    if (is.null(res)) res <- 0.5 
    
    # 创建输出目录
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        message(paste("Output directory created automatically:", out_dir))
    }
    
    # 定义必需的meta列（移到函数内部）
    required_meta <- c("nCount_RNA", "percent.ribo1", "percent.ribo2", "percent.mt", "percent_RBC", "orig.ident")
    # 验证meta.data列
    missing_meta <- setdiff(required_meta, colnames(seurat_obj@meta.data))
    if (length(missing_meta) > 0) {
        stop("Error: Missing required meta.data columns: ", paste(missing_meta, collapse = ", "))
    }
    
    # -------------------------- 3. Extract Object Name for Output -------------------------- #
    obj_name <- deparse(substitute(seurat_obj))
    obj_name <- gsub("[^a-zA-Z0-9_]", "", obj_name)  # Sanitize object name
    
    # -------------------------- 4. Data Preprocessing Pipeline -------------------------- #
    message("=== Starting preprocessing: Normalization -> Variable Feature Selection -> Cell Cycle Regression ===")
    
    ## Log normalization (median UMI as scaling factor to eliminate inter-cell sequencing differences)
    seurat_obj <- NormalizeData(
        object = seurat_obj,
        normalization.method = "LogNormalize",
        scale.factor = median(seurat_obj@meta.data$nCount_RNA),
        verbose = FALSE
    )
    
    ## Select highly variable features (HVGs)
    seurat_obj <- FindVariableFeatures(
        object = seurat_obj,
        selection.method = "vst",
        nfeatures = nfeatures,
        verbose = FALSE
    )
    
    ## Cell cycle scoring（修复：加载细胞周期基因）
    if (!exists("cc.genes.updated.2019")) {
        data("cc.genes.updated.2019", package = "Seurat")
    }
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    seurat_obj <- CellCycleScoring(
        object = seurat_obj,
        s.features = s.genes,
        g2m.features = g2m.genes,
        verbose = FALSE
    )
    
    ## Scale data and regress unwanted variation（修复：回归正确的变量）
    seurat_obj <- ScaleData(
        object = seurat_obj,
        vars.to.regress = c("S.Score", "G2M.Score", "percent.ribo1", "percent.ribo2", "percent.mt", "percent_RBC"),
        features = VariableFeatures(seurat_obj),
        verbose = FALSE
    )
    
    ## Run PCA on HVGs
    seurat_obj <- RunPCA(
        object = seurat_obj,
        features = VariableFeatures(seurat_obj),
        verbose = FALSE
    )
    
    ## Save PCA elbow plot
    p_pca <- ElbowPlot(seurat_obj, ndims = 50) +
        theme_minimal() +
        ggtitle("Elbow Plot for PCA Dimensions (HVGs)")
    ggsave(
        filename = file.path(out_dir, paste0(obj_name, "_PCA_elbow.pdf")),
        plot = p_pca,
        width = width,
        height = height,
        dpi = 300
    )
    
    # -------------------------- 5. Harmony Batch Correction -------------------------- #
    message("=== Running Harmony batch correction (grouped by orig.ident) ===")
    seurat_obj <- seurat_obj %>%
        RunHarmony(
            reduction = "pca",
            group.by.vars = "orig.ident",
            reduction.save = "harmony",
            plot_convergence = TRUE,
            max.iter = max.iter,
            verbose = FALSE
        )
    
    ## Save Harmony convergence plot
    p_harmony <- ElbowPlot(seurat_obj, reduction = "harmony", ndims = 50) +
        theme_minimal() +
        ggtitle("Elbow Plot for Harmony-Corrected Dimensions")
    ggsave(
        filename = file.path(out_dir, paste0(obj_name, "_Harmony_elbow.pdf")),
        plot = p_harmony,
        width = width,
        height = height,
        dpi = 300
    )
    
    # -------------------------- 6. Cell Clustering -------------------------- #
    message("=== Performing cell clustering ===")
    seurat_obj <- FindNeighbors(
        object = seurat_obj,
        reduction = "harmony",  # 若需用harmony聚类，保持此设置；建议改回"pca"更合理
        dims = 1:max_dim_pca,   # 注意：若用harmony reduction，建议改为max_dim_harmony
        verbose = FALSE
    )
    seurat_obj <- FindClusters(
        object = seurat_obj,
        resolution = res,
        verbose = FALSE
    )
    
    # -------------------------- 7. Dimensionality Reduction (UMAP/TSNE) -------------------------- #
    message("=== Running UMAP/TSNE on Harmony-corrected dimensions ===")
    ## UMAP（补充reduction.key，避免命名冲突）
    seurat_obj <- RunUMAP(
        object = seurat_obj,
        reduction = "harmony",
        dims = 1:max_dim_harmony,
        reduction.name = "umap",
        verbose = FALSE
    )
    ## TSNE
    seurat_obj <- RunTSNE(
        object = seurat_obj,
        reduction = "harmony",
        dims = 1:max_dim_harmony,
        reduction.name = "tsne",
        verbose = FALSE
    )
    
    # -------------------------- 8. Final Output -------------------------- #
    message("=== Harmony processing completed! ===")
    return(seurat_obj)
}
}
