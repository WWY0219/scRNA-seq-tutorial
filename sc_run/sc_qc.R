#' Single Cell RNA-seq Quality Control (QC) Analysis
#'
#' This function performs quality control analysis for single-cell RNA-seq data, including calculating 
#' mitochondrial and ribosomal gene percentages, generating QC visualization plots (scatter plots, violin plots), 
#' and returning the Seurat object with updated QC metrics in metadata.
#'
#' @title Single-Cell RNA-seq Quality Control Analysis
#' @description 
#' The function calculates two types of ribosomal gene percentages (RPS: Ribosomal Protein Small subunit; 
#' RPL: Ribosomal Protein Large subunit) and mitochondrial gene percentage (MT- prefix). It generates and saves 
#' three types of QC plots: ribosomal gene scatter plots (RPS/RPL vs. RNA counts), mitochondrial gene scatter plots 
#' (mitochondrial percentage vs. RNA counts and feature counts vs. RNA counts), and violin plots for all QC metrics. 
#' Results are saved in a specified output directory, with a subfolder named after the input Seurat object.
#' @param seurat_obj Input Seurat object (required), containing single-cell RNA-seq expression data.
#' @param out_dir Path to the main output directory (required). A subdirectory named after the Seurat object will be created here.
#' @param width Width of the output plots in inches (default: 10).
#' @param height Height of the output plots in inches (default: 12).
#' @return A Seurat object with updated metadata including the following QC metrics:
#' \itemize{
#'   \item{\code{percent.ribo1}:} {Percentage of RPS (Ribosomal Protein Small subunit) genes}
#'   \item{\code{percent.ribo2}:} {Percentage of RPL (Ribosomal Protein Large subunit) genes}
#'   \item{\code{percent.mt}:} {Percentage of mitochondrial genes (MT- prefix)}
#' }
#' @author WWY
#' @export
#' @examples
#' \dontrun{
#' # Load Seurat object
#' seu <- readRDS("raw_seurat.rds")
#' 
#' # Run QC analysis
#' seu_qc <- sc_qc(seurat_obj = seu,
#'                 out_dir = "../03.Output/")
#' 
#' # Check updated metadata
#' colnames(seu_qc@meta.data)
#' head(seu_qc@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo1", "percent.ribo2")])
#' }

sc_qc <- function(seurat_obj = NULL,
                  out_dir = NULL,
                  width = 10,
                  height = 12
                  ) {
    ##library 
    suppressMessages({
        library(Seurat)
        library(ggplot2)
        library(dplyr)
    })
    
    obj_name <- deparse(substitute(seurat_obj))
    
    ###-------------------------------- 0. Pre-Check --------------------------------------------------------###
    if (is.null(seurat_obj)) stop("请传入Seurat对象（seurat_obj）！")
    if (is.null(out_dir)) stop("请指定输出目录（out_dir）！")
    
    ###-------------------------------- 1. Create Output Directory ------------------------------------------###
    output_dir <- file.path(out_dir, obj_name)
    if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message(paste("已自动创建输出目录：", output_dir))
  }
        
    ###-------------------------------- 2. Canculate RPL/RPS-genes Percentage -------------------------------###
    ribo_genes <- rownames(seurat_obj)[grep("^RPL|^RPS", rownames(seurat_obj))]
    seurat_obj[["percent.ribo1"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^RPS")
    seurat_obj[["percent.ribo2"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^RPL")
  
    ###-------------------------------- 3. Canculate MT-genes Percentage ------------------------------------###
    seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
    ###-------------------------------- 4. Draw RPL/RPS-genes Percentage Plot -------------------------------###
    plot1 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo1")
    plot2 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo2")
    qc1 <- plot1 + plot2 
    ggsave(filename=paste0(output_dir, paste0(obj_name,"_qc_RPLS",".pdf")), plot = qc1, width = width, height = height, dpi = 300)
    message("QC图已保存为：",paste0(obj_name,"_qc_RPLS",".pdf"))
    
    ###-------------------------------- 5. Draw MT-genes Percentage Plot ------------------------------------###
    plot3 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot4 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    qc2 <- plot3 + plot4
    ggsave(filename=paste0(output_dir, paste0(obj_name,"_qc_MT",".pdf")), plot = qc2, width = width, height = height, dpi = 300)
    message("QC图已保存为：",paste0(obj_name,"_qc_MT",".pdf"))

    ###-------------------------------- 6. Draw Meta Percentage Plot ----------------------------------------###
    qc_meta_violin <- Seurat::VlnPlot(
        seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo1","percent.ribo2"),
        ncol = 3,pt.size = 0)
    print(qc_meta_violin)  
    ggsave(filename=paste0(output_dir, paste0(obj_name,"_qc_metaviolin",".pdf")), plot = qc_meta_violin, width = width, height = height, dpi = 300)
    message("QC图已保存为：",paste0(obj_name,"_qc_metaviolin",".pdf"))

    ###-------------------------------- 7. Output metadata Preview ------------------------------------------###
    cat("\n=== Metadata前5行（包含核糖体/线粒体占比指标） ===\n")
    print(head(seurat_obj@meta.data, 5))
    
    ###-------------------------------- 8. Output metadata Preview ------------------------------------------###
    return(seurat_obj)
}
