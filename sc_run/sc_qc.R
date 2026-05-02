#' Single Cell RNA-seq Quality Control (QC) Analysis
#'
#' This function performs quality control analysis for single-cell RNA-seq data, including calculating 
#' mitochondrial and ribosomal gene percentages, generating QC visualization plots (scatter plots, violin plots), 
#' and returning the Seurat object with updated QC metrics in metadata.
#'
#' @title Single-Cell RNA-seq Quality Control Analysis
#' @description 
#' The function calculates two types of ribosomal gene percentages (RPS: Ribosomal Protein Small subunit; 
#' RPL: Ribosomal Protein Large subunit) and mitochondrial gene percentage (MT-/mt- prefix) based on the specified species. 
#' It generates and saves three types of QC plots: ribosomal gene scatter plots (RPS/RPL vs. RNA counts), 
#' mitochondrial gene scatter plots (mitochondrial percentage vs. RNA counts and feature counts vs. RNA counts), 
#' and violin plots for all QC metrics. 
#' Results are saved in a specified output directory, with a subfolder named after the input Seurat object.
#' @param seurat_obj Input Seurat object (required), containing single-cell RNA-seq expression data.
#' @param species Species of the sample, either "human" or "mouse" (default: "human"). Determines the gene symbol patterns.
#' @param out_dir Path to the main output directory (required). A subdirectory named after the Seurat object will be created here.
#' @param width Width of the output plots in inches (default: 10).
#' @param height Height of the output plots in inches (default: 12).
#' @param dpi Resolution of output plots (default: 300).
#' @return A Seurat object with updated metadata including the following QC metrics:
#' \itemize{
#'   \item{\code{percent.ribo1}:} {Percentage of RPS (Ribosomal Protein Small subunit) genes}
#'   \item{\code{percent.ribo2}:} {Percentage of RPL (Ribosomal Protein Large subunit) genes}
#'   \item{\code{percent.mt}:} {Percentage of mitochondrial genes (MT-/mt- prefix)}
#'   \item{\code{percent_RBC}:} {Percentage of hemoglobin genes (HBAB/Hba/Hbb prefix)}
#' }
#' @author WWY
#' @export
#' @examples
#' \dontrun{
#' # Load Seurat object
#' seu <- readRDS("raw_seurat.rds")
#' 
#' # Run QC analysis for human
#' seu_qc_human <- sc_qc(seurat_obj = seu, species = "human", out_dir = "./03.Output/")
#' 
#' # Run QC analysis for mouse
#' seu_qc_mouse <- sc_qc(seurat_obj = seu, species = "mouse", out_dir = "./03.Output/")
#' }

sc_qc <- function(seurat_obj = NULL,
                  species = "human",
                  out_dir = NULL,
                  width = 10,
                  height = 12,
                  dpi = 300) {
    # Load required R packages
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
    
    # Get name of Seurat object
    obj_name <- deparse(substitute(seurat_obj))
    
    # ---------------------------------------- 00. Pre-Check ---------------------------------------------------------
    if (is.null(seurat_obj)) stop("Please provide a Seurat object (seurat_obj)!")
    if (is.null(out_dir)) stop("Please specify the output directory (out_dir)!")
    if (!tolower(species) %in% c("human", "mouse")) stop("The 'species' parameter must be either 'human' or 'mouse'!")
    
    # ---------------------------------------- 01. Create Output Directory -------------------------------------------
    output_dir <- file.path(out_dir, obj_name)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        message(paste("Output directory automatically created:", output_dir))
    }
        
    # ---------------------------------------- 02. Quality Control: Define Patterns & Calculate Metrics -------------
    # Assign regular expressions based on species
    if (tolower(species) == "human") {
        pattern_mt   <- "^MT-"
        pattern_rps  <- "^RPS"
        pattern_rpl  <- "^RPL"
        pattern_rbc  <- "^HB[AB]"
    } else if (tolower(species) == "mouse") {
        pattern_mt   <- "^mt-"
        pattern_rps  <- "^Rps"
        pattern_rpl  <- "^Rpl"
        pattern_rbc  <- "^Hb[ab]"
    }
    
    # Calculate mitochondrial gene percentage
    seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = pattern_mt)

    # Calculate ribosomal gene percentages (RPS/RPL)
    seurat_obj[["percent.ribo1"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = pattern_rps)
    seurat_obj[["percent.ribo2"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = pattern_rpl)

    # Calculate hemoglobin (RBC) gene percentage
    seurat_obj[["percent_RBC"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = pattern_rbc)
  
    # ---------------------------------------- 03. Generate QC Plots ------------------------------------------------
    # RNA count vs. feature count scatter plot
    plot1 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_nRNA.pdf")), 
           plot = plot1, width = width, height = height, dpi = dpi)
    message("QC plot saved as: ", paste0(obj_name, "_qc_nRNA.pdf"))
    
    # Mitochondrial percentage scatter plot
    plot2 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_MT.pdf")), 
           plot = plot2, width = width, height = height, dpi = dpi)
    message("QC plot saved as: ", paste0(obj_name, "_qc_MT.pdf"))

    # Ribosomal gene (RPS/RPL) scatter plots
    plot3 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo1")
    plot4 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.ribo2")
    qc1 <- plot3 + plot4 
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_RPLS.pdf")), 
           plot = qc1, width = width, height = height, dpi = dpi)
    message("QC plot saved as: ", paste0(obj_name, "_qc_RPLS.pdf"))

    # Hemoglobin (RBC) percentage scatter plot
    plot5 <- Seurat::FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent_RBC")
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_RBC.pdf")), 
           plot = plot5, width = width, height = height, dpi = dpi)
    message("QC plot saved as: ", paste0(obj_name, "_qc_RBC.pdf"))
    
    # Violin plots for all QC metrics
    qc_meta_violin <- Seurat::VlnPlot(
        seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo1", "percent.ribo2", "percent_RBC"),
        ncol = 3, pt.size = 0)
    print(qc_meta_violin)  
    ggsave(filename = file.path(output_dir, paste0(obj_name, "_qc_metaviolin.pdf")), 
           plot = qc_meta_violin, width = width, height = height, dpi = dpi)
    message("QC plot saved as: ", paste0(obj_name, "_qc_metaviolin.pdf"))

    # ---------------------------------------- 04. Metadata Preview -----------------------------------------------
    cat(sprintf("\n=== First 5 rows of metadata (Species: %s) ===\n", species))
    print(head(seurat_obj@meta.data, 5))
    
    # ---------------------------------------- 05. Return Updated Seurat Object -----------------------------------
    return(seurat_obj)
}
