#' Batch Processing of Single-Cell Data with Harmony (Retain All Resolutions)
#'
#' This function performs Harmony batch effect correction on single-cell RNA sequencing data, 
#' systematically tests multiple clustering resolutions, and retains all resolution-specific 
#' clustering results in the metadata of the input Seurat object. It also generates comprehensive
#' visualizations (UMAP/tSNE plots) and marker gene analyses for each resolution, facilitating
#' optimal cluster selection for downstream biological interpretation.
#'
#' @title Single-Cell Harmony Batch Processing with Multi-Resolution Clustering
#' @param seurat_obj Input Seurat object containing preprocessed scRNA-seq data (required).
#' @param max_dim_pca Number of principal components (PCs) used for clustering/UMAP before Harmony correction (required, e.g., 30).
#' @param max_dim_harmony Number of Harmony-corrected dimensions used for clustering/UMAP after batch correction (required, e.g., 30).
#' @param resolutions Numeric vector of clustering resolutions to test (required, e.g., c(0.5, 0.8, 1.0, 1.5)).
#' @param max_iter Maximum number of iterations for Harmony algorithm (default: 10).
#' @param top_n Number of top marker genes to retain per cluster (required, e.g., 10).
#' @param width Width of output plots in inches (default: 10).
#' @param height Height of output plots in inches (default: 8).
#' @param dpi Resolution of output plots in dots per inch (default: 300).
#' @param out_dir Path to the main output directory (required). Will be created if it does not exist.
#' @return Seurat object with all resolution-specific clustering results stored in metadata (columns named: cluster_res{resolution}).
#' @author WWY
#' @export
#' @examples
#' \dontrun{
#' # Load preprocessed Seurat object (ensure "orig.ident" is present in metadata)
#' seu <- readRDS("preprocessed_seurat_object.rds")
#' 
#' # Run multi-resolution Harmony processing
#' seu_all_res <- sc_resolutionfinder(
#'   seurat_obj = seu,
#'   max_dim_pca = 30,
#'   max_dim_harmony = 30,
#'   resolutions = c(0.5, 0.8, 1.0, 1.5),
#'   top_n = 10,
#'   out_dir = "../05.Harmony_MultiResolution"
#' )
#' 
#' # View all resolution-specific clustering columns in metadata
#' colnames(seu_all_res@meta.data)[grep("cluster_res", colnames(seu_all_res@meta.data))]
#' }

sc_resolutionfinder <- function(seurat_obj = NULL,
                                max_dim_pca = NULL,
                                max_dim_harmony = NULL,
                                resolutions = NULL,   
                                max_iter = 30,
                                top_n = NULL,
                                width = 10,
                                height = 8,
                                dpi = 300,
                                out_dir = NULL) {
  # Set environment and global options
  Sys.setenv(LANGUAGE = "en")
  options(stringsAsFactors = FALSE, warn = -1)
  
  # Check for required packages and install missing ones
  required_packages <- c("qs", "harmony", "Seurat", "ggplot2", "dplyr", "ggthemes", "openxlsx", "tidyr")
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    stop(paste("The following required packages are not installed:", paste(missing_packages, collapse = ", ")))
  }
  
  # Load required libraries (suppress warning messages)
  suppressMessages({
    library(qs)
    library(harmony)
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(ggthemes)
    library(openxlsx)
    library(tidyr)
  })
  
  # Set random seed for reproducibility
  set.seed(1234) 
  
  # Extract input object name for file naming
  obj_name <- deparse(substitute(seurat_obj))
  
  #--------------------------------------------------------------------- 1. Parameter Validation -------------------------------------------------------------------
  if (is.null(seurat_obj)) stop("Error: Please provide a valid Seurat object to 'seurat_obj'.")
  if (is.null(max_dim_pca)) stop("Error: Please specify the number of PCA dimensions via 'max_dim_pca'.")
  if (is.null(max_dim_harmony)) stop("Error: Please specify the number of Harmony dimensions via 'max_dim_harmony'.")
  if (is.null(resolutions)) stop("Error: Please provide a vector of clustering resolutions via 'resolutions'.")
  if (is.null(top_n)) stop("Error: Please specify the number of top marker genes via 'top_n'.")
  if (is.null(out_dir)) stop("Error: Please specify the output directory via 'out_dir'.")
  
  # Create main output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message(paste("Main output directory created:", out_dir))
  }
  
  #--------------------------------------------------------------------- 2. Global Preprocessing & Harmony Correction ---------------------------------------------------
  message("===== Starting global preprocessing and Harmony batch correction =====")
  
  # Merge RNA assay layers (compatible with Seurat v5+)
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  
  
  # Run Harmony correction if not already performed
  if (!"harmony" %in% names(seurat_obj@reductions)) {
    if (!"pca" %in% names(seurat_obj@reductions)) {
      stop("Error: PCA reduction not found in the Seurat object. Please run RunPCA() before using this function.")
    }
    seurat_obj <- RunHarmony(
      object = seurat_obj,
      group.by.vars = "orig.ident",  # Batch variable (assumed to be stored in "orig.ident")
      dims.use = 1:max_dim_pca,
      max.iter = max_iter,
      verbose = FALSE
    )
    message("Harmony batch correction completed successfully.")
  } else {
    message("Using existing Harmony reduction from the input Seurat object.")
  }
  
  #--------------------------------------------------------------------- 3. Multi-Resolution Clustering Pipeline -------------------------------------------------------
  for (res in resolutions) {
    message(paste("\n===== Processing resolution:", res, "====="))
    
    # Create resolution-specific subdirectory
    res_out_dir <- file.path(out_dir, paste0(obj_name, "_resolution_", res))
    if (!dir.exists(res_out_dir)) {
      dir.create(res_out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Define column name for current resolution's clustering results
    cluster_col <- paste0("cluster_res", res)
    
    # ---------------------- Pre-Harmony: Clustering & Visualization (for comparison) ----------------------
    seurat_temp_pre <- seurat_obj
    # Run clustering on raw PCA dimensions
    seurat_temp_pre <- FindNeighbors(seurat_temp_pre, reduction = "pca", dims = 1:max_dim_pca, verbose = FALSE)
    seurat_temp_pre <- FindClusters(seurat_temp_pre, resolution = res, verbose = FALSE)
    # Run UMAP/tSNE for visualization
    seurat_temp_pre <- RunUMAP(seurat_temp_pre, reduction = "pca", dims = 1:max_dim_pca, 
                               reduction.name = "umap_naive_pre", verbose = FALSE)
    seurat_temp_pre <- RunTSNE(seurat_temp_pre, reduction = "pca", dims = 1:max_dim_pca, 
                               reduction.name = "tsne_naive_pre", verbose = FALSE)
    
    # Save pre-Harmony visualization plots
    # UMAP/tSNE with cluster labels
    p_umap_pre <- DimPlot(seurat_temp_pre, reduction = "umap_naive_pre", label = TRUE, pt.size = 0.3) + 
      coord_fixed() + 
      ggtitle(paste0("UMAP before Harmony (Resolution = ", res, ")")) +
      theme_minimal()
    ggsave(p_umap_pre, file = file.path(res_out_dir, paste0(obj_name, "_umap_before_Harmony.pdf")), 
           width = width, height = height, dpi = dpi)
    p_tsne_pre <- DimPlot(seurat_temp_pre, reduction = "tsne_naive_pre", label = TRUE, pt.size = 0.3) + 
      coord_fixed() + 
      ggtitle(paste0("tSNE before Harmony (Resolution = ", res, ")")) +
      theme_minimal()
    ggsave(p_tsne_pre, file = file.path(res_out_dir, paste0(obj_name, "_tsne_before_Harmony.pdf")), 
           width = width, height = height, dpi = dpi)
    
    # UMAP/tSNE colored by sample (orig.ident)
    p_umap_pre_sample <- DimPlot(seurat_temp_pre, reduction = "umap_naive_pre", group.by = "orig.ident", 
                                 label = FALSE, pt.size = 0.2) + 
      coord_fixed() + 
      ggtitle(paste0("UMAP by Sample (Before Harmony, Res=", res, ")")) +
      theme_minimal()
    ggsave(p_umap_pre_sample, file = file.path(res_out_dir, paste0(obj_name, "_umap_sample_before_Harmony.pdf")), 
           width = width, height = height, dpi = dpi)
      p_tsne_pre_sample <- DimPlot(seurat_temp_pre, reduction = "tsne_naive_pre", group.by = "orig.ident", 
                                 label = FALSE, pt.size = 0.2) + 
      coord_fixed() + 
      ggtitle(paste0("tSNE by Sample (Before Harmony, Res=", res, ")")) +
      theme_minimal()
    ggsave(p_tsne_pre_sample, file = file.path(res_out_dir, paste0(obj_name, "_tsne_sample_before_Harmony.pdf")), 
           width = width, height = height, dpi = dpi)
    
    # UMAP/tSNE split by sample
    p_umap_pre_sample_split <- DimPlot(seurat_temp_pre, reduction = "umap_naive_pre", split.by = "orig.ident", 
                                       label = FALSE, pt.size = 0.2) + 
      coord_fixed() + 
      ggtitle(paste0("UMAP Split by Sample (Before Harmony, Res=", res, ")")) +
      theme_minimal()
    ggsave(p_umap_pre_sample_split, file = file.path(res_out_dir, paste0(obj_name, "_umap_sample_split_before_Harmony.pdf")), 
           width = 15, height = 8, dpi = dpi)
      p_tsne_pre_sample_split <- DimPlot(seurat_temp_pre, reduction = "tsne_naive_pre", split.by = "orig.ident", 
                                       label = FALSE, pt.size = 0.2) + 
      coord_fixed() + 
      ggtitle(paste0("tSNE Split by Sample (Before Harmony, Res=", res, ")")) +
      theme_minimal()
    ggsave(p_tsne_pre_sample_split, file = file.path(res_out_dir, paste0(obj_name, "_tsne_sample_split_before_Harmony.pdf")), 
           width = 15, height = 8, dpi = dpi)
    
    # ---------------------- Post-Harmony: Clustering (Save to Original Object) ----------------------
    # Run clustering on Harmony-corrected dimensions
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:max_dim_harmony, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = cluster_col, verbose = FALSE)
    
    # Count number of clusters for current resolution
    n_cluster <- length(unique(seurat_obj[[cluster_col]][, 1]))
    message(paste("Post-Harmony clustering completed:", n_cluster, "clusters (Metadata column:", cluster_col, ")"))
    
    # ---------------------- Post-Harmony: Visualization & Marker Gene Analysis ----------------------
    seurat_temp_post <- seurat_obj
    # Run resolution-specific UMAP/tSNE (avoid overwriting global reductions)
    seurat_temp_post <- RunUMAP(seurat_temp_post, reduction = "harmony", dims = 1:max_dim_harmony, 
                                reduction.name = paste0("umap_harmony_res", res), verbose = FALSE)
    seurat_temp_post <- RunTSNE(seurat_temp_post, reduction = "harmony", dims = 1:max_dim_harmony, 
                                reduction.name = paste0("tsne_harmony_res", res), verbose = FALSE)
    
    # Save Harmony reduction quality plots
    # Elbow plot for Harmony dimensions
    p_elbow <- ElbowPlot(seurat_obj, reduction = "harmony", ndims = 60) + theme_minimal()
    ggsave(p_elbow, file = file.path(res_out_dir, paste0(obj_name, "_harmony_elbow_plot.pdf")), 
           width = width, height = height, dpi = dpi)
    
    # Post-Harmony UMAP/tSNE with cluster labels
    p_umap_post <- DimPlot(seurat_temp_post, reduction = paste0("umap_harmony_res", res), 
                           group.by = cluster_col, label = TRUE, pt.size = 0.3) + 
      coord_fixed() + 
      ggtitle(paste0("UMAP after Harmony (Resolution = ", res, ")")) +
      theme_minimal()
    ggsave(p_umap_post, file = file.path(res_out_dir, paste0(obj_name, "_umap_after_Harmony.pdf")), 
           width = width, height = height, dpi = dpi)
     p_tsne_post <- DimPlot(seurat_temp_post, reduction = paste0("tsne_harmony_res", res), 
                           group.by = cluster_col, label = TRUE, pt.size = 0.3) + 
      coord_fixed() + 
      ggtitle(paste0("tSNE after Harmony (Resolution = ", res, ")")) +
      theme_minimal()
    ggsave(p_tsne_post, file = file.path(res_out_dir, paste0(obj_name, "_tsne_after_Harmony.pdf")), 
           width = width, height = height, dpi = dpi)
    
    # Post-Harmony UMAP/tSNE split by sample
      p_umap_post_sample_split <- DimPlot(seurat_temp_post, reduction = paste0("umap_harmony_res", res), 
                                        split.by = "orig.ident", label = FALSE, pt.size = 0.3) + 
      coord_fixed() + 
      ggtitle(paste0("UMAP Split by Sample (After Harmony, Res=", res, ")")) +
      theme_minimal()
    ggsave(p_umap_post_sample_split, file = file.path(res_out_dir, paste0(obj_name, "_umap_sample_split_after_Harmony.pdf")), 
           width = 15, height = 8, dpi = dpi)
    p_tsne_post_sample_split <- DimPlot(seurat_temp_post, reduction = paste0("tsne_harmony_res", res), 
                                        split.by = "orig.ident", label = FALSE, pt.size = 0.3) + 
      coord_fixed() + 
      ggtitle(paste0("tSNE Split by Sample (After Harmony, Res=", res, ")")) +
      theme_minimal()
    ggsave(p_tsne_post_sample_split, file = file.path(res_out_dir, paste0(obj_name, "_tsne_sample_split_after_Harmony.pdf")), 
           width = 15, height = 8, dpi = dpi)
    
    # ---------------------- Marker Gene Identification & Export ----------------------
    message(paste("Calculating marker genes for resolution", res, "..."))
    # Identify positive marker genes for each cluster
    markers <- FindAllMarkers(seurat_obj,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              group.by = cluster_col,
                              verbose = FALSE)
    
    # Filter significant markers (adjusted p-value < 0.01 & log2FC > 0.25)
    markers_significant <- markers %>%
      filter(avg_log2FC > 0.25 & p_val_adj < 0.01)
    
    # Export significant markers
    write.xlsx(markers_significant, 
               file = file.path(res_out_dir, paste0(obj_name, "_significant_marker_genes.xlsx")),
               overwrite = TRUE)
    message(paste("Significant marker genes saved:", nrow(markers_significant), "genes"))
    
    # Extract top N marker genes per cluster
    markers_top <- markers_significant %>%
      group_by(cluster) %>%
      slice_max(n = top_n, order_by = avg_log2FC) %>%
      ungroup() %>%
      select(gene, cluster) %>%
      mutate(cluster = paste0("Cluster_", cluster))
    
    # Ensure each cluster has exactly top_n entries (fill with NA if insufficient)
    markers_top_complete <- expand.grid(
      cluster = unique(markers_top$cluster),
      rank = 1:top_n
    ) %>%
      left_join(
        markers_top %>%
          group_by(cluster) %>%
          mutate(rank = row_number()) %>%
          ungroup(),
        by = c("cluster", "rank")
      ) %>%
      select(-rank)
    
    # Reshape to wide format and export
    markers_top_wide <- markers_top_complete %>%
      group_by(cluster) %>%
      mutate(column_index = row_number()) %>%
      pivot_wider(names_from = cluster, values_from = gene) %>%
      select(-column_index)
    
    write.csv(markers_top_wide, 
              file = file.path(res_out_dir, paste0(obj_name, "_top", top_n, "_marker_genes.csv")),
              row.names = FALSE, na = "")
    write.xlsx(markers_top_wide, 
               file = file.path(res_out_dir, paste0(obj_name, "_top", top_n, "_marker_genes.xlsx")),
               overwrite = TRUE, na.string = "")
    
    # Save temporary Seurat object for current resolution
    qsave(seurat_temp_post, 
          file = file.path(res_out_dir, paste0(obj_name, "_seurat_resolution_", res, ".qs")))
    
    # Clean up temporary objects to free memory
    rm(seurat_temp_pre, seurat_temp_post)
    gc()
  }
  
  #--------------------------------------------------------------------- 4. Final Output & Return -------------------------------------------------------------------
  # Save the final Seurat object with all resolutions
  final_obj_path <- file.path(out_dir, paste0(obj_name, "_multi_resolution_processed.qs"))
  qsave(seurat_obj, file = final_obj_path)
  
  # Print completion message
  message("\n===== All resolutions processed successfully! =====")
  message("Clustering results are stored in the Seurat object's metadata with columns:")
  message(paste("- Column naming convention:", "cluster_res{resolution} (e.g., cluster_res0.5)"))
  message(paste("- Example to view columns:", "colnames(seurat_obj@meta.data)[grep('cluster_res', colnames(seurat_obj@meta.data))]"))
  message(paste("- Final Seurat object saved to:", final_obj_path))
  
  return(seurat_obj)
}
