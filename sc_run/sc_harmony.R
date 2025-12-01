#' Batch Processing of Single-Cell Data with Harmony (Retain All Resolutions)
#'
#' This function performs Harmony batch correction, batch-tests multiple clustering resolutions, 
#' and retains ALL resolution-specific clustering results in the input Seurat object's metadata.
#'
#' @title Single-Cell Harmony Processing (Retain All Resolutions)
#' @param seurat_object Input Seurat object (required).
#' @param max.dim Number of PCA dimensions used for clustering/UMAP BEFORE Harmony (required, e.g., 30).
#' @param max.dim_harmony Number of Harmony dimensions used for clustering/UMAP AFTER Harmony (required, e.g., 30).
#' @param resolutions Vector of resolutions for batch clustering (required, e.g., c(0.5, 0.8, 1, 1.5)).
#' @param max.iter Maximum iterations for Harmony (default: 10).
#' @param top_n Number of top marker genes to retain per cluster (required).
#' @param out_dir Main output directory (required).
#' @return Seurat object with ALL resolution clustering results (columns: clusters_res{res}).
#' @author WWY
#' @export
#' @examples
#' \dontrun{
#' seu <- readRDS("preprocessed_seu.rds") # Seurat object with "orig.ident" sample metadata
#' seu_all_res <- sc_harmony(
#'   seurat_object = seu,
#'   max.dim = 30,
#'   max.dim_harmony = 30,
#'   resolutions = c(0.5, 0.8, 1, 1.5),
#'   top_n = 10,
#'   out_dir = "../05.Harmony_AllRes"
#' )
#' # View all resolution-specific columns in metadata
#' colnames(seu_all_res@meta.data)[grep("clusters_res", colnames(seu_all_res@meta.data))]
#' }

sc_harmony <- function(seurat_object = NULL,
                       max.dim = NULL,
                       max.dim_harmony = NULL, 
                       resolutions = NULL,   
                       max.iter = 10, 
                       top_n = NULL, 
                       out_dir = NULL) {
    ## Load required libraries (suppress warning messages)
    suppressMessages({
        library(DoubletFinder)
        library(Seurat)
        library(ggplot2)
        library(dplyr)
        library(ggthemes)
        library(openxlsx)
    })
    set.seed(1234) # Set random seed for reproducibility
  
    # Get the name of the input Seurat object for file naming
    obj_name <- deparse(substitute(seurat_object))
  
    ###--------------------------------------------------------------------- 1. Parameter Validation -------------------------------------------------------------------###
    if (is.null(seurat_object)) stop("Please provide a Seurat object (seurat_object)!")
    if (is.null(max.dim)) stop("Please specify the number of PCA dimensions for clustering/UMAP BEFORE Harmony (max.dim)!")
    if (is.null(max.dim_harmony)) stop("Please specify the number of Harmony dimensions for clustering/UMAP AFTER Harmony (max.dim_harmony)!")
    if (is.null(resolutions)) stop("Please specify a vector of resolutions for batch clustering (resolutions)!")
    if (is.null(top_n)) stop("Please specify the number of top marker genes to retain per cluster (top_n)!")
    if (is.null(out_dir)) stop("Please specify the main output directory (out_dir)!")
  
    # Create main output directory if it doesn't exist
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        message("Main output directory created:", out_dir)
    }

    ###--------------------------------------------------------------------- 2. Global Preprocessing --------------------------------------------------------------###
    message("===== Starting global preprocessing (Normalization → FindVariableFeatures → PCA → Harmony) =====")
    # Normalization, variable feature selection, and PCA (before Harmony)
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
    seurat_object <- ScaleData(seurat_object, verbose = FALSE)
    seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose = FALSE)
  
    # Perform Harmony batch correction (run once)
    seurat_object <- seurat_object %>% RunHarmony(
        reduction = "pca",
        group.by.vars = "orig.ident",
        reduction.save = "harmony",    
        plot_convergence = TRUE,
        max.iter = max.iter,
        verbose = FALSE
    )
    # Save Harmony convergence plot
    ggsave(file.path(out_dir, paste0(obj_name, "_Harmony_convergence.pdf")), width = 8, height = 6, dpi = 300)
    message("===== Global preprocessing completed =====")

    ###--------------------------------------------------------------------- 3. Batch Process Each Resolution (Retain All Results) -------------------------------------------------------------------###
    for (res in resolutions) {
        message("\n===== Starting processing for resolution:", res, "=====")
        # Create resolution-specific subdirectory
        res_out_dir <- file.path(out_dir, paste0(obj_name, "_resolution_", res))
        if (!dir.exists(res_out_dir)) {
            dir.create(res_out_dir, recursive = TRUE)
        }
    
        # Define column name for current resolution's clustering results (e.g., clusters_res0.5)
        cluster_col <- paste0("clusters_res", res)
    
        # ---------------------- Before Harmony: Clustering + UMAP (for visualization comparison) ----------------------
        # Clustering before Harmony
        seurat_temp_pre <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:max.dim, verbose = FALSE)
        seurat_temp_pre <- FindClusters(seurat_temp_pre, resolution = res, verbose = FALSE)
        # UMAP before Harmony
        seurat_temp_pre <- RunUMAP(seurat_temp_pre, reduction = "pca", dims = 1:max.dim, reduction.name = "umap_naive_pre", verbose = FALSE)
    
        # Save visualization results before Harmony
        # UMAP plot with clusters
        p_umap_pre <- DimPlot(seurat_temp_pre, reduction = "umap_naive_pre", label = TRUE, pt.size = 0.2) + 
            coord_fixed() + ggtitle(paste0("UMAP before Harmony (resolution = ", res, ")"))
        ggsave(p_umap_pre, file = file.path(res_out_dir, paste0(obj_name, "_umap_before_Harmony.pdf")), width = 10, height = 8, dpi = 300)
    
        # UMAP plot colored by sample (orig.ident)
        p_umap_pre_orig1 <- DimPlot(seurat_temp_pre, reduction = "umap_naive_pre", group.by = "orig.ident", label = FALSE, pt.size = 0.2) + 
            coord_fixed() + ggtitle(paste0("UMAP by orig.ident (before Harmony, res=", res, ")"))
        ggsave(p_umap_pre_orig1, file = file.path(res_out_dir, paste0(obj_name, "_umap_orig.ident_before_Harmony.pdf")), width = 10, height = 8, dpi = 300)
    
        # UMAP plot split by sample (orig.ident)
        p_umap_pre_orig2 <- DimPlot(seurat_temp_pre, reduction = "umap_naive_pre", split.by = "orig.ident", label = FALSE, pt.size = 0.2) + 
            coord_fixed() + ggtitle(paste0("UMAP by SampleID (before Harmony, res=", res, ")"))
        ggsave(p_umap_pre_orig2, file = file.path(res_out_dir, paste0(obj_name, "_umap_SampleID_before_Harmony.pdf")), width = 15, height = 8, dpi = 300)
    
        # ---------------------- After Harmony: Clustering (Save Results to Original Object) ----------------------
        # Clustering based on Harmony dimensions; save results to original object's metadata
        seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:max.dim_harmony, verbose = FALSE)
        seurat_object <- FindClusters(seurat_object, 
                                      resolution = res, 
                                      cluster.name = cluster_col,  # Cluster column name includes resolution
                                      verbose = FALSE)
        # Count number of clusters for current resolution
        n_cluster <- length(unique(seurat_object[[cluster_col]]))
        message("Post-Harmony clustering completed for resolution", res, ":", n_cluster, "clusters (result column:", cluster_col, ")")
    
        # ---------------------- After Harmony: UMAP (Temporary Visualization, No Overwrite) ----------------------
        # Create temporary object for UMAP (avoid overwriting original object's UMAP)
        seurat_temp_post <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:max.dim_harmony, 
                                    reduction.name = paste0("umap_harmony_res", res), verbose = FALSE)
    
        # Save visualization results after Harmony
        # Elbow plot and dimension heatmap
        p_elo <- ElbowPlot(seurat_object, reduction = "harmony", ndims = 60)
        ggsave(p_elo, file = file.path(res_out_dir, paste0(obj_name, "_ElbowPlot_after_Harmony.pdf")), width = 10, height = 8, dpi = 300)
    
        p_dimheatmap <- DimHeatmap(seurat_object, reduction = "harmony", dims = 1:20, cells = 500, balanced = TRUE)
        ggsave(p_dimheatmap, file = file.path(res_out_dir, paste0(obj_name, "_DimHeatmap_after_Harmony.pdf")), width = 10, height = 8, dpi = 300)
    
        # UMAP plot with clusters (using current resolution's cluster column)
        p_umap_post <- DimPlot(seurat_temp_post, reduction = paste0("umap_harmony_res", res), 
                               group.by = cluster_col, label = TRUE, pt.size = 0.2) + 
            coord_fixed() + ggtitle(paste0("UMAP after Harmony (resolution = ", res, ")"))
        ggsave(p_umap_post, file = file.path(res_out_dir, paste0(obj_name, "_umap_after_Harmony.pdf")), width = 10, height = 8, dpi = 300)
    
        # UMAP plot colored by sample (orig.ident)
        p_umap_post_orig1 <- DimPlot(seurat_temp_post, reduction = paste0("umap_harmony_res", res), 
                                     group.by = "orig.ident", label = FALSE, pt.size = 0.2) + 
            coord_fixed() + ggtitle(paste0("UMAP by orig.ident (after Harmony, res=", res, ")"))
        ggsave(p_umap_post_orig1, file = file.path(res_out_dir, paste0(obj_name, "_umap_orig.ident_after_Harmony.pdf")), width = 12, height = 8, dpi = 300)
    
        # UMAP plot split by sample (orig.ident)
        p_umap_post_orig2 <- DimPlot(seurat_temp_post, reduction = paste0("umap_harmony_res", res), 
                                     split.by = "orig.ident", label = FALSE, pt.size = 0.2) + 
            coord_fixed() + ggtitle(paste0("UMAP by SampleID (after Harmony, res=", res, ")"))
        ggsave(p_umap_post_orig2, file = file.path(res_out_dir, paste0(obj_name, "_umap_SampleID_after_Harmony.pdf")), width = 15, height = 8, dpi = 300)
    
        # ---------------------- Marker Gene Analysis ----------------------
        message("Calculating marker genes for resolution", res, "...")
        seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
        # Find marker genes based on current resolution's cluster column
        markers <- FindAllMarkers(seurat_object, 
                                  only.pos = TRUE,  # Retain only positive markers
                                  min.pct = 0.25,   # Minimum percentage of cells expressing the gene in the cluster
                                  logfc.threshold = 0.25,  # Minimum log2 fold change
                                  group.by = cluster_col,  # Use current resolution's cluster column
                                  verbose = FALSE)
    
        # Filter significant markers and save
        markers_sig <- markers[markers$avg_log2FC > 0.25 & markers$p_val_adj < 0.01, ]
        write.xlsx(markers_sig, 
                   file = file.path(res_out_dir, paste0(obj_name, "_MarkerGene_sig.xlsx")),
                   overwrite = TRUE)
        message("Significant marker genes saved (", nrow(markers_sig), "genes)")
    
        # Extract top_n marker genes per cluster and save
        markers_top <- markers_sig %>%
            group_by(cluster) %>%
            top_n(n = top_n, wt = avg_log2FC)  # Rank by average log2 fold change  
    
        # Handle cases where clusters have fewer/more than top_n genes
        top_genes <- markers_top %>% select(gene, cluster)
        top_genes$cluster <- paste0("cluster_", top_genes$cluster)
        cluster_counts <- table(top_genes$cluster)
    
        # Fill NA for clusters with fewer than top_n genes
        clusters_under <- names(cluster_counts[cluster_counts < top_n])
        if (length(clusters_under) > 0) {
            for (clu in clusters_under) {
                need <- top_n - cluster_counts[clu]
                supple <- data.frame(gene = rep(NA, need), cluster = rep(clu, need), stringsAsFactors = FALSE)
                top_genes <- rbind(top_genes, supple)
            }
        }
    
        # Truncate clusters with more than top_n genes
        clusters_over <- names(cluster_counts[cluster_counts > top_n])
        if (length(clusters_over) > 0) {
            for (clu in clusters_over) {
                top_genes <- top_genes %>% filter(!(cluster == clu & row_number() > top_n))
            }
        }
    
        # Reshape and save top marker genes
        top_genes_unstacked <- unstack(top_genes, gene ~ cluster)
        write.csv(top_genes_unstacked, 
                  file = file.path(res_out_dir, paste0(obj_name, "_MarkerGene_top", top_n, ".csv")),
                  row.names = FALSE)
        write.xlsx(top_genes_unstacked, 
                   file = file.path(res_out_dir, paste0(obj_name, "_MarkerGene_top", top_n, ".xlsx")),
                   overwrite = TRUE)
        message("Top", top_n, "marker genes per cluster saved")
    
        # Save temporary Seurat object for current resolution (optional, for separate inspection)
        saveRDS(seurat_temp_post, 
                file = file.path(res_out_dir, paste0(obj_name, "_seurat_res_", res, ".rds")))
    }
  
    message("\n===== All resolutions processed successfully! =====")
    message("Note: The original Seurat object contains all resolution-specific clustering results. Column names follow the format 'clusters_res{resolution}'.")
    message("Example: To view all resolution columns → colnames(seurat_object@meta.data)[grep('clusters_res', colnames(seurat_object@meta.data))]")
  
    ###--------------------------------------------------------------------- 4. Return Object with All Resolutions -------------------------------------------------------------------###
    # Save Seurat object containing all resolution-specific clustering results
    saveRDS(seurat_object, 
            file = file.path(out_dir, paste0(obj_name, "_Allres.rds")) )  

    return(seurat_object)
    
}
