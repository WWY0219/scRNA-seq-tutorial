#' DoubleFinder-based Doublet Removal with Preprocessing
#'
#' This function performs data preprocessing, detects and removes doublets using DoubletFinder,
#' records removed doublet cell IDs, and saves results.
#'
#' @title Doublet Removal with Preprocessing and Result Saving
#' @description
#' Steps: 1. Create output directory; 2. Data preprocessing; 3. UMAP dimensionality reduction;
#' 4. Cell clustering; 5. Doublet detection with DoubletFinder; 6. Visualize doublet results;
#' 7. Filter singlet cells and record removed doublet IDs; 8. Save and return clean Seurat object.
#' @param seurat_obj Input Seurat object (required).
#' @param max.dim Number of PCA dimensions used for UMAP and clustering.
#' @param max.pcs Number of PCA dimensions used for DoubletFinder.
#' @param res Clustering resolution for FindClusters.
#' @param dbrate Expected doublet rate per cell.
#' @param pN pN parameter for DoubletFinder.
#' @param methods Normalization method, "SCT" or "LogNormalize".
#' @param width Width of output plots.
#' @param height Height of output plots.
#' @param out_dir Path to main output directory.
#' @return A clean Seurat object containing only singlet cells.
#' @author WWY
#' @export

sc_doublefinder <- function(seurat_obj = NULL,
                            max.dim = NULL,
                            max.pcs = NULL,
                            res = 1,
                            dbrate = 8 * 1e-6,
                            pN = 0.25,
                            methods = "SCT",
                            width = 8,
                            height = 6,
                            out_dir = NULL) {

    ## Load required libraries
    suppressMessages({
      library(DoubletFinder)
      library(Seurat)
      library(ggplot2)
      library(dplyr)
      library(qs)
    })

    set.seed(1234)

    seurat_temp <- seurat_obj

    ###-------------------------------- 0. Parameter Check -------------------------------------------###
    if (is.null(seurat_obj)) stop("Error: 'seurat_obj' (input Seurat object) is required!")
    if (is.null(max.dim)) stop("Error: 'max.dim' (PCA dimensions for UMAP) is required!")
    if (is.null(max.pcs)) stop("Error: 'max.pcs' (PCA dimensions for DoubletFinder) is required!")
    if (is.null(out_dir)) stop("Error: 'out_dir' (output directory path) is required!")
    methods <- match.arg(methods, choices = c("SCT", "LogNormalize"))

    ###-------------------------------- 1. Create Output Directory -----------------------------------###
    obj_name <- deparse(substitute(seurat_obj))
    dir_name <- sub("_[^_]*$", "", obj_name)
    output_dir <- file.path(out_dir, dir_name)

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        message(paste("Output directory created automatically:", output_dir))
    }

    ###-------------------------------- 2. Normalization -> PCA --------------------------------------###
    if (methods == "SCT") {
        message("Starting SCTransform -> PCA...")
        DefaultAssay(seurat_obj) <- "RNA"
        seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay.name = "SCT",
                                  variable.features.n = 3000, vst.flavor = "v2", verbose = FALSE)
        DefaultAssay(seurat_obj) <- "SCT"
        seurat_obj <- RunPCA(seurat_obj, assay = "SCT",
                             features = VariableFeatures(seurat_obj), verbose = FALSE)
        use_sct <- TRUE
    } else {
        message("Starting LogNormalize -> FindVariableFeatures -> ScaleData -> PCA...")
        DefaultAssay(seurat_obj) <- "RNA"
        seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                                    scale.factor = median(seurat_obj@meta.data$nCount_RNA))
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
        seurat_obj <- ScaleData(seurat_obj)
        seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
        use_sct <- FALSE
    }

    message("Normalization method: ", methods, " | DoubletFinder sct = ", use_sct)

    ###-------------------------------- 3. UMAP Dimensionality Reduction ------------------------------###
    message(paste("Running UMAP based on PCA dimensions 1:", max.dim, "..."))
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:max.dim,
                          reduction.name = "umap_naive", verbose = FALSE)

    ###-------------------------------- 4. Cell Clustering --------------------------------------------###
    message("Starting Cell Clustering (resolution=", res, ")...")
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:max.dim)
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    message("Clustering is completed, obtained ", length(unique(seurat_obj$seurat_clusters)), " clusters")

    ###-------------------------------- 5. Doublet Detection with DoubletFinder -----------------------###
    # Identify optimal pK value
    message("Calculating optimal pK value...")
    sweep.res.list <- paramSweep(seurat_obj, PCs = 1:max.pcs, sct = use_sct)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    opt_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    message("Optimal pK value: ", opt_pK)

    # Estimate homotypic doublet proportion
    homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)
    message("Homotypic doublet proportion: ", round(homotypic.prop, 3))

    doublet.rate <- ncol(seurat_obj) * dbrate
    nExp_poi <- round(doublet.rate * nrow(seurat_obj@meta.data))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    ## Run DoubletFinder (High confidence)
    message("Starting DoubletFinder - High Confidence...")
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:max.pcs, pN = pN, pK = opt_pK,
                                nExp = nExp_poi.adj, reuse.pANN = NULL, sct = use_sct)

    pANN_col <- paste0("pANN_", pN, "_", opt_pK, "_", nExp_poi.adj)

    ## Run DoubletFinder (Low confidence)
    message("Starting DoubletFinder - Low Confidence...")
    seurat_obj <- doubletFinder(seurat_obj, PCs = 1:max.pcs, pN = pN, pK = opt_pK,
                                nExp = nExp_poi, reuse.pANN = pANN_col, sct = use_sct)

    ###-------------------------------- 6. Visualization of DoubletFinder Results ---------------------###
    ## UMAP plot for high confidence doublets
    df_col_high <- paste0("DF.classifications_", pN, "_", opt_pK, "_", nExp_poi.adj)

    p_double.adj <- DimPlot(seurat_obj, reduction = "umap_naive",
                            group.by = df_col_high, raster = FALSE) +
        ggtitle("High Confidence Doublets vs Singlets")

    print(p_double.adj)

    ggsave(filename = file.path(output_dir, paste0("Doublets_HighConfidence_", dir_name, ".pdf")),
           plot = p_double.adj, width = width, height = height, dpi = 300)

    message("High confidence doublet UMAP plot saved")

    ## UMAP plot for low confidence doublets
    df_col_low <- paste0("DF.classifications_", pN, "_", opt_pK, "_", nExp_poi)

    p_double <- DimPlot(seurat_obj, reduction = "umap_naive",
                        group.by = df_col_low, raster = FALSE) +
        ggtitle("Low Confidence Doublets vs Singlets")

    print(p_double)

    ggsave(filename = file.path(output_dir, paste0("Doublets_LowConfidence_", dir_name, ".pdf")),
           plot = p_double, width = width, height = height, dpi = 300)

    message("Low confidence doublet UMAP plot saved")

    ## Violin plot for pANN (high confidence)
    pann_feature.adj <- paste0("pANN_", pN, "_", opt_pK, "_", nExp_poi.adj)

    p_pann.adj <- VlnPlot(seurat_obj, features = pann_feature.adj, pt.size = 0) +
        ggtitle("pANN Distribution (High Confidence)")

    print(p_pann.adj)

    ggsave(filename = file.path(output_dir, paste0("pANN_HighConfidence_", dir_name, ".pdf")),
           plot = p_pann.adj, width = width, height = height, dpi = 300)

    message("High confidence pANN violin plot saved")

    ## Violin plot for pANN (low confidence)
    pann_feature <- paste0("pANN_", pN, "_", opt_pK, "_", nExp_poi)

    p_pann <- VlnPlot(seurat_obj, features = pann_feature, pt.size = 0) +
        ggtitle("pANN Distribution (Low Confidence)")

    print(p_pann)

    ggsave(filename = file.path(output_dir, paste0("pANN_LowConfidence_", dir_name, ".pdf")),
           plot = p_pann, width = width, height = height, dpi = 300)

    message("Low confidence pANN violin plot saved")

    ###-------------------------------- 7. Select Singlet Cells --------------------------------------###
    seurat_obj@meta.data[["DF_hi.lo"]] <- seurat_obj@meta.data[[df_col_low]]

    seurat_obj@meta.data[["DF_hi.lo"]][
        seurat_obj@meta.data[[df_col_low]] == "Doublet" &
        seurat_obj@meta.data[[df_col_high]] == "Singlet"
    ] <- "Doublet-Low Confidence"

    seurat_obj@meta.data[["DF_hi.lo"]][
        seurat_obj@meta.data[[df_col_low]] == "Doublet" &
        seurat_obj@meta.data[[df_col_high]] == "Doublet"
    ] <- "Doublet-High Confidence"

    cat("Doublet confidence distribution:\n")
    print(table(seurat_obj@meta.data[["DF_hi.lo"]]))

    # UMAP plot for doublet confidence
    p_confidence <- DimPlot(seurat_obj, reduction = "umap_naive", group.by = "DF_hi.lo",
                            cols = c("Singlet" = "skyblue",
                                     "Doublet-High Confidence" = "red",
                                     "Doublet-Low Confidence" = "gold"),
                            label = FALSE) +
        ggtitle("Doublet Confidence Classification")

    print(p_confidence)

    ggsave(filename = file.path(output_dir, paste0("Doublet_Confidence_", dir_name, ".pdf")),
           plot = p_confidence, width = width, height = height, dpi = 300)

    message("Doublet confidence UMAP plot saved")

    ## Extract singlet cells
    seurat_singlet <- subset(seurat_obj, subset = DF_hi.lo == "Singlet")

    ## Record removed doublet cell IDs
    all_cell_ids <- colnames(seurat_obj)
    singlet_cell_ids <- colnames(seurat_singlet)
    removed_doublet_ids <- setdiff(all_cell_ids, singlet_cell_ids)

    write.table(data.frame(removed_doublet_ids = removed_doublet_ids),
                file.path(output_dir, paste0(dir_name, "_removed_doublet_ids.txt")),
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    message("Removed doublet cell IDs saved to: ",
            file.path(output_dir, paste0(dir_name, "_removed_doublet_ids.txt")))

    ## Subset original Seurat object to retain only singlets
    seurat_temp <- subset(seurat_temp, cells = singlet_cell_ids)

    message("Doublet removal completed! Retained cells: ", ncol(seurat_singlet),
            ", Original cells: ", ncol(seurat_obj))

    ###-------------------------------- 8. Save Clean Seurat Object ----------------------------------###
    qsave(seurat_temp, file = file.path(output_dir, paste0(dir_name, "_clean.qs")))
    message("Clean Seurat object saved to: ",
            file.path(output_dir, paste0(dir_name, "_clean.qs")))

    ###-------------------------------- 9. Final Message ----------------------------------------------###
    message("\n===== All processes completed! Results saved to: ", output_dir, " =====")

    ###-------------------------------- 10. Return Clean Seurat Object --------------------------------###
    return(seurat_temp)
}
