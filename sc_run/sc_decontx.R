#' Single-cell RNA-seq Data Decontamination with decontX (Single Sample Only)
#'
#' This function removes environmental RNA contamination from single-sample Seurat objects using decontX.
#' It outputs contamination statistics, visualizations, filters high-contamination cells, and saves results.
#'
#' @title Single-sample scRNA-seq Decontamination (decontX)
#' @param seurat_obj Input single-sample Seurat object (required).
#' @param out_dir Output directory (required; will be created if it doesn't exist).
#' @param savecell Maximum contamination threshold for retaining cells (default: 0.2, i.e., retain cells with contamination < 0.2).
#' @param plot_width Width of visualization plots (inches, default: 8).
#' @param plot_height Height of visualization plots (inches, default: 6).
#' @return Filtered Seurat object (without contamination metadata) and saves removed cell IDs.
#' @author Custom
#' @export
#' @examples
#' \dontrun{
#' # Load preprocessed single-sample Seurat object
#' seu <- readRDS("preprocessed_single_sample_seu.rds")
#' 
#' # Run decontamination
#' seu_decont <- sc_decontx(
#'   seurat_obj = seu,
#'   out_dir = "../03.DecontX_Results",
#'   savecell = 0.2  # Retain cells with contamination < 0.2
#' )
#' }

sc_decontx <- function(seurat_obj = NULL,
                       out_dir = NULL, 
                       savecell = 0.2,
                       plot_width = 8,
                       plot_height = 6) {
    
    suppressMessages({
        library(Seurat)
        library(decontX)
        library(ggplot2)
        library(dplyr)
    })
    set.seed(1234)

    ###--------------------------------------------------------------------- 1. Parameter Validation -------------------------------------------------------------------###
    if (is.null(seurat_obj)) stop("Please provide a Seurat object (seurat_obj)!")
    if (savecell < 0 | savecell > 1) stop("Error: savecell (contamination threshold) must be between 0 and 1!")
    if (is.null(out_dir)) stop("Please specify an output directory (out_dir)!")
    
    ###--------------------------------------------------------------------- 2. Create Output Directory -------------------------------------------------------------------###
    obj_name <- deparse(substitute(seurat_obj))
    dir_name <- sub("_[^_]*$", "", obj_name)
    out_dir <- file.path("../03.Output", dir_name)
    
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
        message(paste("Automatically created output directory:", out_dir))
    }

    ###--------------------------------------------------------------------- 3. Run decontX -------------------------------------------------------------------###
    message("===== Beginning DecontX! =====")
    
    # Extract raw count matrix
    count <- GetAssayData(object = seurat_obj, assay = "RNA", layer = "counts")
    decontX_results <- decontX(count) 
    
    # Add contamination scores to metadata (temporarily for filtering)
    seurat_obj$Contamination <- decontX_results$contamination

    message("===== DecontX Completed! =====")

    ###--------------------------------------------------------------------- 4. Contamination Statistics -------------------------------------------------------------------###
    message("\n===== Contamination Statistics =====")
    # Calculate key statistics
    contam_stats <- seurat_obj$Contamination %>%
        summarise(
            Min = min(.),
            Q25 = quantile(., 0.25),
            Median = quantile(., 0.5),
            Q75 = quantile(., 0.75),
            Max = max(.),
            Mean = mean(.),
            SD = sd(.),
            Q95=(.,0.95)
        ) %>%
        round(4)
  
    # Print statistics to console
    print(contam_stats, row.names = FALSE)
  
    # Save statistics
    write.csv(contam_stats, 
              file = file.path(out_dir, paste0(obj_name, "_Contamination_Stats.csv")),
              row.names = FALSE)
    message(paste("Note: Statistics saved to:", out_dir))

    ###--------------------------------------------------------------------- 5. Contamination Visualization -------------------------------------------------------------------###
    message("\n===== Generating Contamination Visualization =====")
    contam_plot_data <- data.frame(Contamination = seurat_obj$Contamination)
  
    # Create violin + boxplot
    p_contam <- ggplot(contam_plot_data, aes(x = "", y = Contamination)) +
        geom_violin(fill = "#87CEEB", alpha = 0.7, linewidth = 0.5) +
        geom_boxplot(width = 0.2, fill = "white", linewidth = 0.8, outlier.size = 1) +
        geom_hline(yintercept = savecell, color = "red", linetype = "dashed", linewidth = 0.8) +
        annotate("text", x = 1.2, y = savecell, color = "red", 
                 label = paste("Retention Threshold =", savecell), hjust = 0) +
        labs(
            x = "", 
            y = "Contamination Ratio", 
            title = paste(obj_name, "Sample - Cell Contamination Distribution")
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 10),
            panel.grid = element_blank()
        )
  
    # Save plot
    ggsave(plot = p_contam,file = file.path(out_dir, paste0(obj_name, "_Contamination_VlnBoxPlot.pdf")), width = plot_width, height = plot_height, dpi = 300)
    message(paste("Note: Visualization saved to:", out_dir))
    
    ###--------------------------------------------------------------------- 6. Filter High-Contamination Cells -------------------------------------------------------------------###
    message("\n===== Filtering Cells by Contamination Threshold =====")
    # Identify removed cells
    removed_cell_ids <- colnames(seurat_obj)[seurat_obj$Contamination >= savecell]
    
    # Filter cells (keep contamination < savecell)
    seurat_filtered <- seurat_obj[, seurat_obj$Contamination < savecell]
    
    # Remove contamination metadata from filtered object
    seurat_filtered$Contamination <- NULL
  
    # Output filtering statistics
    message(paste("Original cell count:", ncol(seurat_obj)))
    message(paste("Retained cell count (contamination <", savecell, "):", ncol(seurat_filtered)))
    message(paste("Removed cell count:", length(removed_cell_ids)))
    message(paste("Cell retention rate:", round(ncol(seurat_filtered)/ncol(seurat_obj)*100, 2), "%"))

    ###--------------------------------------------------------------------- 7. Save Results -------------------------------------------------------------------###
    # Save removed cell IDs
    write.table(removed_cell_ids,
                file = file.path(out_dir, paste0(obj_name, "_Removed_CellIDs.txt")),
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    ###--------------------------------------------------------------------- 8. Completion Message -------------------------------------------------------------------###
    message("\n===== Decontamination Analysis Completed! =====")
    message("Output files:")
    message(paste("1. Original object with contamination data:", paste0(obj_name, "_DecontX_WithContam.rds")))
    message(paste("2. Filtered object (recommended for downstream use):", paste0(obj_name, "_DecontX_Filtered.rds")))
    message(paste("3. Raw decontX results:", paste0(obj_name, "_DecontX_RawResults.rds")))
    message(paste("4. Contamination statistics:", paste0(obj_name, "_Contamination_Stats.csv")))
    message(paste("5. Visualization:", paste0(obj_name, "_Contamination_VlnBoxPlot.pdf")))
    message(paste("6. Removed cell IDs:", paste0(obj_name, "_Removed_CellIDs.txt")))

    return(seurat_filtered)
}
