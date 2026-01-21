# ============================================ Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(qs) 
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(stringr)
library(ggsci)
library(SCOP)
library(pheatmap)
library(ggrepel)
library(fastCNV)
set.seed(1234)
list.files()
dir.create("./fastCNV/")


# ====================================================== Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")

# ====================================================== Load Data ===================================================
seurat_obj <- fastCNV(seuratObj = seurat_obj, 
                      sampleName = "seurat_obj", 
                      referenceVar = "celltype", 
                      referenceLabel = c("T/NK", "Myeloid", "B","Mast", "Plasma"), 
                      printPlot = T)

# ====================================================== CNV Fraction =================================================== 
common_theme <- theme(
  plot.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6)
)
FeaturePlot(seurat_obj, features = "cnv_fraction", reduction = "umap" ) & common_theme |
  DimPlot(scColon1, reduction = "umap", group.by =  "celltype") & common_theme
