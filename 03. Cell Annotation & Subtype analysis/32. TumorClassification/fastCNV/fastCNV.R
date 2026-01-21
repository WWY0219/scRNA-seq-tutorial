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
  DimPlot(seurat_obj, reduction = "umap", group.by =  "celltype") & common_theme

## boxplot
ggplot(FetchData(seurat_obj, vars = c("celltype", "cnv_fraction")), 
       aes(celltype, cnv_fraction, fill = celltype)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

## 按染色体臂划分的CNV
library(scales)
FeaturePlot(seurat_obj, features = "20.p_CNV")  +
  scale_color_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1), 
                       rescaler = function(x, to = c(0, 1), from = NULL) {
                         rescale_mid(x, to = to, mid = 0)
                       }) +
  common_theme |
FeaturePlot(seurat_obj, features = "X.q_CNV") +
  scale_color_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1), 
                       rescaler = function(x, to = c(0, 1), from = NULL) {
                         rescale_mid(x, to = to, mid = 0)
                       }) +
  common_theme


# ====================================================== CNV classification =================================================== 
seurat_obj <- CNVClassification(seurat_obj)
 
DimPlot(seurat_obj, group.by = "20.p_CNV_classification") &
  scale_color_manual(values = c(gain = "red", no_alteration = "grey", loss = "blue")) &
  common_theme |
DimPlot(seurat_obj, group.by = "X.q_CNV_classification") &
  scale_color_manual(values = c(gain = "red", no_alteration = "grey", loss = "blue")) &
  common_theme


# ====================================================== CNV clusters =================================================== 
DimPlot(seurat_obj, group.by = "cnv_clusters") + common_theme

library(SeuratObject)
ggplot(FetchData(seurat_obj, vars = c("cnv_clusters", "celltype")), aes(annot, fill = cnv_clusters)) +
  geom_bar(position = "fill") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))

# ====================================================== CNV tree =================================================== 
tree_data <- CNVTree(seurat_obj, 
                     values = "calls", 
                     cnv_thresh = 0.09, 
                     healthyClusters = "1")







