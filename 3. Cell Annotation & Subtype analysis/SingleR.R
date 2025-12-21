# =================================== Prepare Environment ===================================================
## Global Settings
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")  # replace your workspace
getwd()
## Loading Packages
library(qs)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(SingleR)
set.seed(1234)
dir.create("../03.Output/SingleR/")
load(x.Rdata)
# =================================== Load scData with celltype_single  ===================================================
seurat_obj <- readRDS("../03.Output/seurat_obj_harmony/seurat_obj_harmony.rds")
Idents(seurat_obj) <- "clusters_res0.5"
DimPlot(seurat_obj,  reduction = 'umap', 
        group.by="clusters_res0.5", label = TRUE, pt.size = 0.5) + NoLegend()

# human
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
DICE <- DatabaseImmuneCellExpressionData() 
NHD <- NovershternHematopoieticData() 
MID <- MonacoImmuneData()
# mouse
MRD <- MouseRNAseqData() 
IGD <- ImmGenData()

sce1 <- sce ###赋值给其他变量，避免修改原变量。
sce_for_SingleR <- GetAssayData(sce1, slot="data")
sce_for_SingleR
clusters <- sce1@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, 
                     ref = hpca.se, 
                     labels = hpca.se$label.main,
                     method = "cluster", 
                     clusters = clusters, 
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")

table(pred.hesc$labels)
MedBioInfoCloud: table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 
sce1@meta.data$singleR = celltype[match(clusters,celltype$ClusterID),'celltype']
P9 <- DimPlot(sce1, reduction = "tsne", group.by = "singleR")
P9
