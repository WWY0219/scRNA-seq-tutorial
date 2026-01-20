# =========================================================== Prepare Environment ===================================================
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
library(patchwork)
library(RColorBrewer)
library(RcppML)
library(readr)
library(irGSEA)
library(ggsci)
library(doMC)
library(ggunchull)
set.seed(1234)
dir.create("./irGSVA/")

# ================================================ Load scData with celltype_major  ==================================================
seurat_obj <- qread("./03.Output/seurat_obj_harmony/seurat_obj_harmony.qs")
DimPlot(seurat_obj,  reduction = 'umap', 
        group.by="clusters_res0.5", label = TRUE, pt.size = 0.5) + NoLegend()

# ===============================================irGSEA通过内置的msigdbr包进行打分======================================================
## Hallmark基因集打分
seurat_obj <- irGSEA.score(object = seurat_obj,
                           assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 1,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T, 
                           species = "Homo sapiens", category = "H",  
                           subcategory = NULL, geneid = "symbol",
                           method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')
Seurat::Assays(seurat_obj)

## KEGG基因集打分 
seurat_obj <- irGSEA.score(object = seurat_obj,
                           assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 1,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T, 
                           species = "Homo sapiens", category = "C2",  
                           subcategory = "CP:KEGG", geneid = "symbol",
                           method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')

## GO-BP基因集打分 
seurat_obj <- irGSEA.score(object = seurat_obj,
                           assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 1,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T, 
                           species = "Homo sapiens", category = "C5",  
                           subcategory = "GO:BP", geneid = "symbol",
                           method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')


