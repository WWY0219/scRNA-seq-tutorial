# =================================== Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")  #replace your workspace
getwd()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
set.seed(1234)

# =================================== Load scData with celltype_major  ===================================================
surat_obj <- readRDS("../03.Output/seurat_obj_harmony/seurat_obj_harmony.rds")
Idents(seurat_obj) <- "clusters_res0.5"



# =================================== Head cellmarker For celltype-major ====================================
cellmarkers <- list(
  #间充质细胞
  MSC =c("CD34","PDGFRA"),
  #免疫细胞
  immunocell= c("PTPRC"),
  T_cell = c("CD3E", "CD3D", "CD4", "CD8A"),
  B_cell = c("MS4A1", "CD19", "CD79A"),
  Monocyte_Macrophage = c("CD14", "CD68", "CSF1R", "FCGR3A", "LYZ","ITGAM"),
  Dendritic_cell = c("CD1C", "CLEC9A", "HLA-DRA", "CD83"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "NCAM1"),  # NCAM1即CD56
  Mast_cell = c("CPA3", "TPSAB1", "KIT"),
  # 基质细胞
  Fibroblast = c("COL1A1", "COL3A1", "THY1", "DCN", "FAP"),  # THY1即CD90
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),  # PECAM1即CD31
  Smooth_muscle_cell = c("ACTA2", "MYH11", "TAGLN"),
  SMC_cancercell =c("ACTG2","PRLR","SFRP4"),
  ESC =c("PLN","RGS5","SUSD2"),
  #增值细胞 
  Proliferating_cell = c("MKI67", "PCNA", "TOP2A", "CCNB1"),
  # 上皮细胞
  Epithelial_cell = c("EPCAM", "KRT8", "KRT18", "CDH1")
)

p_cellmarkers <- DotPlot(object = seurat_obj,
             features = cellmarkers,  
             cols = c("grey", "red"),
             cluster.idents = TRUE
            ) +
RotatedAxis() +
theme(panel.border = element_rect(color = "black", fill = NA),  
      panel.spacing = unit(1, "mm"),
      strip.text = element_text(margin = margin(b = 3, unit = "mm")),
      strip.placement = 'outlet',
      axis.line = element_blank()
     ) +
labs(x = "", y = "")
print(p_cellmarkers)
ggsave("../03.Output/seurat_obj_cellmarkers_res0.5.pdf", plot = p_cellmarkers, width = 20, height = 12, dpi = 300)

# =================================== Major celltype Annotation =====================================
















