# =======================Methods 1=================================
cell_major_colors <- c("Tumor"= "#D1352B",
                       "MSC" = "#D2EBC8",
                       "Fibroblast" = "#7DBFA7",
                       "SMC"="#EE934E"，
                       "T/NK"="#3C77AF",
                       "Mono/Macro" = "#AECDE1",
                       "Endothelial" = "#8FA4AE",
                       "Epithelial" = "#BBDD78",
                       "B" = "#F5D2A8")
DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype_major",
  col = cell_major_colors,
  pt.size = 1.2,        
  alpha = 0.8,          
  label = TRUE,         
  label.size = 5,       
  repel = TRUE         
) + 
  ggtitle("UMAP of Major Cell Types") +  
  theme(plot.title = element_text(hjust = 0.5))  



# =======================Methods 2=================================
cell_major_colors <- c("#D1352B", "#D2EBC8", "#7DBFA7", "#EE934E", "#3C77AF",
                       "#AECDE1", "#8FA4AE", "#BBDD78", "#F5D2A8", "#9B5B33",
                      "#B383B9", "EE934E")
names(cell_major_colors) <- levels(factor(
  x = levels(seurat_obj$celltype), 
  levels = levels(seurat_obj$celltype), 
  ordered = TRUE
))
cell_major <- seurat_obj@meta.data$celltype_major
## ！！按照细胞分群来分配颜色！！
names(cell_major_colors) <- cell_major



# ================================= Colors ===============================================
## 张泽民 60 配色
AcaZZM60colors <- c('#4b6aa8', '#3ca0cf', '#c376a7', '#ad98c3', '#cea5c7', '#53738c', '#a5a9b0', '#a78982', '#696a6c', '#92699e',
                    '#d69971', '#df5734', '#6c408e', '#ac6894', '#d4c2db', '#537eb7', '#83ab8e', '#ece399', '#405993', '#cc7f73',
                    '#b95055', '#d5bb72', '#bc9a7f', '#e0cfda', '#d8a0c0', '#e6b884', '#b05545', '#d69a55', '#64a776', '#cbdaa9',
                    '#efd2c9', '#da6f6d', '#ebb1a4', '#a44e89', '#a9c2cb', '#b85292', '#6d6fa0', '#8d689d', '#c8c7e1', '#d25774',
                    '#c49abc', '#927c9a', '#3674a2', '#9f8d89', '#72567a', '#63a3b8', '#c4daec', '#61bada', '#b7deea', '#e29eaf',
                    '#4490c4', '#e6e2a3', '#de8b36', '#c4612f', '#9a70a8', '#76a2be', '#408444', '#c6adb0', '#9d3b62', '#2d3462')




