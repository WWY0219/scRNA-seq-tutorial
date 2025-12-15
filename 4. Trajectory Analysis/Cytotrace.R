# ============================== Prepare Environment ===============================
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
library(ggpubr)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(scRNAtoolVis)
library(CytoTRACE2)
library(monocle3)
set.seed(1234)
# 查看工作路径下的文件
list.files()

# ============================== Load scData for TraceAnalysis =========================
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)

clusterCornerAxes(object = seurat_obj,
                  reduction = 'umap',
                  clusterCol = "subcelltype",  #replace your explore subtype
                  noSplit = T,
                  cellLabel = T,
                  cellLabelSize = 5)
DimPlot(seurat_obj, redcution="umap", goup.by="subcelltype",pt.size=2 )
Idents(seurat_obj)<-"subcelltype"
levels(Idents(seurat_obj))

# ============================== Run CytoTRACE2 ========================================
seurat_obj_cytotrace <- cytotrace2(seurat_obj,
                        is_seurat = TRUE,
                        slot_type = "counts",
                        species = "human",
                        seed = 1234)
## replace UMAP/TSNE
seurat_obj_cytotrace@reductions[["umap"]] <- seurat_obj@reductions[["tsne"]] 

# ============================== Plot CytoTRACE2 ========================================
## making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame( phenotype=seurat_obj_cytotrace@meta.data$subcelltype) %>% set_rownames(., colnames(seurat_obj_cytotrace))
## cytoTRACE2 plot
plots <- plotData(cytotrace2_result = seurat_obj_cytotrace, 
                  annotation = annotation, 
                  is_seurat = TRUE,
                  seed = 1234)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
cyto_fig  <- (p1+p2+p3+p4) + plot_layout(ncol = 2)
cyto_fig
ggsave("../03.Output/cytoTRACE2_primaryfig.pdf", plot = cyto_fig, width = 10, height = 10, dpi =300)

## Other Figures 
cyto_featureplot <- FeaturePlot(seurat_obj_cytotrace, "CytoTRACE2_Relative", pt.size = 1.5) + 
  scale_colour_gradientn(colours = (c("#9E0142", "#F46D43", "#FEE08B", 
                                      "#E6F598", "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)
print(cyto_featureplot)
ggsave("../03.Output/cytoTRACE2_featureplot.pdf", plot = cyto_featureplot, width = 10, height = 10, dpi =300)

## Other Fig
p <- ggboxplot(result_LM@meta.data, x="subcelltype", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",                               #轮廓颜色
                fill="subcelltype",                            #填充
                palette = "npg",
                xlab = F,                                      #不显示x轴的标签
                bxp.errorbar=T,                                #显示误差条
                bxp.errorbar.width=0.5,                        #误差条大小
                size=1,                                        #箱型图边线的粗细
                outlier.shape=NA,                              #不显示outlier
                legend = "right")                               
###指定组比较(replace your group)
my_comparisons <- list(c("SMC_C1", "SMC_C2"), c("SMC_C1", "SMC_C3"),c("SMC_C2", "SMC_C3"))
p_box <- p + tat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")
ggsave("../03.Output/cytoTRACE2_boxplot.pdf", plot = p_box, width = 10, height = 10, dpi =300)

# ============================== Run Monocle3 ========================================
Idents(seurat_obj) <- seurat_obj$subcelltype
levels(Idents(seurat_obj))                  
# 根据CytoTRACE2/先验知识进行重设等级
seurat_obj$celltype <- factor(seurat_obj$subcelltype,levels = c("FC_C0","FC_C1","FC_C2","SMC_C1",
                                                   "SMC_C0","SMC_C2"))
Idents(seurat_obj) <- seurat_obj$subcelltype

## 提取数据
expression_matrix <- GetAssayData(seurat_obj, assay = 'RNA',slot = 'counts')
cell_metadata <- seurat_obj@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
##
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# 归一化/预处理数据
# 这一步使用的PCA分析，dim数代表纳入的PCA数量
cds <- preprocess_cds(cds, num_dim = 25)
# 这个函数用于确认设定的dim数是否足够代表主要变异
plot_pc_variance_explained(cds)

# 可选(去批次处理)
#cds <- align_cds(cds, num_dim = 100, alignment_group = "GSE_num")

# 降维聚类，可选择UMAP、PCA或者TSNE
cds <- reduce_dimension(cds,reduction_method='UMAP',preprocess_method = 'PCA')
plot_cells(cds, label_groups_by_cluster=T ,  color_cells_by = "subcelltype")
