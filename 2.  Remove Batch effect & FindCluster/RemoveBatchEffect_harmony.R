# ============================ 加载配置和工作环境 ============================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
list.files()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(harmony)
library(RColorBrewer)
set.seed(1234)


# ===============================================Load scData after QC-tutorial===================================================
seurat_obj <- readRDS("../03.Output/seurat_obj_merge_QC.rds")
seurat_obj[["RNA"]]=JoinLayers(seurat_obj[["RNA"]])
print(seurat_obj)
table(seurat_obj@meta.data$orig.ident)

# ===============================================Load scData after QC-tutorial===============================================
## 用细胞总 UMI 计数的中位数作为缩放因子消除细胞间测序差异
seurat_obj <- NormalizeData(seurat_obj, normalization.method ="LogNormalize", 
                            scale.factor = median(seurat_obj@meta.data$nCount_RNA))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000) 

## 计算细胞周期评分
cc.genes.updated.2019 <- cc.genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes                        
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes) 

## 回归掉不感兴趣的变量
seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo1",
                                                        "percent_ribo2","percent_mt","percent_RBC"))
## 使用HVG去跑PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj,ndims = 50) 

## Perform Harmony batch correction
seurat_obj <- seurat_obj %>% RunHarmony(
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony",    
  plot_convergence = TRUE,
  max.iter = 30,
  verbose = FALSE
)
## Save Harmony convergence plot
ElbowPlot(seurat_obj, reduction = "harmony", ndims = 50) + 
  theme_minimal() + 
  ggtitle("Elbow Plot for Harmony-Corrected Dimensions")
ggsave(file.path(out_dir, paste0(obj_name, "_Harmony_convergence.pdf")), width = 8, height = 6, dpi = 300)


# ===============================================Run sc_resolutionfinder.R===============================================
## Loading sc_resolutionfinder.R
source("sc_resolutionfinder.R")
res <- seq(0.1, 1, by = 0.2)
res
seurat_obj <- sc_resolutionfinder(seurat_obj = seurat_obj,
                                  max_dim_pca = 30,
                                  max_dim_harmony = 30,
                                  resolutions = 0.5, 
                                  top_n = 20,
                                  out_dir = "../03.Output/")
##Output profiles-seurat_obj only contains resolution data, not containing Reduction Data(UMAP/tSNE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30, 
                                reduction.name = "umap", verbose = FALSE)
seurat_obj <- RunTSNE(seurat_obj, reduction = "harmony", dims = 1:30, 
                                reduction.name = "tsne", verbose = FALSE)
qsave(seurat_obj, "../03.Output/seurat_obj_harmony.qs")


# ===============================================Find Best Resolution===============================================
table(seurat_obj@meta.data$seurat_clusters)
## clustree method
p_clustree <- clustree(seurat_obj, prefix = "cluster_res") + coord_flip()
p_clustree 
ggsave("../03.Output/clustree.pdf",width =30 ,height =30,dpi =300)

## other methods (待补充）


## Determined your resolution
levels(Idents(seurat_obj))
Idents(seurat_obj) <- "cluster_res0.5"
plot1=DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T)
plot2=DimPlot(seurat_obj,reduction = "umap",group.by = "cluster_res0.5",label = T)+NoLegend()
plot = plot1+plot2
ggsave("../03.Output/USOO_umap.pdf",width = 20, height = 15,dpi =300)


# ===============================================Cellmarkers Figures===============================================
dir.create("../03.Output/cellmarkers_fig/")

## AllcellMarkers-1
cell_markers <- list(
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
  Endothelial_cell = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),  # PECAM1即CD31
  Smooth_muscle_cell = c("ACTA2", "MYH11", "TAGLN"),
  SMC_cancercell =c("ACTG2","PRLR","SFRP4"),
  ESC =c("PLN","RGS5","SUSD2"),
  #增值细胞 
  Proliferating_cell = c("MKI67", "PCNA", "TOP2A", "CCNB1"),
  # 其他常见细胞
  Epithelial_cell = c("EPCAM", "KRT8", "KRT18", "CDH1")
)

p1 <- DotPlot(object = seurat_obj,
              features = cell_markers,  
              cols = c("grey", "red"),
              luster.idents = TRUE) +
RotatedAxis() +
theme( panel.border = element_rect(color = "black", fill = NA),  
       panel.spacing = unit(1, "mm"),
       strip.text = element_text(margin = margin(b = 3, unit = "mm")),
       strip.placement = 'outlet',
       axis.line = element_blank()
     ) +
labs(x = "", y = "")
print(p1)
ggsave("../03.Output/usoooo_cellmarker_0.5.pdf", plot = p, width = 20, height = 12, dpi = 300)

##Cell-Markers2
known_markers=list("T/NK"=c('CD3E','CD4','CD8A'),
                    B=c("CD79A","MS4A1"),
                    Plasma=c("MZB1"),
                    Epithelium=c('EPCAM',"KRT8","KRT18"),
                    Tumor=c("MKI67","FOXJ1","SOX2","SOX9"),
                    "Mono/Macro"=c("CD14","CD68","CD163"),
                    Mast=c("KIT"),
                    Fibroblast=c("COL5A2",'PDGFRB'),
                    Pericyte=c("CSPG4","RGS5"),
                    Endothelium=c("PECAM1")) 
plot3=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()

p_know <-plot2|plot3
ggsave("../03.Output/USOO_p_know.pdf", plot = p_know, width = 20, height = 12, dpi = 300)

# ===============================================Major Celltype Annotation===============================================
top10 <- seurat_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()                            # Doheatmap图
VlnPlot(seurat_obj, features = top10$gene[1:5], pt.size=0)                           #小提琴图观察基因分布
marker2 <- FindMarkers(object =USOO, ident.1 = 2)
marker_DEG <- FindMarkers(object = USOO, ident.1 = 10, ident.2 = 11,                 # 第二个cluster（对比组）
                                                      logfc.threshold = 0.25,        # 最小log2倍数变化（过滤微小差异）
                                                      min.pct = 0.1,                 # 基因在至少10%的细胞中表达（过滤低表达基因）
                                                      test.use = "wilcox",           # 统计检验方法（默认wilcox，适用于单细胞数据）
                                                      assay = "RNA" )                  # 使用的assay（默认"RNA"）


## CellAnnotationsMethods--1
meta_supp = data.frame(seurat_cluster = 0:(length(unique(seurat_obj$cluster_res0.5)) - 1), celltype = NA)
meta_supp[meta_supp$seurat_cluster %in% c(0), 'celltype'] = 'B'
meta_supp[meta_supp$seurat_cluster %in% c(4), 'celltype'] = 'Plasma'
meta_supp[meta_supp$seurat_cluster %in% c(1), 'celltype'] = 'T/NK'
meta_supp[meta_supp$seurat_cluster %in% c(5), 'celltype'] = 'Stromal'
meta_supp[meta_supp$seurat_cluster %in% c(3), 'celltype'] = 'Monocyte'
meta_supp[meta_supp$seurat_cluster %in% c(2), 'celltype'] = 'Tumor'
meta_supp[meta_supp$seurat_cluster %in% c(2), 'celltype'] = 'db-like'

for (i in 1:nrow(meta_supp)) {
  seurat_obj@meta.data[which(seurat_obj$cluster_res0.5 == meta_supp$seurat_cluster[i]), 'celltype_major'] = meta_supp$celltype[i]
}
Idents(seurat_obj) <- 'celltype_major'

## remove db-like cell
seurat_obj_filtered <- subset(seurat_obj, subset = celltype_major != "db-like")

# 看看注释情况
plot4=DimPlot(seurat_obj_filtered,group.by = "celltype_major",label = T)&NoLegend()
plot5=DotPlot(object = seurat_obj_filtered,
              features = known_markers,
              scale=T,
              group.by = "celltype_major")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot4|plot5
table(seurat_obj_filtered@meta.data$celltype_major,useNA = "always")  #If NA, plot5 will error

## save .qs
qsave(seurat_obj_filtered, file = 'seurat_obj_post_annotation.qs')


## CellAnnotationsMethods--2
new.cluster.ids <- c("0"= "Tumor cell",
                    "1"=" Endothelial_cell",
                    "2"= "Tumor cell",
                    "3"= "MSC",
                    ...)
seurat_obj <- RenameIdents(USOO, new.cluster.ids)
p_umap <- DimPlot(seurat_obj,
                  reduction = "umap",  # 若用t-SNE则改为"tsne"
                  label = TRUE,        # 显示细胞类型标签
                  label.size = 5,      # 标签字体大小
                  repel = TRUE,        # 避免标签重叠
                  cols = NULL          # 自动分配颜色（或自定义：cols = c("细胞类型1"="red", ...)）
) + 
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

print(p_umap)

# 保存图片
ggsave(
  "../03.Output/USOO_celltype_umap.pdf",
  plot = p_umap,
  width = 10,
  height = 8,
  dpi = 300
)
USOO@meta.data$celltype <- Idents(USOO)
# ============================ 假设细胞类型注释存储在 meta.data$celltype 中===================================
# 统计每个 orig.ident × celltype 的细胞数量
celltype_counts <- table(USOO$orig.ident, USOO$celltype) %>% 
  as.data.frame() %>% 
  rename(orig.ident = Var1, celltype = Var2, count = Freq)

# 计算每个 orig.ident 中各细胞类型的占比（百分比）
celltype_proportions <- celltype_counts %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100)
p_prop <- ggplot(celltype_proportions, aes(x = orig.ident, y = percentage, fill = celltype)) +
  geom_col(position = "fill") +  # 百分比堆积
  scale_y_continuous(labels = scales::percent_format()) +  # y轴显示百分比
  labs(x = "orig.ident", y = "celltype percentage", fill = "celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_prop)
ggsave("../03.Output/USOO_origident_proportion.pdf", p_prop, width = 10, height = 6, dpi = 300)

# ============================ ggplot2绘制UMAP图===================================
library(ggunchull)   # 绘制细胞群轮廓线（基于密度的凸包）
library(ggrepel)     # 避免标签重叠
library(ggplot2)     # 基础绘图包
library(dplyr)       # 数据处理包

##定义色板
cols <- c(
  # 蓝色系（3种渐变）
  colorRampPalette(c("#1f78b4", "#a6cee3"))(3),
  # 绿色系（3种渐变）
  colorRampPalette(c("#33a02c", "#b2df8a"))(3),
  # 粉色
  "#df65b0",  
  # 橙色系（3种渐变）
  colorRampPalette(c("#fdbf6f", "#ff7f00"))(3), 
  # 基础单色
  "#696969", "#d8bfd8", "#008b00", "#fb9a99", "#8b5a2b", "#da70d6",
  # 青绿色系（2种渐变）
  colorRampPalette(c("#66cdaa", "#5f9ea0"))(2), 
  # 紫色系（6种渐变）
  colorRampPalette(c("#cab2d6", "#6a3d9a"))(6),
  # 补充单色
  "#1c6597", "#cc6805", "#c6abd4", "#569395", "#dc65d8", "#f99190", "#057605","#575757"
)

### 2. 准备绘图数据（所有单细胞的metadata+UMAP坐标）
# 提取UMAP_harmony坐标并合并metadata
meta <- cbind(
  seurat_obj@meta.data,  # 含celltype的所有细胞metadata
  seurat_obj@reductions$umap_harmony@cell.embeddings  
)

# *查看UMAP坐标的实际列名（避免列名写错，比如是"UMAP_1"还是"umapharmony_1"）*
cat("UMAP坐标列名：", colnames(obj@reductions$umap_harmony@cell.embeddings), "\n")
# 把列名赋值给变量（后续用变量引用，避免硬编码错误）
umap_col1 <- colnames(obj@reductions$umap_harmony@cell.embeddings)[1]  # 第一列（x轴）
umap_col2 <- colnames(obj@reductions$umap_harmony@cell.embeddings)[2]  # 第二列（y轴）

### 为celltype分配颜色
names(cols) <- levels(factor(
  x = levels(meta$celltype), 
  levels = levels(meta$celltype), 
  ordered = TRUE
))

### 4. 计算每个celltype标签的位置（用变量引用列名，避免错误）
main_type_med <- meta %>% 
  group_by(celltype) %>% 
  summarise(
    # 用变量引用UMAP列名，确保和meta中的列名一致
    x = median(.data[[umap_col1]]) - 1,  # 标签x坐标（中位数-偏移）
    y = median(.data[[umap_col2]]) - 1,  # 标签y坐标（中位数-偏移）
    .groups = "drop"
  )

### 5. 绘制UMAP图（全程用变量引用列名，彻底避免列名错误）
umap_plot <- ggplot(
  data = meta,
  aes(x = .data[[umap_col1]], y = .data[[umap_col2]])  # 用变量引用x/y轴列名
) +
  # 2. 绘制单细胞点
  geom_point(aes(color = celltype), size = 0.2, show.legend = FALSE) +
  # 3. 隐藏坐标轴刻度
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  # 4. 应用色板
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  # 5. 添加celltype标签（用main_type_med的x/y列，而非原UMAP列名）
  geom_text_repel(
    data = main_type_med,
    aes(x = x, y = y, label = celltype, color = celltype),  # 用标签数据框的x/y列
    fontface = "bold", size = 6, show.legend = FALSE,
    box.padding = 0.5, max.overlaps = 30
  ) +
  # 6. 设置坐标轴标题
  labs(
    x = "UMAP_harmony 1",
    y = "UMAP_harmony 2"
  ) +
  # 7. 调整主题
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(hjust = 0.05, size = 12)
  )
print(umap_plot)
ggsave("../03.Output/umap.pdf", plot = umap_plot,
    width = 8,
    height = 8,
    dpi = 300,
  device = "pdf",
  limitsize = FALSE
  )
p <- p + 
    stat_unchull(
    aes(fill = celltype, color = celltype), 
    alpha = 0.1,        # 提高透明度（降低渲染负载）
    size = 0.3,         # 减细轮廓线（避免线条过重）
    lty = 1,            # 用实线（虚线计算更复杂）
    delta = 1,          # 增大平滑度（关键！降低密度计算复杂度）
    k = 5,              # 减少邻居数量（默认k=10，减小后更快）
    show.legend = FALSE
  )
p
ggsave("../03.Output/umap.pdf", plot = p,
    width = 8,
    height = 8,
    dpi = 300
  )



