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
library(infercnv)
library(CopyKAT)
library(dplyr)
library(ggsignif)          
library(circlize)
library(ComplexHeatmap)
set.seed(1234)

# =================================== Load scData with celltype_major  ===================================================
surat_obj <- qread("../seurat_obj.qs")
Idents(seurat_obj) <- "clusters_res0.5"
##绘制样本和cluster的umap图


# =================================== Run CopyKAT  ===================================================
copykat.test <- copykat(rawmat=expr_matrix,   # 表达矩阵
                        id.type="S",          # symbol or ensembl
                        ngene.chr=5,            # 对应细胞的基因至少在每个染色体上包含5个基因，这个细胞才会被包咯iu
                        win.size=25,     # 评分滑窗的基因数
                        KS.cut=0.1,      # 控制敏感度
                        sam.name="test",             # sample name
                        distance="euclidean",        # cluster-method
                        norm.cell.names="",          # normal cell for control
                        output.seg="FLASE", # 是否输出seg文件
                        plot.genes="TRUE",          # 是否在热图中绘制基因
                        genome="hg20",              # hg20--hg38，mm10--mm10
                        n.cores=4)                  # 开启4核心的并行计算

# =================================== CopyKAT-data Operate  ===================================================
## CNV标签
pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[which(pred.test$copykat.pred %in% c("aneuploid","diploid")),]  # keep defined cells
head(pred.test)
table(pred.test$copykat.pred,useNA = "always"))

CNA.test <- data.frame(copykat.test$CNAmat)
CNA.test[1:4,1:4]

##add pred.test into metadata
cnv_map <- pred.test[, c("cell.names", "copykat.pred")]
rownames(cnv_map) <- cnv_map$cell.names
cnv_map <- cnv_map[, "copykat.pred", drop = FALSE]

usoo_barcodes <- colnames(USOO)  # !!seurat_obj_all name!!
cnv_annot <- data.frame(
  copykat.pred = rep(NA, length(usoo_barcodes)),
  row.names = usoo_barcodes
)
matched_barcodes <- intersect(rownames(cnv_map), usoo_barcodes)
cnv_annot[matched_barcodes, "copykat.pred"] <- cnv_map[matched_barcodes, "copykat.pred"]

USOO <- AddMetaData(USOO, metadata = cnv_annot)
table(USOO$copykat.pred, useNA = "always")

# =================================== CopyKAT-data Visulization  ===================================================
DimPlot(object = USOO, reduction = "umap", group.by = "copykat.pred",  
        cols = c(
          "diploid" = "#4DAF4A",    # 绿色（二倍体/正常细胞）
          "aneuploid" = "#E41A1C",  # 红色（非整倍体/肿瘤细胞）
          "NA" = "lightgray"        # 浅灰色（无CNV结果的细胞）
        ),
        na.value = "lightgray",     # 强制NA为浅灰色（兜底）
        label = TRUE,               # 显示分类标签（diploid/aneuploid）
        label.size = 5,             # 标签字体大小
        repel = TRUE,               # 标签不重叠
        pt.size = 0.8,              # 点的大小（适配单细胞数量）
        shuffle = TRUE              # 打乱点的绘制顺序（避免NA层覆盖有标签的点）
       ) +
ggtitle("CNV Classification (copykat)") +
xlab("UMAP_1") + ylab("UMAP_2") +
theme_bw() +
theme(
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中/加粗
  axis.title = element_text(size = 14),                              # 坐标轴标题大小
  axis.text = element_text(size = 12),                               # 坐标轴刻度大小
  legend.title = element_blank(),                                    # 隐藏图例标题
  legend.text = element_text(size = 12),                             # 图例文字大小
  panel.grid = element_blank()                                       # 去除网格线
  )







