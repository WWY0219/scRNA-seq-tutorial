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

# ============================================= Load singleR Database  ===================================================
# human
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()
DICE <- DatabaseImmuneCellExpressionData() 
NHD <- NovershternHematopoieticData() 
MID <- MonacoImmuneData()
# mouse
MRD <- MouseRNAseqData() 
IGD <- ImmGenData()

# ================================================ Run SingleR  ===========================================================
sce <- seurat_obj                                        # 赋值给其他变量，避免修改原变量
sce_for_SingleR <- GetAssayData(sce, slot="data")
head(sce_for_SingleR)
clusters <- sce@meta.data$seurat_clusters
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
sce@meta.data$singleR = celltype[match(clusters,celltype$ClusterID),'celltype']
p <- DimPlot(sce, reduction = "tsne", group.by = "singleR")
p
plotScoreHeatmap(pred.hesc)

# 基于 per-cell “deltas”诊断，Delta值低，说明注释结果不是很明确
plotDeltaDistribution(pred.hesc, ncol = 3)


# ================================================ 自定义数据注释集 ===========================================================
myref <- sce2##这里为了检测，我们将参考数据集与目标数据集用同一个数据进行测试
myref$celltype <- Idents(myref)
table(Idents(myref))
# 读入参考数据集 -------
Refassay <- log1p(AverageExpression(myref, verbose = FALSE)$RNA)#求细胞的平均表达
#Ref <- textshape::column_to_rownames(Ref, loc = 1)#另一种得到参考矩阵的办法
head(Refassay)#看看表达矩阵长啥样
ref_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=Refassay))
#参考数据集需要构建成一个SingleCellExperiment对象
# if(!require(scater))BiocManager::install('scater')
library(scater)
ref_sce <- scater::logNormCounts(ref_sce)
colData(ref_sce)$Type <- colnames(Refassay)
###提取自己的单细胞矩阵##########
testdata <- GetAssayData(sce, slot="data")
pred <- SingleR(test=testdata, 
                ref=ref_sce, 
                labels=ref_sce$Type,
                #clusters = scRNA@active.ident
                )
table(pred$labels)
cellType=data.frame(seurat= sce@meta.data$seurat_clusters,
                    predict=pred$labels)#得到seurat中编号与预测标签之间的关系

lalala <- as.data.frame(table(cellType[,1:2]))
finalmap <- lalala %>% group_by(seurat) %>% top_n(n = 1, wt = Freq)#找出每种seurat_cluster注释比例最高的对应类型
finalmap <-finalmap[order(finalmap$seurat),]$predict#找到seurat中0：12的对应预测细胞类型
print(finalmap)

testname <- sce
new.cluster.ids <- as.character(finalmap)
names(new.cluster.ids) <- levels(testname)
testname <- RenameIdents(testname, new.cluster.ids)
P11 <- DimPlot(sce,label = T)
P12 <- DimPlot(testname,label = T)#比较一下测试数据与参考数据集之间有没有偏差
P11 + P12 + P9

