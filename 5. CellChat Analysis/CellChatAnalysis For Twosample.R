# ============================== Prepare Environment ===============================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(qs)
library(CellChat)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(SCP)

set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("../03.Output/")


# ================================================================ Load scData for cellchat ==========================================================
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)
DimPlot(seurat_obj, pt.size = 0.8,group.by = "celltype_major",label = T)
table(seurat_obj@meta.data$celltype_major)

# 自己分析的时候记得改名，这里是scRNA_left
data.input <- GetAssayData(scRNA_left, slot = 'data') # normalized data matrix
meta <- scRNA_left@meta.data[,c("orig.ident","celltype_major")]
colnames(meta) <-  c("samples","labels")
table(meta$labels)
table(meta$labels)
identical(rownames(meta),colnames(data.input))
## 根据研究情况进行细胞排序
celltype_order <- c(
  "Tumor", "Fibroblast", "SMC", "T/NK", 
  "LEC", "Monocyte", "Mastcell", "Neutrophil", 
  "Endothelial","Epithelial")
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)
# 根据 meta$labels 的顺序进行排序
ordered_indices <- order(meta$labels)
# 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))

# 构建cellchat
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human  
CellChatDB.use <- subsetDB(CellChatDB)
#CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole databa
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net_left <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_left <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
## !!重复以上操作对第二个数据集！！
## merge Cellchat_obj
object.list <- list(left=cellchat_left, right=cellchat_right)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

# ================================================================ Load scData for cellchat ==========================================================
## ------------
gg1 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(1,2),measure = "weight")
gg1 + gg2

gg3 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg4 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg3 + gg4






















