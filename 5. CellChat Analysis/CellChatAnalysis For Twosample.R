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
library(ComplexHeatmap)
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

# ================================================================ CellChat Visulization ==========================================================
## ------------柱状图---------------
gg1 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, 
                           group = c(1,2),measure = "weight")
gg1 + gg2

gg3 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg4 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg3 + gg4

## ------------circle---------------
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

## ------------heatmap---------------
gg1 <- netVisual_heatmap(cellchat)
gg1
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

## ------------violin------------------
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("left", "right")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

## ------------violin------------------
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

patchwork::wrap_plots(plots = gg)

## 可视化从left样本到right样本的发出与接收信号差异变化
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T/NK", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tumor", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))

## -------------------------------
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
ht1
ht2

## -----------------------
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:7),  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = 4, 
                        targets.use = c(1),  
                        comparison = c(1, 2), 
                        max.dataset = 2, 
                        title.name = "Increased signaling in Left", angle.x = 45, 
                        remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, 
                        sources.use = 4, 
                        targets.use = c(1),  
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = "Decreased signaling in right", angle.x = 45, 
                        remove.isolate = T)
gg1 + gg2


## ------------------Chord diagram-------------------
par(mfrow = c(1,2), xpd=TRUE, mar = c(1,1,2,1))
netVisual_chord_gene(object.list[[2]], 
                     sources.use = 4, 
                     targets.use = c(1), 
                     slot.name = 'net', 
                     net = net.up, 
                     lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], 
                     sources.use = 4, 
                     targets.use = c(1), 
                     slot.name = 'net', 
                     net = net.down, 
                     lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

## ----------------------词云图----------------------
par(mfrow = c(1,2), xpd=TRUE, mar = c(1,1,2,1))
library(wordcloud)
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

## --------------------------------------------------
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

## --------------------------------------------------
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

# Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
ht[[1]]

## -----------------------------Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# Chord diagram 的另外一种形式
# group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
# names(group.cellType) <- levels(object.list[[1]]@idents)
# pathways.show <- c("CXCL") 
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
# }

par(mfrow = c(1, 3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 4, targets.use = c(5:6), lab.cex = 0.5, title.name = paste0("Signaling from Tm - ", names(object.list)[i]))
}

# compare all the interactions sending from fibroblast to inflamatory immune cells
# par(mfrow = c(1, 2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_chord_gene(object.list[[i]], sources.use = c(1,2, 3, 4), targets.use = c(8,10),  title.name = paste0("Signaling received by Inflam.DC and .TC - ", names(object.list)[i]), legend.pos.x = 10)
# }













