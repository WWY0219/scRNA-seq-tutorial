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


# ============================== Prepare scData for cellchat =========================
data.input <- GetAssayData(seurat_obj, layer = 'data')           # normalized data matrix
meta <- seurat_obj@meta.data[,c("orig.ident","celltype_major")]  # your cellType
colnames(meta) <-  c("samples","labels")
table(meta$labels)
identical(rownames(meta),colnames(data.input))                   # 严格判断两个向量是否一致
celltype_order <- c(
  "Tumor", "Fibroblast", "SMC", "T/NK", 
  "LEC", "Monocyte", "Mastcell", "Neutrophil", 
  "Endothelial","Epithelial")
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)
ordered_indices <- order(meta$labels)
## 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))

# ============================== Create cellchat_obj =========================
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels")
levels(cellchat@idents)
## Loading Database
CellChatDB <- CellChatDB.human          # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)        # Display database category-All
dplyr::glimpse(CellChatDB$interaction)  # Display cell
## subset use category
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use           # load database into cellchat_obj
## expressiondata preoperate
cellchat <- subsetData(cellchat,features = NULL)
cellchat <- identifyOverExpressedGenes(cellchat)  # find high-expression genes
# 默认情况下,cellchat使用object@data.signaling进行网络推断
# 同时也提供了smoothData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
cellchat <- smoothData(cellchat, adj = PPI.human)
# 该分析的关键参数是类型，即计算每个细胞组的平均基因表达的方法。默认情况下，type = “triMean”，产生较少但更强的交互。当设置 type = “truncatedMean” 时，应为trim分配一个值，从而产生更多交互。请详细检查上述计算每个细胞组平均基因表达的方法。
# 使用的是投射到PPI网络的模式时候需要用FALSE。如果使用了raw data就需要设置为TRUE
cellchat <- computeCommunProb(cellchat, 
                              type = "triMean",
                              raw.use = FALSE)   # use rawdata please choose TRUE
# 如果所研究的信号没有被测到，可以采用如下函数进行探查，trim设为0.1或者0.05
# computeAveExpr(cellchat, features = c("CXCL12","CXCR4"),type =  "truncatedMean",trim = 0.1)
# 如果发现修改参数之后所研究的信号被测到了，那就修改代码如下
# cellchat <- computeCommunProb(cellchat, type =  "truncatedMean",trim = 0.1,raw.use = FALSE) 

## min.cells是设置阈值，最小是需要10个细胞参与通讯推断(可以自定义)
cellchat <- filterCommunication(cellchat, min.cells = 10)
# CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平上的通信概率。 
# NB:推断的每个配体-受体对的细胞间通信网络和每个信号通路分别存储在槽'net'和'netP'中。
cellchat <- computeCommunProbPathway(cellchat)

## subset communications
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)
DT::datatable(df.net)

## save cellchat_obj
qsave(cellchat,"cellchat.qs")
save(df.net,file = "df.net.Rdata")
write.csv(df.net,"df.net.csv")
# 计算聚合细胞-细胞通信网络
# 互作网络整合,可以设置soure和target，不设置就是默认全部
cellchat <- aggregateNet(cellchat)

# ================================================================ CellChat Visulization ==========================================================
## ---------------CircleFig----------------
groupSize <- as.numeric(table(cellchat@idents)) 
par(mfrow = c(1,3), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 targets.use = 'your target celltype')   #please use your targeted celltype

## ---------------互作数量与重要性图---------------------
#指定顺序和指定颜色
celltype_order <- c(
  "Tumor", "Fibroblast", "SMC", "T/NK", "LEC", 
 "Monocyte", 
  "Mastcell", 
  "Neutrophil", 
  "Endothelial",
  "Epithelial"
)
mat <- as.data.frame(cellchat@net$weight)
mat <- mat[celltype_order,]#行排序
mat <- mat[,celltype_order] %>% as.matrix()
# 生成颜色向量（例如使用彩虹色）
color.use <- rainbow(nrow(mat))
# 将颜色向量命名为矩阵的行名
names(color.use) <- rownames(mat)
# 如果图片显示不全,需要考虑是不是重新设置mfrow
par(mfrow = c(5,4), xpd=TRUE,mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}


## --------------------层次图--------------------------
cellchat@netP$pathways
levels(cellchat@idents) 
pathways.show <- "PDGF"  # pathways.show <- df.net$pathway_name#计算到的所有通路
stopifnot(pathways.show %in% cellchat@netP$pathways)
# Hierarchy plot
# vertex.receiver定义层次图的左边细胞
vertex.receiver = seq(1:7)  # 两侧细胞加起来等于总细胞数
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  vertex.receiver = vertex.receiver,
  layout = "hierarchy",
  top = 0.5,
  remove.isolate = FALSE,
  vertex.label.cex = 0.7,
  edge.width.max = 6
)


## -------------------左侧显示只当细胞的的通讯，右侧显示剩余细胞的通讯情况--------------------------
### circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
### Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

## -------------------分组弦图--------------------------
levels(cellchat@idents)
group.cellType <- c(rep("T/NK", 9), "Tumor","Fibroblast","SMC","LEC",
                    "Monocyte","Mastcell","Neutrophil","Endothelial","Epithelial" )  #T/NK细胞含有亚型写法，否则不需重复
names(group.cellType) <- levels(cellchat@idents)
print(head(group.cellType, 10))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,              #分组展示的顺序
                     title.name = paste0(pathways.show, " signaling network"))

## -------------------Heatmap--------------------------
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")



# ================================================================ 配受体展示 ==========================================================
## -------------------贡献柱状图--------------------------
p1 <- netAnalysis_contribution(cellchat, signaling = pathways.show,
                               title =  pathways.show)                 # 展现对特定通路的贡献程度
p2 <- netAnalysis_contribution(cellchat, signaling = df.net$pathway_name)
cowplot::plot_grid(p1, p2, align = "h",ncol=2)

## ---------------可视化由单个配体-受体对介导的细胞间通讯--------------------
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)
pairLR
LR.show <- pairLR[2,]    # show one ligand-receptor pair
### Hierarchy plot
vertex.receiver = seq(1,6) # a numeric vector
p1 <- netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver,
                     layout = "hierarchy")
# Circle plot
p2 <- netVisual_individual(cellchat, signaling = pathways.show, 
                     pairLR.use = LR.show, layout = "circle")  # "chrod"
cowplot::plot_grid(p1, p2 ,align = "h",ncol=2)

## ---------------通路展示-------------------
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
dir.create('04.pathwat.show')
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  # netVisual(cellchat, signaling = pathways.show.all[i], 
  #           vertex.receiver = vertex.receiver, layout = "hierarchy")#直接出pdf与svg
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0('04.pathwat.show/',pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}
#循环存出连线图
gg

## ---------------Bubble plot-------------------
### sources.use是发出信号的细胞系,target.use是接受信号的细胞系
levels(cellchat@idents) 
netVisual_bubble(cellchat, 
                 sources.use = seq(1:5), 
                 targets.use = c(7:9), 
                 remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)

### 增加signaling参数用于展示特定的配受体
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = seq(1:9), 
                 targets.use = c(7), 
                 signaling = c("CXCL"),
                 remove.isolate = FALSE)
ggsave("bubbleplot2.pdf",width = 5,height = 10)

# 指定信号通路展示配受体
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL","SPP1"))
netVisual_bubble(cellchat, sources.use = c(1:9),
                 targets.use = c(7), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("bubbleplot-LR.pdf",width = 5,height = 10)

# 可以通过增加下面的参数去设置X轴上的顺序
# sort.by.target = T
# sort.by.source = T
# sort.by.source = T, sort.by.target = T
# sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE

## -------使用小提琴/点图绘制信号转导基因表达分布-----------
# CellChat可以使用Seurat包装函数plotGeneExpression绘制与L-R对或信号通路相关的信号转导基因的基因表达分布。
# 该功能提供 “violin”、“dot”、“bar” 三种类型的可视化。
# 或用户可以使用 extractEnrichedLR 提取与推断的 L-R 对或信号通路相关的信号转导基因，然后使用Seurat或其他软件包绘制基因表达。
plotGeneExpression(cellchat, signaling = "VEGF", 
                   enriched.only = TRUE, 
                   type = "violin")

## --------计算并可视化网络中心性得分--------------------
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)


## --------在二维空间中可视化占优势的发送者(源)和接收者(目标)---------
### 从所有信号通路对聚合细胞-细胞通信网络的信号作用分析
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
### 对特定细胞间通讯网络的信号作用分析
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL"))
gg2
gg1 + gg2

## ---------识别对某些细胞群的输出或输入信号贡献最大的信号------------
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
# ht1
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
# ht2
# ht1 + ht2
# class(ht1)

# 特定的signaling
cellchat@netP$pathways
htout <- netAnalysis_signalingRole_heatmap(cellchat, 
                                        pattern = "outgoing",
                                        signaling = c("ICAM","TGFb"))
htout

htcome <- netAnalysis_signalingRole_heatmap(cellchat, 
                                        pattern = "incoming",
                                        signaling = c("ICAM","TGFb"))
htcome


## ----------识别整体通信模式/以探索多种细胞类型和信号通路如何协调运作-----------
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

# 当输出模式的数量为5时，Cophenetic值和Silhouette值都开始突然下降。
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", 
                                          k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "incoming")
# 当输出模式的数量为6时，Cophenetic值和Silhouette值都开始突然下降。
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", 
                                          k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")
