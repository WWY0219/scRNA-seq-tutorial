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
library(patchwork)
library(RColorBrewer)
library(scRNAtoolVis)
library(CytoTRACE2)
library(monocle)
library(ggpubr)
library(argparse)
library(tidyverse)
set.seed(1234)
# 查看工作路径下的文件
list.files()
outdir <- getwd()
outputDir <- file.path(outdir,"monocle2")
if (!dir.exists(paste0(outputDir))){
  dir.create(paste0(outputDir))
}

# ============================================================ Load scData for TraceAnalysis ==================================================
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

# ============================================================= Color =============================================================================
MYCOLOR<-c("#FF6700", "#9370DB", "#F8D568","#00AD43", "#89CFF0", "#BA160C",
           "#FF91AF", "#A6A6A6", "#006DB0","#C154C1", "#D99A6C", "#96C8A2",
           "#FBEC5D", "#77B5FE", "#E29CD2","#007F5C", "#ACBF60", "#7B68EE",
           "#00FFCD", "#D3AF37", "#50C878","#1974D2", "#FF0000", "#9F00FF",
           "#91A3B0", "#8B4513", "#4166F5","#C19A6B", "#6EAEA1", "#39FF14",
           "#0247FE", "#AF6E4D", "#FF66CC","#9400D3", "#A57164", "#00FFFF",
           "#DF00FF", "#FFAA1D", "#20B2AA","#CB6D51", "#5218FA", "#BC8F8F","#000000")
MYCOLOR2=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF",
            "#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF",
            "#d1cc78","#b2387e","#E64B354C","#4DBBD54C",
            "#00A0874C","#3C54884C","#F39B7F4C","#8491B44C",
            "#91D1C24C","#DC00004C","#fcec7c","#ff92e0",
            "#E64B357F","#4DBBD57F","#00A0877F","#3C54887F",
            "#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F",
  "#eae8c1","#ef4db1","#E64B35B2","#4DBBD5B2",
  "#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2",
  "#91D1C2B2","#DC0000B2","#edd600","#c64ea4" )
#其他配色
cols <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3', 'darkmagenta', 'black', 'grey','lightpink', 'red', 'orange', 'gold', 'darkgreen')

#state 配色
colsp <-c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF','#FD7446FF','#D5E4A2FF','#197EC0FF','#F05C3BFF','#46732EFF',
          '#71D0F5FF','#370335FF','#075149FF','#C80813FF','#91331FFF','#1A9993FF','#FD8CC1FF','#FF6700','#9370DB',
          '#F8D568','#00AD43','#89CFF0','#BA160C','#FF91AF','#A6A6A6','#006DB0','#C154C1','#D99A6C','#96C8A2','#FBEC5D')


# ============================================================= Prepare InputData =============================================================================
## 提取表达数据
count_matrix <- GetAssayData(smc, assay = "RNA", layer = "counts")
count_count<- as(as.matrix(count_matrix),'spareMatrix')
## 提取表型数据
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
## 提取基因信息
fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
fd <- new('AnnotatedDataFrame', data = fData)

## ------------------------------------------构建CDS文件--------------------------------------------------
### 将FPKM/TPM数据转换为UMI数据（read count）rpc_matrix <- relative2abs(cd)
cds <- newCellDataSet(count_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()       # expressionFamily: 数据为TPM/FPKM时设置为tobit(Lower = 0.1)，数据为count时设置为negbinomial.size())
                     )

## ------------------------------------------估计size factor和离散度--------------------------------------
### size factors有助于消除细胞间mRNA捕获的差异；dispersions用于后续的差异表达分析
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# ============================================================= 轨迹构建 =============================================================================
## ------------------------------Selecting genes with high dispersion across cells---------------
## ---------挑选细胞间高度离散的基因 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值------------

### Method-1(使用seurat选择的高变基因）
express_genes <- VariableFeatures(seurat_obj)
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

### Method-2(使用clusters差异表达基因）
deg.cluster <- FindAllMarkers(seurat_obj)
express_genes <- subset(deg.cluster, p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

### Method-3（使用monocle选择的高变基因）
disp_table <- dispersionTable(cds)
ordering_genes_temp <- subset(disp_table, mean_expression >= 0.1  & dispersion_empirical >= 1 * dispersion_fit) 
ordering_genes<-ordering_genes_temp$gene_id
cds <- setOrderingFilter(gbm_cds, ordering_genes)
write.table(ordering_genes, file = paste0(outputDir,"/","dispersion.ordering_genes.xls"), row.names = T, quote = F,sep='\t',col.names=T)
pdf(file = paste(outputDir,"ordering_genes.pdf",sep='/'), width = 9, height = 5)
print(plot_ordering_genes(cds))
dev.off()



## ----------------------------------------------------逆时间轴轨迹构建和在逆时间内排列细胞-----------------------------------------------
### 用DDRTree构建trajectory
tryCatch({
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
}, error = function(e) {
  print("Error in reduceDimension(). Try to apply auto_param_selection = F")
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree',auto_param_selection = F)
})

### order cells along the trajectory
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 1)  # 自定义root_state

cds$barcode <- colnames(cds)
rep <- as.data.frame(pData(cds))
write.table(pData(cds), "cell_Pseudotime.txt", row.names = T, quote = F,sep=',')

### 对state_levels画图的大小进行定义
state_levels <- levels(pData(cds)$State)
if (length(state_levels) <= 2){
  widths_state = 8
  heights_state = 5
  nrows_state = 1
}else if(1 < length(state_levels) & length(state_levels) <= 8){
  widths_state = length(state_levels)*1.5
  heights_state = length(state_levels)*1.5
  nrows_state = 2
}else{
  widths_state = length(state_levels)*1.5
  heights_state = length(state_levels)*1.5
  nrows_state = 3
}


## --------------------------------------------------------------Visulazition----------------------------------------------------------------
### 对"celltype" 画图的大小进行定义!!!需要修改！！！
celltype_levels <-levels(cds$celltype)
if (length(celltype_levels) <= 2){
  widths_celltype = 8
  heights_celltype = 5
  nrows_celltype = 1
}else if(1 < length(celltype_levels) & length(celltype_levels) <= 8){
  widths_celltype = length(celltype_levels)*1.5
  heights_celltype = length(celltype_levels)*1.5
  nrows_celltype = 2
}else{
  widths = length(celltype_levels)*1.5
  heights = length(celltype_levels)*1.5
  nrows = 3
}

### 对"sample" 画图的大小进行定义
sample_levels <-levels(cds$orig.ident)
if (length(celltype_levels) <= 2){
  widths_orig = 8
  heights_orig = 5
  nrows_orig = 1
}else if(1 < length(celltype_levels) & length(celltype_levels) <= 8){
  widths_orig = length(celltype_levels)*1.5
  heights_orig = length(celltype_levels)*1.5
  nrows_orig = 2
}else{
  widths_orig = length(celltype_levels)*1.5
  heights_orig = length(celltype_levels)*1.5
  nrows_orig = 3
}

### Pseudotime_Score
p1 <-plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
p1
ggsave(p1,file = paste0(outputDir,"/","1_Pseudotime_trajectory.pdf"), width = 8, height = 6, dpi=300, limitsize = FALSE)
p1 <-plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")+facet_wrap(~orig.ident)
p1
ggsave(p1,file = paste0(outputDir,"/","2_Pseudotime_trajectory_splited_by_orig.pdf"), width = widths_orig, height = heights_orig,limitsize = FALSE)

### state_Score
p2 <-plot_cell_trajectory(cds, color_by = "State") +scale_color_manual(values=colsp)+
  ggtitle("State") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p2
p2 <-plot_cell_trajectory(cds, color_by = "State", cell_size = 0.5,cell_link_size = 0.5) +
  facet_wrap(~State, nrow = nrows_state, scales = "free")+
  scale_color_manual(values=colsp)+ggtitle("State") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p2
p2 <-plot_cell_trajectory(cds, color_by = "State", cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~orig.ident, nrow = nrows_orig,scales = "free")+scale_color_manual(values=colsp)+ggtitle("State") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p2
p2 <- ggplot(rep, aes(Pseudotime, fill = State)) +geom_density() +facet_wrap(~State,nrow =nrows_state, scales = "free_y") +
  theme_bw() +RotatedAxis() +
  theme(
    strip.text = element_blank(),
    strip.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank()
  ) +scale_fill_manual(values = colsp)
p2

### celltype
p3 <-plot_cell_trajectory(gbm_cds, color_by = 'celltype')+scale_color_manual(values=MYCOLOR)+
  ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p3
ggsave(p3 ,file =paste0(outputDir,"/","7_Celltype_trajectory.pdf"), width = 8, height = 6,limitsize = FALSE)

# celltype splited by orig.ident
p3 <-plot_cell_trajectory(gbm_cds, color_by = 'celltype')+scale_color_manual(values=MYCOLOR)+facet_wrap(~orig.ident)+
  ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p3
ggsave(p3 ,file =paste0(outputDir,"/","8_Celltype_trajectory_splited_by_orig.pdf"), width = widths_orig, height = heights_orig,limitsize = FALSE)

p3 <-plot_cell_trajectory(gbm_cds, color_by = 'celltype', cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~celltype, nrow = nrows_celltype, scales = "free")+scale_color_manual(values=MYCOLOR)+ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p3
ggsave(p3 ,file =paste0(outputDir,"/","9_Celltype_trajectory_splited_by_celltype.pdf"), width = widths_celltype, height = heights_celltype,limitsize = FALSE)

p3 <- ggplot(rep, aes(Pseudotime, fill = celltype)) +geom_density() +facet_grid(~celltype, scales = "free_y") +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()
) +scale_fill_manual(values = MYCOLOR)
p3
ggsave(p3 ,file =paste0(outputDir,"/","10_Density_trajectory_splited_by_celltype.pdf"), height = heights_celltype,  width = widths_celltype,limitsize = FALSE)

p3 <- ggplot(rep, aes(Pseudotime, fill = celltype)) +geom_density() +facet_wrap(~orig.ident) +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()) +scale_fill_manual(values = MYCOLOR)
p3
ggsave(p3 ,file =paste0(outputDir,"/","11_Density_trajectory_splited_by_orig.pdf"), height = heights_orig,  width = widths_orig,limitsize = FALSE)

### Sample
p4 <-plot_cell_trajectory(gbm_cds, color_by = 'orig.ident') +scale_color_manual(values=MYCOLOR)+ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p4
ggsave(p4 ,file =paste0(outputDir,"/","12_Sample_trajectory.pdf"), width = 8, height = 6,limitsize = FALSE)

p4 <-plot_cell_trajectory(gbm_cds, color_by = 'orig.ident', cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~celltype, nrow = nrows_celltype, scales = "free")+scale_color_manual(values=MYCOLOR) + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p4
ggsave(p4 ,file =paste0(outputDir,"/","12_Sample_trajectory_splited_by_orig.pdf"), width = widths_celltype, height = heights_celltype,limitsize = FALSE)

p4 <- ggplot(rep, aes(Pseudotime, fill = orig.ident)) +geom_density() +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()) +scale_fill_manual(values = MYCOLOR2)
p4
ggsave(p4 ,file =paste0(outputDir,"/","13_Density_trajectory_splited_by_orig.pdf"), height = 8, width = 6,limitsize = FALSE)

p4 <- ggplot(rep, aes(Pseudotime, fill = orig.ident)) +geom_density() +facet_wrap(~celltype) +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()) +scale_fill_manual(values = MYCOLOR2)
p4
ggsave(p4 ,file =paste0(outputDir,"/","14_Density_trajectory_splited_by_orig.pdf"), height = heights_celltype, width = widths_celltype,limitsize = FALSE)

### 
## ----------------------------------------------------逆时间轴轨迹构建和在逆时间内排列细胞-----------------------------------------------



































