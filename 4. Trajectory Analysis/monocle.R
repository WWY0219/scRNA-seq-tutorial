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

## Method-1(使用seurat选择的高变基因）
express_genes <- VariableFeatures(seurat_obj)
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

## Method-2(使用clusters差异表达基因）
deg.cluster <- FindAllMarkers(seurat_obj)
express_genes <- subset(deg.cluster, p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

## Method-3（使用monocle选择的高变基因）
disp_table <- dispersionTable(cds)
ordering_genes_temp <- subset(disp_table, mean_expression >= 0.1  & dispersion_empirical >= 1 * dispersion_fit) 
ordering_genes<-ordering_genes_temp$gene_id
cds <- setOrderingFilter(gbm_cds, ordering_genes)
write.table(ordering_genes, file = paste0(outputDir,"/","dispersion.ordering_genes.xls"), row.names = T, quote = F,sep='\t',col.names=T)
pdf(file = paste(outputDir,"ordering_genes.pdf",sep='/'), width = 9, height = 5)
print(plot_ordering_genes(cds))
dev.off()
















