# ============================================ Prepare Environment ===================================================
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
library(clustree)
library(cowplot)
library(stringr)
library(ggsci)
library(SCP)
library(pheatmap)
library(ggrepel)
library(msigdbr)
library(GSVA)
set.seed(1234)
list.files()
dir.create("../03.Output/")


# ============================================ Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")


Idents(LMS) <- "celltype"
TC <- subset(LMS, idents = c("Tumor cell"))
DimPlot(TC)
##计算平均表达量
TC.average <- AverageExpression(TC)
TC.average$RNA
mymatrix  <- as.data.frame(TC@assays$RNA$data["KMT2C",])
mymatrix$celltype <- TC$celltype
colnames(mymatrix)[1]<- "KMT2C expression"
library(ggpubr)
mymatrix <- mymatrix %>% filter(mymatrix$`KMT2C expression` > 0)
p <- ggboxplot(mymatrix, x = "celltype", y = "KMT2C expression", 
               color = "celltype",
               palette = "jco",  # 可选的调色板，也可自定义
               add = "jitter",  # 添加散点，展示数据分布
               xlab = "Cell Type",
               ylab = "KMT2C Expression",
               title = "KMT2C Expression in metastatic Tumorcell")+
  stat_compare_means(
    method = "t.test",  # 两组比较用t检验；多组可用"anova"或"kruskal.test"
    label = "p.signif", # 显示格式化的P值（如p=0.032）
    label.x = 1.5,      # P值标签的x位置（可根据需要调整）
    color = "black"     # P值标签颜色
  )
p
p1 <- ggviolin(mymatrix, x="celltype", y="KMT2C Expression",fill = "celltype",
               palette = c("#00AFBB", "#FC4E07"), 
               add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(method = "t.test",  # 两组比较用t检验；多组可用"anova"或"kruskal.test"
                     label = "p.signif", #
                     label.x = 1.5,      # P值标签的x位置（可根据需要调整）
                     color = "black") +  
  stat_compare_means(label.y = 50) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1))  # breaks 设置刻度间隔
p1
