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



# ============================================ GSVA Preparation ===================================================
## windwo处理
genesets <- msigdbr(species = "Homo sapiens", category = "C2") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

Idents(seurat_obj) <- seurat_obj$celltype assays = "RNA", slot = "data")[[1]]
vaPar <- gsvaParam(expr, genesets,maxDiff = TRUE)
gsvaPar 
gsva.res <- gsva(gsvaPar)
dim(gsva.res)
head(gsva.res)[1:5,1:5]
write.csv(gsva.res,"gsva.res.csv")

gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
gsva_d = gsva.res[sample(nrow(gsva.res),30),]
pheatmap::pheatmap(gsva_d, show_colnames = T, 
                   scale = "row",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

#数据来源是上面的
# 气泡图
library(reshape2)
gsva_long <- melt(gsva_d, id.vars = "Genesets")

# 创建气泡图
ggplot(gsva_long, aes(x = Var2, y = Var1, size = value, color = value)) +
  geom_point(alpha = 0.7) +  # 使用散点图层绘制气泡，alpha设置点的透明度
  scale_size_continuous(range = c(1, 6)) +  # 设置气泡大小的范围
  theme_bw() + 
  scale_color_gradient(low = "#336699", high =  "tomato") +
  labs(x = "Gene Set", y = "Sample", size = "GSVA Score")+
  ggtitle("GSVA analysis") +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5),
        plot.title = element_text(hjust = 0.5))
