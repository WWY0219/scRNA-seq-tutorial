####gsva####

library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(ggsci) 
library(patchwork) 
library(ggpubr)
library(RColorBrewer) 
library(msigdbr)
library(qs)
library(GSVA)
sce <- smc
genesets <- msigdbr(species = "Homo sapiens", category = "C2") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
#  这里的分组会决定接下来的分析哦！
Idents(sce) <- sce$subcelltype
expr <- AverageExpression(sce, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
head(expr)
# gsva默认开启全部线程计算
gsvaPar <- gsvaParam(expr, genesets,maxDiff = TRUE)
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
