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
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(patchwork) 
library(ggsci)


setwd("E:/#WWY_Project/课题二子宫肌瘤发生发展/Biotech/02.scRNA-seq")
smc <- qread("./01.Data/seurat_obj_SMC_Fibro.qs")
# 相关性分析
av <-AggregateExpression(smc,
                         group.by = c("orig.ident","subcelltype"),
                         assays = "RNA",
                         #layer = "counts",
                         return.seurat = FALSE)  # 返回总的计数 
av=as.data.frame(av[[1]])
head(av)[1:3,1:3] # 可以看到是整数矩阵
# 相关性分析
# 找到sd最显著的排名前1000基因
cg=names(tail(sort(apply(log(av+1), 1, sd)),1000)) 
# 提前这些基因矩阵并进行log处理后计算不同样本对应的细胞的相关性值
df =cor(as.matrix(log(av[cg,]+1)))
colnames(df)
ac=as.data.frame(str_split(colnames(df),'_',simplify = T))
rownames(ac)=colnames(df)
ac$V1
ac$group = ifelse(grepl('GSM4942397|GSM4942396|GSM4942398|GSM4942399|GSM5023319',ac$V1),"MM","LM")
table(ac$group)
head(ac)
pheatmap::pheatmap(df ,
                   show_colnames = F,
                   show_rownames = F,
                   annotation_col = ac ) 
dev.off()
# 3.对于每个细胞亚群进行主成分分析##################
av <-AggregateExpression(smc ,
                         group.by = c("orig.ident","subcelltype"),
                         assays = "RNA") 
av=as.data.frame(av[[1]])
df=log(av +1) 
head(ac)
celltp = unique(ac$V2);celltp
library(patchwork) 
library("FactoMineR") 
library("factoextra")  
library(ggstatsplot)
library(pheatmap)
pca_list <- lapply(celltp, function(x){
  # x=celltp[1]
  x
  #～～～主成分分析图p2～～～
  exp <- df[,rownames(ac[ac$V2==x,])]  
  cg=names(tail(sort(apply(exp, 1, sd)),1000)) 
  exp=exp[cg,]
  dat.pca <- PCA(as.data.frame(t(exp)) , graph = FALSE)
  group_list=ac[ac$V2==x,'group']
  this_title <- paste0(x,'_PCA')
  p2 <- fviz_pca_ind(dat.pca,
                     geom.ind = "point", # show points only (nbut not "text")
                     col.ind = group_list, # color by groups
                     palette = "Dark2",
                     addEllipses = TRUE, # Concentration ellipses
                     legend.title = "Groups")+
    ggtitle(this_title)+
    theme_ggstatsplot()+
    theme(plot.title = element_text(size=12,hjust = 0.5))
  
  p2
})
wrap_plots(pca_list, byrow = T, nrow = 2 )




#####GSEA####
av <-AggregateExpression(smc,
                         group.by = c("orig.ident","group"),
                         assays = "RNA",
                         slot = "counts",
                         return.seurat = FALSE)  # 返回总的计数 
av=as.data.frame(av[[1]])
head(av)[1:3,1:3] # 可以看到是整数矩阵
library(tibble)
library(DESeq2)
# 1. Get counts matrix
counts_res <- av
head(counts_res)

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_res))
colData <- colData %>%
  mutate(condition = ifelse(grepl('Normal', samples), 'Normal', 'Leiomyoma')) %>%
  column_to_rownames(var = 'samples')

# 3. perform DESeq2 --------
# Create DESeq2 object  
dds <- DESeqDataSetFromMatrix(countData = counts_res,
                              colData = colData,
                              design = ~ condition) # condition 表示差异分析将基于colData的condition 变量进行

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
# run DESeq2
dds <- DESeq(dds)
# Check the coefficients for the comparison
resultsNames(dds)
# Generate results object
res <- results(dds, name = "condition_Normal_vs_Leiomyoma")
res
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG_deseq2 = na.omit(DEG)

#添加上下调信息
DEG_deseq2 <- DEG_deseq2 %>%
  mutate(Type = if_else(padj > 0.05, "stable",
                        if_else(abs(log2FoldChange) < 1, "stable",
                                if_else(log2FoldChange >= 1, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Symbol")
# ggplot绘图
ggplot(DEG_deseq2, aes(log2FoldChange,-log10(padj))) +
  geom_point(size = 3.5, alpha = 0.8,
             aes(color = Type),show.legend = T)  +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'black',lwd=0.8) + 
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'black',lwd=0.8)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
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
