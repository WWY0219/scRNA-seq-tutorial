setwd("G:/#Biotechnologic Analysis###/GSE237638（ULMS-scRNA）")
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(harmony)
makecore <- function(workcore,memory){
  if(!require(Seurat))install.packages('Seurat')
  if(!require(future))install.packages('future')
  plan("multisession", workers = workcore)
  options(future.globals.maxSize= memory*1024*1024**2)
}
makecore(8,8)                              #这里以八线程,8GB为例

colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

LMS <- readRDS("./01.Data/LMS_seurat_res_0.5.rds")
set.seed(1234)
Idents(LMS) <- "clusters_res0.5"
## LMS[["RNA"]] <- JoinLayers(LMS[["RNA"]]) ##h合并layers（数据整理）

cell_markers <- list(
  
  MSC =c("CD34","PDGFRA"),
  #免疫细胞
  immunocell= c("PTPRC"),
  T_cell = c("CD3E", "CD3D", "CD4", "CD8A"),
  B_cell = c("MS4A1", "CD19", "CD79A"),
  Monocyte_Macrophage = c("CD14", "CD68", "CSF1R", "FCGR3A", "LYZ","ITGAM"),
  Dendritic_cell = c("CD1C", "CLEC9A", "HLA-DRA", "CD83"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "NCAM1"),  # NCAM1即CD56
  Mast_cell = c("CPA3", "TPSAB1", "KIT"),
  # 基质细胞
  Fibroblast = c("COL1A1", "COL3A1", "THY1", "DCN", "FAP"),  # THY1即CD90
  Endothelial_cell = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),  # PECAM1即CD31
  Smooth_muscle_cell = c("ACTA2", "MYH11", "TAGLN"),
  SMC_cancercell =c("ACTG2","PRLR","SFRP4"),
  ESC =c("PLN","RGS5","SUSD2"),
  #增值细胞 
  Proliferating_cell = c("MKI67", "PCNA", "TOP2A", "CCNB1"),
  
  # 其他常见细胞
  Epithelial_cell = c("EPCAM", "KRT8", "KRT18", "CDH1")
)
p <- DotPlot(
  object = LMS,
  features = cell_markers,  # 修正：使用正确的marker列表名称
  cols = c("grey", "red"),
  cluster.idents = TRUE
) +
  RotatedAxis() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),  # 修正：添加fill=NA避免覆盖
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin = margin(b = 3, unit = "mm")),
    strip.placement = 'outlet',
    axis.line = element_blank()
  ) +
  labs(x = "", y = "")
p

new.cluster.ids <- c("0"= "Tumor cell",
                     "1"="Endothelial cell",
                     "2"= "CAF",
                     "3"= "Monocyte",
                     "4"="Tumor cell",
                     "5"="Metastatic Tumor cell",
                     "6"="Metastatic Tumor cell",
                     "7"= "Metastatic Tumor cell", 
                     "8"= "T/NK cell",
                     "9"="Vascular Smooth Muscle Cell",
                     "10"="Monocyte",
                     "11"="Metastatic Tumor cell",
                     "12"= "Tumor cell",
                     "13"= "CAF",
                     "14"="Neutrophils",
                     "15"="Mast cell",
                     "16"="Endothelial cell",
                     "17"= "db-like",
                     "18"= "CAF",
                     "19"="Metastatic Tumor cell",
                     "20"="Lymphatic Endothelial Cell",
                     "21"="db-like",
                     "22"= "Metastatic Tumor cell",
                     "23"= "db-like")

FeaturePlot(LMS, reduction = "umap_harmony_res0.5", features = c("THY1","PDPN"))
FeaturePlot(LMS, reduction = "umap_harmony_res0.5", features = c("PECAM1","CDH5"))

LMS <- RenameIdents(LMS, new.cluster.ids)
LMS@meta.data$celltype <- Idents(LMS)

LMS <- subset(LMS, idents = c("db-like"), invert = T)
unique(Idents(LMS))
DimPlot(LMS, reduction = "umap_harmony_res0.5", label = TRUE, pt.size = 1.5)




####decontx####
library(decontX)
counts <- GetAssayData(LMS,layer = "data",assay = 'RNA')
decontX_results <- decontX(counts)
LMS$Contamination =decontX_results$contamination
FeaturePlot(LMS, 
            features = 'Contamination', 
            raster=FALSE       # 细胞过多时候需要加这个参数
) + 
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('scVI_UMAP_1')+
  ylab('scVI_UMAP_2')
LMS = LMS[,LMS$Contamination < 0.2]


LMS.markers <- FindAllMarkers(LMS, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
top5 <- LMS.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(LMS, features = top5$gene,angle = 90,size=3) +
  # 更换颜色（例如从蓝色到红色的渐变，可自定义）
  scale_fill_gradientn(colors = c("#3182BD", "white", "#E6550D"))  # 蓝-白-红
ggsave("./doheatmap.pdf",width = 10,height = 11,dpi = 300)
ggsave("./doheatmap.png",width = 10,height = 11,dpi = 300)
# 假设细胞类型注释存储在 meta.data$celltype 中
# 统计每个 orig.ident × celltype 的细胞数量
LMS$celltype <- droplevels(LMS$celltype)
Idents(LMS) <- "celltype"
celltype_counts <- table(LMS$orig.ident, LMS$celltype) %>% 
  as.data.frame() %>% 
  rename(orig.ident = Var1, celltype = Var2, count = Freq)

# 计算每个 orig.ident 中各细胞类型的占比（百分比）
celltype_proportions <- celltype_counts %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100)
p_prop <- ggplot(celltype_proportions, aes(x = orig.ident, y = percentage, fill = celltype)) +
  geom_col(position = "fill") +  # 百分比堆积
  scale_y_continuous(labels = scales::percent_format()) +  # y轴显示百分比
  labs(x = "orig.ident", y = "celltype percentage", fill = "celltype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16))

print(p_prop)

#样本分组
LMS@meta.data$group <- ifelse(
  LMS@meta.data$orig.ident %in% c("A32_T", "C10_T"),         # 对照样本ID
  "Primary",                                             # 分组结果：原发
  "Metastatic"                                            # 分组结果：转移
)
DimPlot(LMS, reduction = 'umap_harmony_res0.5', 
        group.by = "group",
        label = F, 
        pt.size =2)

#设置group顺序
LMS@meta.data$group <- factor(LMS@meta.data$group, 
                              levels=c("Primary", "Metastatic") )
#设置celltype的顺序
LMS@meta.data$celltype <- factor(LMS@meta.data$celltype, levels=c("Tumor cell", "Metastatic Tumor cell",  "CAF","Endothelial cell", "Monocyte", "T/NK cell","Mast cell", 
                                                                 'Neutrophils',"Lymphatic Endothelial Cell","Vascular Smooth Muscle Cell"))

####scRNAtools####
library(scRNAtoolVis)
clusterCornerAxes(object = LMS,
                  reduction = 'umap_harmony_res0.5',
                  clusterCol = "celltype",
                  noSplit = T,
                  cellLabel = T,
                  cellLabelSize = 5)
# 1. 确保细胞类型标识正确（设置Idents为celltype）
Idents(LMS) <- "celltype"

# 2. 获取所有细胞类型
celltypes <- levels(Idents(LMS))

# 3. 初始化列表存储每个细胞类型的Top5标记基因
top5_markers_list <- list()

# 4. 循环计算每个细胞类型的标记基因，并提取Top5
for (ct in celltypes) {
  # 计算当前细胞类型与其他所有细胞的差异基因（默认按avg_log2FC降序）
  markers <- FindMarkers(
    object = LMS,
    ident.1 = ct,          # 目标细胞类型
    ident.2 = NULL,        # 与所有其他细胞比较
    logfc.threshold = 0.25, # 最小log2倍数变化（可调整）
    min.pct = 0.1,         # 最小表达比例（可调整）
    verbose = FALSE
  )
  
  # 为结果添加基因名列（从行名转换），并提取Top5
  markers <- markers %>%
    rownames_to_column("gene") %>%  # 行名转为gene列
    arrange(desc(avg_log2FC)) %>%   # 按avg_log2FC降序排序
    slice_head(n = 5) %>%           # 取前5行（Top5）
    column_to_rownames("gene")      # 恢复基因为行名
  
  # 存储到列表
  top5_markers_list[[ct]] <- markers
}

# 5. 查看结果（示例：打印第一个细胞类型的Top5标记基因）
print("第一个细胞类型的Top5标记基因：")
print(top5_markers_list[[1]])

# 6. （可选）将所有结果合并为一个数据框
top5_markers_df <- do.call(rbind, lapply(
  names(top5_markers_list), 
  function(ct) {
    df <- top5_markers_list[[ct]] %>%
      rownames_to_column("gene") %>%
      mutate(celltype = ct)  # 添加细胞类型列
    return(df)
  }
))

markerVolcano(markers = cell_markers,
              topn = 5,
              labelCol = ggsci::pal_npg()(9))







####KMT2C出图####
Idents(LMS) <- "celltype"
TC <- subset(LMS, idents = c("Metastatic Tumor cell","Tumor cell"))
DimPlot(TC)
##计算平均表达量
TC.average <- AverageExpression(TC)
TC.average$RNA


mydeg <- FindMarkers(TC,ident.1 = 'Metastatic Tumor cell',ident.2 = 'Tumor cell', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
mydeg$SYMBOL <- rownames(mydeg)
mygene <- mydeg %>% filter(p_val_adj != 0) %>% top_n(n = -100, wt = p_val_adj) %>% rownames()
mygene <- mydeg %>% filter(p_val_adj != 0) %>% filter(p_val_adj >-  0.05) %>% rownames()
mygene[1:10]

library(org.Hs.eg.db)  ##org.Hs.eg.db-用于人类的基因注释包
library(clusterProfiler)
gene.df <- mapIds(org.Hs.eg.db, keys = mygene, keytype = "SYMBOL", column = "ENTREZID")
go_enrichment <- enrichGO(gene = gene.df,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "All",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
go_enrichment1  <- DOSE::setReadable(go_enrichment, OrgDb='org.Hs.eg.db',keyType='ENTREZID') 
barplot(go_enrichment1, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
dotplot(go_enrichment1, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
KEGG_enrichment  <- enrichKEGG(gene = gene.df, organism = 'hsa', pvalueCutoff = 0.05)
barplot(KEGG_enrichment, showCategory = 30, title = "KEGG Enrichment Analysis") 
eKEGG <- setReadable(KEGG_enrichment, OrgDb = org.Hs.eg.db, keyType="ENTREZID")       ##将geneID转换为gene symbol
browseKEGG(KEGG_enrichment, 'hsa04820')

write.csv(eKEGG@result,"KEGG-TJUS-R.csv")
df <- read_csv("KEGG-TJUS-R.csv") %>% 
  select(2,5,9,11,12) %>%
  filter(category %in% c("Environmental Information Processing","Human Diseases","Organismal Systems")) %>% 
  mutate(`-log10(p.adjust)`=-log10(p.adjust)) %>% 
  arrange(category,desc(`-log10(p.adjust)`))

dff <- df %>%  mutate(col=case_when(category=="Human Diseases" ~ "#7294D4",
                                    category=="Environmental Information Processing" ~ "#00A08A",
                                    category=="Organismal Systems" ~ "#F98400",
                                    TRUE ~ "white"))

col_mapping <- dff %>% select(Description,col) %>% deframe()

dff$Description <- factor(dff$Description,levels = dff$Description %>% rev())

p_KEGG <- ggplot(dff,aes(`-log10(p.adjust)`,Description))+
  geom_bar(stat = "identity",aes(fill=category),show.legend = F,width=0.5)+
  geom_text(data=dff,
            aes(label=Count),vjust=0.5,hjust=0,size=3)+
  scale_x_continuous(expand= expansion(mult = c(0,0.05)))+
  scale_y_discrete(expand= expansion(mult = c(0.03,0)))+
  scale_fill_manual(values = c("#00A08A","#7294D4","#F98400"),
                    na.translate = FALSE)+
  scale_color_manual(values = c("#7294D4","#00A08A","#F98400"))+
  labs(y=NULL)+
  theme(axis.text.y=element_text(color=rev(col_mapping)),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        axis.line.x=element_line(color="black"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),unit="cm"))+
  annotate(geom="rect",xmin=-0.15,xmax=-0.01,ymin=0.2,ymax=4.5,
           fill="#F98400",alpha=0.5)+
  annotate(geom="rect",xmin=-0.15,xmax=-0.01,ymin=4.5,ymax=31.5,
           fill="#7294D4",alpha=0.5)+
  annotate(geom="rect",xmin=-0.15,xmax=-0.01,ymin=31.5,ymax=38.3,
           fill="#00A08A",alpha=0.5)+
  annotate(geom="text",x=-0.1,y=32.5,label="Environmental Information Processing",size=3,color="black",
           angle=90,fontface="bold")+
  annotate(geom="text",x=-0.1,y=18,label="Human Diseases",size=3,color="black",
           angle=90,fontface="bold")+
  annotate(geom="text",x=-0.1,y=3,label="Organismal Systems",size=3,color="black",
           angle=90,fontface="bold")
p_KEGG
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














#寻找cluster的marker
LMS.markers <- FindAllMarkers(LMS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
LMS.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- LMS.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(LMS, features = top10$gene) + NoLegend()                           # Doheatmap图
VlnPlot(LMS, features = top10$gene[1:20], pt.size=0)                         #小提琴图观察基因分布
DimPlot(LMS,reduction = 'umap',label = T)
DimPlot(LMS,group.by = 'orig.ident',label=T)
marker2 <- FindMarkers(object = LMS, ident.1 = 2)
marker_DEG <- FindMarkers(object = LMS, ident.1 = 10, ident.2 = 11,           # 第二个cluster（对比组）
                        logfc.threshold = 0.25,  # 最小log2倍数变化（过滤微小差异）
                        min.pct = 0.1,           # 基因在至少10%的细胞中表达（过滤低表达基因）
                        test.use = "wilcox",     # 统计检验方法（默认wilcox，适用于单细胞数据）
                        assay = "RNA"            # 使用的assay（默认"RNA"）
)
FeaturePlot(LMS, features = c("THY1"))

#细胞群表达点状图#
markers <- c("PDGFRA","CDH5","VWF","PECAM1","CPA3","MS4A1","ASPM",
             "PDPN","CD3E","CD14","HIGD1B",
             "THY1","CSF3R","CD68","CXCR4","PTPRC","S100A8","MKI67","TOP2A")
p <- DotPlot(LMS,features = markers)+coord_flip()
p + theme(
  axis.text.x = element_text(
    angle = 90,                 # 旋转角度，90 度为垂直
    hjust = 1,                  # 水平对齐（1 表示右对齐，避免文字溢出）
    vjust = 0.5                 # 垂直对齐
  )
)

DimPlot(LMS, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
names(LMS@meta.data)
DimPlot(LMS,split.by = 'orig.ident')


#细胞组分柱状图
sample_table <- as.data.frame(table(LMS@meta.data$orig.ident,LMS@meta.data$celltype))
names(sample_table) <- c("Samples","celltype","CellNumber")
plot_sample<-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage") + RotatedAxis()
plot_sample

#分组柱状图
group_table <- as.data.frame(table(LMS@meta.data$group,LMS@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
plot_group<-ggplot(sce2@meta.data,aes(x=group,fill=celltype))+
  geom_bar(position="fill")+
  scale_fill_manual(values=colour) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")
plot_group

plot_sample + plot_group + 
  plot_layout(widths = c(2, 1))

####寻找差异基因####
LMS <- readRDS("./03.Output/LMS_anno.rds")
##自动计算差异基因
cellfordeg<-levels(LMS$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(LMS, ident.1 = paste0(cellfordeg[i],"_lms001_T"), ident.2 = paste0(cellfordeg[i],"_A32_T"), verbose = FALSE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}
list.files()

#自动计算后
top10 <- CELLDEG  %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10
LMS <- ScaleData(LMS, features =  rownames(LMS))                     #对单细胞数据中所有基因进行标准化和中心化 
DoHeatmap(LMS,features = top10,size=3)

#根据celltype绘制小提琴图
Idents(LMS) <- "celltype"
pvlnplot <- VlnPlot(LMS,features = c("KMT2C","KMT2D"),idents = c("Tumor cell","Metastatic Tumor cell"))
findPvalue <- FindMarkers(LMS, 
                         ident.1 = "Tumor cell",         # 第一类细胞（待比较组）
                         ident.2 = "Metastatic Tumor cell",   # 第二类细胞（参照组））
                         features = c("KMT2C", "KMT2D")      # 目标基因
  )
p_values <- data.frame(
  gene = rownames(findPvalue),
  p_val = sprintf("p = %.3e", findPvalue$p_val)                # 科学计数法显示（如p = 2.345e-05）
)
p_val_text <- sprintf("p = %.3e", findPvalue$p_val)
p7 <- pvlnplot +
  annotate(
    "text",
    x = 1.5,                                          # x轴位置（两组中间）
    y = max(ggplot_build(pvlnplot)$data[[1]]$y) * 0.9,  # y轴位置（接近最大值的90%处）
    label = p_val_text,
    size = 4,  # 字体大小
    fontface = "bold"  # 加粗
  ) +
  ggtitle("KMT2C")  
p7

#点状图
DotPlot(LMS, features = c("KMT2C","KMT2D"),split.by ='group',cols = c('blue','yellow','pink')) 



##提取表达量（KMT2C&KMT2D）
genes <- c("KMT2C","KMT2D")
mymatrix  <- as.data.frame(LMS@assays$RNA$data[genes,])
mymatrix2 <-t(mymatrix)%>%as.data.frame()
mymatrix2$celltype<-LMS$celltype
mymatrix2[,ncol(mymatrix2)+1]<-LMS$group
colnames(mymatrix2)[ncol(mymatrix2)] <- "group"

#绘图
library(ggplot2)
p8 <- ggplot2::ggplot(mymatrix2,aes(x=celltype,y=KMT2C,fill=group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange','pink'))
p8


#手动寻找差异基因
diff_genes <- FindMarkers(LMS, ident.1 = "metastatic Tumorcell", ident.2 = "Malignant cell", logfc.threshold = 0.25,
                          test.use = "wilcox",
                          layer = "data",
                          min.pct = 0.1,
                          min.diff.pct = -Inf)
head(diff_genes)
diff <- subset(diff_genes,diff_genes$p_val>0 & diff_genes$p_val < 0.05)
#寻找celltype.group之间的差异基因
Idents(LMS) <- "celltype.group"
mydeg <- FindMarkers(LMS,ident.1 = 'Malignant_lms001',ident.2 = 'SMC_lms001', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)

#提取未完成
Idents(LMS) <- "celltype"
TC <- as.matrix(GetAssayData(LMS,layer = 'data'))[,WhichCells(LMS,idents= c("metastatic Tumorcell","Malignant cell"))]
data <- LayerData(LMS, assay = "RNA", layer = "counts",idents= c("metastatic Tumorcell","Malignant cell"))
data <- as(data, "sparseMatrix")          #稀疏矩阵 
data <- as.data.frame(data)














