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
library(harmony)
library(patchwork)
library(RColorBrewer)
library(scRNAtoolVis)
library(CytoTRACE2)
library(monocle3)
library(ggpubr)
set.seed(1234)
# 查看工作路径下的文件
list.files()

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

# ============================================================ Run CytoTRACE2 ===================================================================
seurat_obj_cytotrace <- cytotrace2(seurat_obj,
                        is_seurat = TRUE,
                        slot_type = "counts",
                        species = "human",
                        seed = 1234)
## replace UMAP/TSNE
seurat_obj_cytotrace@reductions[["umap"]] <- seurat_obj@reductions[["tsne"]] 

# ============================================================ Plot CytoTRACE2 =================================================================
## making an annotation dataframe that matches input requirements for plotData function
annotation <- data.frame( phenotype=seurat_obj_cytotrace@meta.data$subcelltype) %>% set_rownames(., colnames(seurat_obj_cytotrace))
## cytoTRACE2 plot
plots <- plotData(cytotrace2_result = seurat_obj_cytotrace, 
                  annotation = annotation, 
                  is_seurat = TRUE,
                  seed = 1234)
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno
cyto_fig  <- (p1+p2+p3+p4) + plot_layout(ncol = 2)
cyto_fig
ggsave("../03.Output/cytoTRACE2_primaryfig.pdf", plot = cyto_fig, width = 10, height = 10, dpi =300)


# ============================================================ plotData 非UMAP降维版本 =================================================================
seurat <- seurat_obj_cytotrace                                                                           # 你的Seurat对象（含CytoTRACE2结果和降维）
reduction_name <- "umap"                                                                                 # 降维类型（umap/tsne）
labels <- c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent", "Totipotent")
colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")                            # 潜能颜色
rel_colors <- c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")            # 相对顺序颜色

# 提取降维坐标的范围（用于固定坐标轴，避免图形偏移）
embeddings <- seurat@reductions[[reduction_name]]@cell.embeddings
x_limits <- range(embeddings[, 1], na.rm = TRUE)
y_limits <- range(embeddings[, 2], na.rm = TRUE)

# 预处理：计算CytoTRACE2_Score_clipped（用于第一张图的颜色映射）
seurat@meta.data[["CytoTRACE2_Score_clipped"]] <- 5.5 - 6 * seurat@meta.data[["CytoTRACE2_Score"]]
seurat@meta.data[["CytoTRACE2_Score_clipped"]] <- -pmax(pmin(seurat@meta.data[["CytoTRACE2_Score_clipped"]], 5), 0)

## -------------------------------------------------------Potency score-------------------------------------------------------------
p1 <- FeaturePlot(
  seurat, 
  features = "CytoTRACE2_Score_clipped",  # 绘图的特征列
  reduction = reduction_name              # 使用的降维类型
) +
  # 颜色映射：连续渐变，对应潜能标签
  scale_colour_gradientn(
    colours = rev(colors),                # 颜色反转（与潜能匹配）
    na.value = "transparent",             # 缺失值透明
    labels = labels,                      # 图例标签
    limits = c(-5, 0),                    # 颜色范围（与clipped值匹配）
    name = "Potency score \n"             # 图例标题
  ) +
  # 图例样式：黑色边框和刻度
  guides(colour = guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black"
  )) +
  # 坐标轴标签（动态匹配降维类型）
  xlab(paste0(reduction_name, "1")) +
  ylab(paste0(reduction_name, "2")) +
  # 标题和主题样式
  ggtitle("CytoTRACE 2") +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1  # 固定宽高比为1:1
  ) +
  # 固定坐标轴范围（避免不同图的坐标偏移）
  coord_cartesian(xlim = x_limits, ylim = y_limits)

# 显示图1
print(p1)
ggsave("p1.pdf",plot =p1,width=6,heigh=6,dpi=300)


## -------------------------------------------------------Potency category-------------------------------------------------------------
p2 <- DimPlot(
  seurat, 
  reduction = reduction_name,             # 使用的降维类型
  group.by = "CytoTRACE2_Potency",        # 分组依据：潜能分类
  label = FALSE                           # 不显示分组标签
) +
  # 颜色映射：手动指定潜能分类的颜色
  scale_color_manual(
    values = colors,
    name = "Potency category",            # 图例标题
    breaks = rev(labels)                  # 图例顺序反转（与颜色匹配）
  ) +
  # 坐标轴标签
  xlab(paste0(reduction_name, "1")) +
  ylab(paste0(reduction_name, "2")) +
  # 标题和主题样式
  ggtitle("CytoTRACE 2") +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1
  ) +
  # 固定坐标轴范围
  coord_cartesian(xlim = x_limits, ylim = y_limits)

# 显示图2
print(p2)
ggsave("p2.pdf",plot =p2,width=6,heigh=6,dpi=300)


## -------------------------------------------------------CytoTRACE2_Relative_UMAP-------------------------------------------------------------
p3 <- FeaturePlot(
  seurat, 
  features = "CytoTRACE2_Relative",       # 绘图的特征列
  reduction = reduction_name              # 使用的降维类型
) +
  # 颜色映射：连续渐变，对应相对分化顺序
  scale_colour_gradientn(
    colours = rel_colors,                 # 相对顺序的颜色系
    na.value = "transparent",
    limits = c(0, 1),                     # 颜色范围（0~1）
    breaks = seq(0, 1, by = 0.2),         # 颜色刻度
    labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),  # 刻度标签
    name = "Relative\norder \n"           # 图例标题（换行）
  ) +
  # 图例样式
  guides(colour = guide_colorbar(
    frame.colour = "black",
    ticks.colour = "black"
  )) +
  # 坐标轴标签和标题
  xlab(paste0(reduction_name, "1")) +
  ylab(paste0(reduction_name, "2")) +
  ggtitle("CytoTRACE 2") +
  # 主题样式
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    aspect.ratio = 1
  ) +
  # 固定坐标轴范围
  coord_cartesian(xlim = x_limits, ylim = y_limits)

# 显示图3
print(p3)
ggsave("p3.pdf",plot =p3,width=6,heigh=6,dpi=300)


## -------------------------------------------------------CytoTRACE2_Boxplot_byPheno-------------------------------------------------------------
# 步骤1：预处理数据（按表型计算中位数，用于排序和颜色映射）
mtd <- seurat@meta.data[, c("Phenotype", "CytoTRACE2_Score")]  # 提取表型和评分
# 按表型计算中位数并排序
medians <- mtd %>%
  group_by(Phenotype) %>%
  summarise(CytoTRACE2_median_per_pheno = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
  arrange(desc(CytoTRACE2_median_per_pheno))
# 合并中位数到原数据，用于排序表型
mtd <- mtd %>%
  inner_join(medians, by = "Phenotype")
# 将表型转为因子（按中位数排序，确保箱线图顺序）
mtd$Phenotype <- factor(mtd$Phenotype, levels = medians$Phenotype)

# 步骤2：绘制箱线图
p4 <- ggplot(
  mtd[!is.na(mtd$Phenotype), ],  # 过滤缺失表型的行
  aes(x = Phenotype, y = CytoTRACE2_Score)  # x轴：表型，y轴：CytoTRACE2评分
) +
  # 箱线图：按中位数填充颜色，隐藏离群点（用jitter替代）
  geom_boxplot(
    aes(fill = CytoTRACE2_median_per_pheno),
    width = 0.8,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  # 散点图：添加每个细胞的点，增加细节
  geom_jitter(
    aes(fill = CytoTRACE2_median_per_pheno),
    width = 0.05,  # x轴抖动范围（避免点重叠）
    height = 0,    # y轴不抖动
    alpha = 0.5,
    shape = 21,    # 圆形带边框
    stroke = 0.1,  # 边框粗细
    size = 1       # 点大小
  ) +
  # 主题：经典样式（无背景网格）
  theme_classic() +
  # y轴刻度：0~1，间隔0.2
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    # 右侧次坐标轴：显示潜能分类标签
    sec.axis = sec_axis(
      trans = ~.,  # 与主坐标轴一致
      breaks = seq(0, 1, by = 1/12),  # 间隔匹配潜能分类
      labels = c("", "Differentiated", "", "Unipotent", "", "Oligopotent", "", "Multipotent", "", "Pluripotent", "", "Totipotent", "")
    )
  ) +
  # 填充颜色映射：与潜能颜色一致
  scale_fill_gradientn(
    colors = rev(colors),
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    labels = labels
  ) +
  # 点的颜色映射（与填充一致，可选）
  scale_color_gradientn(
    colors = rev(colors),
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    labels = labels
  ) +
  # x轴标签：自动换行（避免长标签重叠）
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  # 坐标轴标签和标题
  labs(x = "Phenotype", y = "Potency score") +
  ggtitle("Developmental potential by phenotype") +
  # 主题样式：隐藏图例，调整坐标轴和刻度
  theme(
    legend.position = "None",  # 隐藏图例（箱线图的颜色由中位数决定，无需图例）
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    # 右侧次坐标轴的刻度线（仅显示潜能分类对应的刻度）
    axis.ticks.y.right = element_line(
      color = c("black", NA, "black", NA, "black", NA, "black", NA, "black", NA, "black", NA, "black")
    ),
    axis.ticks.length.y.right = unit(0.3, "cm"),  # 右侧刻度长度
    aspect.ratio = 0.8  # 箱线图宽高比（可调整）
  )

# 显示图4
print(p4)




















## Other Figures 
cyto_featureplot <- FeaturePlot(seurat_obj_cytotrace, "CytoTRACE2_Relative", pt.size = 1.5) + 
  scale_colour_gradientn(colours = (c("#9E0142", "#F46D43", "#FEE08B", 
                                      "#E6F598", "#66C2A5", "#5E4FA2")), 
                         na.value = "transparent", 
                         limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = c("0.0 (More diff.)", 
                                    "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + 
  ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  theme(legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5, 
                                  margin = margin(b = 20))) + 
  theme(aspect.ratio = 1)
print(cyto_featureplot)
ggsave("../03.Output/cytoTRACE2_featureplot.pdf", plot = cyto_featureplot, width = 10, height = 10, dpi =300)

## Other Fig
p <- ggboxplot(result_LM@meta.data, x="subcelltype", y="CytoTRACE2_Score", width = 0.6, 
                color = "black",                               #轮廓颜色
                fill="subcelltype",                            #填充
                palette = "npg",
                xlab = F,                                      #不显示x轴的标签
                bxp.errorbar=T,                                #显示误差条
                bxp.errorbar.width=0.5,                        #误差条大小
                size=1,                                        #箱型图边线的粗细
                outlier.shape=NA,                              #不显示outlier
                legend = "right")                               
###指定组比较(replace your group)
my_comparisons <- list(c("SMC_C1", "SMC_C2"), c("SMC_C1", "SMC_C3"),c("SMC_C2", "SMC_C3"))
p_box <- p + tat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")
ggsave("../03.Output/cytoTRACE2_boxplot.pdf", plot = p_box, width = 10, height = 10, dpi =300)

# ============================== Run Monocle3 ========================================
Idents(seurat_obj) <- seurat_obj$subcelltype
levels(Idents(seurat_obj))                  
# 根据CytoTRACE2/先验知识进行重设等级
seurat_obj$celltype <- factor(seurat_obj$subcelltype,levels = c("FC_C0","FC_C1","FC_C2","SMC_C1",
                                                   "SMC_C0","SMC_C2"))
Idents(seurat_obj) <- seurat_obj$subcelltype

## 提取数据
expression_matrix <- GetAssayData(seurat_obj, assay = 'RNA',slot = 'counts')
cell_metadata <- seurat_obj@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)
##
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# 归一化/预处理数据
# 这一步使用的PCA分析，dim数代表纳入的PCA数量
cds <- preprocess_cds(cds, num_dim = 25)
# 这个函数用于确认设定的dim数是否足够代表主要变异
plot_pc_variance_explained(cds)

# 可选(去批次处理)
#cds <- align_cds(cds, num_dim = 100, alignment_group = "GSE_num")

# 降维聚类，可选择UMAP、PCA或者TSNE
cds <- reduce_dimension(cds,reduction_method='UMAP',preprocess_method = 'PCA')
plot_cells(cds, label_groups_by_cluster=T ,  color_cells_by = "subcelltype")
