# =================================== Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")  #replace your workspace
getwd()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(infercnv)
library(CopyKAT)
library(dplyr)
library(ggsignif)          
library(circlize)
library(ComplexHeatmap)
set.seed(1234)

# =================================== Load scData with celltype_major  ===================================================
surat_obj <- qread("../seurat_obj.qs")
Idents(seurat_obj) <- "clusters_res0.5"
##绘制样本和cluster的umap图

# =================================== Run Infercnv  ======================================================================
###构建矩阵表达数据
counts <- GetAssayData(seurat_obj, assay ="RNA", layer = 'counts')  
###细胞注释文件（分组信息）
anno <_ data.frame(
    celltype.group = seurat_obj$seurat_cluster,
    row.names =rownames(USOO_pbject@meta.data)
    )
head(anno)
gene_order <- "../01.Data/Infercnv/genev37.txt"
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,               #矩阵表达数据
                                    annotations_file = anno,                  #细胞注释文件
                                    delim="\t",
                                    gene_order_file = gene_order,             #基因注释文件
                                    min_max_counts_per_cell = c(100, +Inf),
                                    ref_group_names = c("1", "16"))           #参考细胞类型

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,                                   #use 1 for smart-seq, 0.1 for 10x-genomics 
                             out_dir = "../03.Output/", 
                             cluster_by_groups = T,                          #控制是否按预先定义的细胞分组（如细胞类型 / 样本）进行聚类
                             HMM = FALSE,
                             denoise = TRUE,                                 #决定是否对数据进行去噪处理
                             num_threads = 8,
                             rite_expr_matrix = T                            #输出处理后的矩阵
                             )

save(infercnv_obj, file = "infercnv_obj.rdata")

# =================================== Infercnv Score ======================================================================
###计算CNV score：(表达偏差-1)^2的均值（偏差越大，平方后值越大）
expr2 <- (expr - 1)^2  # 偏差平方（中心化到1，消除正负偏差抵消）
CNV_score <- as.data.frame(colMeans(expr2))  # 每个细胞的平均偏差平方（即CNV score）
colnames(CNV_score) <- "CNV_score"  # 列名改为CNV_score
CNV_score$cell <- rownames(CNV_score)  # 添加细胞ID列

# 2. 合并CNV score与聚类结果（便于后续分析）
kmeans_df$cell <- rownames(kmeans_df)  # 聚类结果添加细胞ID列
CNV_score <- CNV_score %>% inner_join(kmeans_df, by = "cell")  # 按细胞ID合并

##按细胞类别看 CNV score差异
p_cnv <- ggplot(CNV_score, aes(x = class, y = CNV_score, fill = class) ) +
geom_violin(alpha = 0.4, scale = "width") +                                # scale="width"确保小提琴宽度一致（更美观）
  stat_boxplot(
    geom = "errorbar",
    position = position_dodge(width = 0.1),
    width = 0.1,
    color = "black"  
  ) +
geom_boxplot(
    alpha = 0.5,
    outlier.size = 0,  # 隐藏离群点（避免干扰）
    size = 0.3,
    width = 0.3,
    color = "black"   
  ) +
scale_fill_manual(
    values = color_vec  
  ) +
geom_signif(
    comparisons = lapply(setdiff(unique(CNV_score$class), "normal"), function(cls) c("normal", cls)),
                         step_increase = 0.1, map_signif_level = TRUE,textsize = 3,vjust = 0.5                                   
    ) + theme_bw() +
    
  labs(
    x = "Cell Type",
    y = "CNV Score",
    title = "CNV Score Distribution Across Cell Types"  
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # x轴标签旋转+调整大小
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),  # 标题居中
    legend.title = element_blank(),  # 隐藏图例标题
    legend.position = "none"  # 若类别标签清晰，可隐藏图例
  )
print(p_cnv)
ggsave( "../03.Output/InferCNV_analysis/CNV_score_violin.pdf",plot = p_cnv,width = 10, height = 6, dpi = 300)
##按细胞类别看 CNV score差异
p_cnv <- ggplot(CNV_score, aes(x = class, y = CNV_score, fill = class) ) +
  geom_violin(alpha = 0.4, scale = "width") +                                # scale="width"确保小提琴宽度一致（更美观）
  stat_boxplot(
    geom = "errorbar",
    position = position_dodge(width = 0.1),
    width = 0.1,
    color = "black"  
  ) +
  geom_boxplot(
    alpha = 0.5,
    outlier.size = 0,  # 隐藏离群点（避免干扰）
    size = 0.3,
    width = 0.3,
    color = "black"   
  ) +
  scale_fill_manual(
    values = color_vec  
  ) +
  geom_signif(
    comparisons = lapply(setdiff(unique(CNV_score$class), "normal"), function(cls) c("normal", cls)),
    step_increase = 0.1,                          # 垂直间距（避免标记重叠）
    map_signif_level = TRUE,                      # 显示*（建议开启，更直观）
    textsize = 3,                                 # 显著性文字大小
    vjust = 0.5                                   # 调整标记位置（避免与图形重叠）
  ) +
  theme_bw() +  # 白色背景主题
  labs(
    x = "Cell Type",
    y = "CNV Score",
    title = "CNV Score Distribution Across Cell Types"  
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # x轴标签旋转+调整大小
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),  # 标题居中
    legend.title = element_blank(),  # 隐藏图例标题
    legend.position = "none"  # 若类别标签清晰，可隐藏图例
  )
print(p_cnv)
ggsave( "../03.Output/08.infercnv_analysis/CNV_score_violin.pdf",plot = p_cnv,width = 10, height = 6, dpi = 300)
    
##按聚类看CNV score差异
p_cnv <- ggplot(CNV_score, aes(x = kmeans_class, y = CNV_score, fill = kmeans_class)) +
  geom_violin(alpha = 0.4) +
  stat_boxplot(
    geom = "errorbar",
    position = position_dodge(width = 0.1),
    width = 0.1
  ) +
  geom_boxplot(
    alpha = 0.5,
    outlier.size = 0,
    size = 0.3,
    width = 0.3
  ) +
  scale_fill_manual(
    values = ggsci::pal_npg()(10)  # 聚类配色（与热图一致）
  ) +
  theme_bw() +
  labs(
    x = "K-means Cluster",
    y = "CNV Score"
  ) +
  theme(
    legend.title = element_blank()
  )

print(p_cnv)
##聚类与细胞类别的关联：验证恶性聚类
# 计算每个聚类中各细胞类别的比例
cell.prop <- as.data.frame(
  prop.table(table(CNV_score$kmeans_class, CNV_score$class))  # 按行计算比例（聚类内各类别占比）
)
colnames(cell.prop) <- c("kmeans_class", "class", "proportion")  # 重命名列

# 绘制堆叠条形图（按比例填充）
p_prop <- ggplot(cell.prop, aes(x = kmeans_class, y = proportion, fill = class)) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill"确保每列高度为1（比例之和100%）
  ggtitle("") +  # 隐藏标题
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.5, "cm"),  # 调整坐标轴刻度长度
    guides(fill = guide_legend(title = NULL)),  # 隐藏图例标题
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # 旋转x轴标签（修正参数名axis.text.x.bottom为标准axis.text.x）
  ) +
  scale_fill_manual(
    values = ggsci::pal_igv()(length(unique(cell.prop$class)))  # 颜色数量与class类别数一致（避免警告）
  ) +
  labs(
    x = "K-means Cluster",
    y = "Proportion"
  )

# 显示图形
print(p_prop)
#定义恶性细胞：整合所有信息标记细胞类型
# 1. 添加恶性细胞标签（Type列）
CNV_score <- CNV_score %>% 
  dplyr::mutate(
    Type = case_when(
      # 条件1：非6、10聚类（正常聚类）+ 肿瘤类别 → 恶性上皮（Malignant Epi）
      !(kmeans_class %in% c(6, 10)) & 
        class %in% c("Basal Epithelial_Tumor", "Intermediate_Tumor", "luminal Epithelial_Tumor") ~ "Malignant Epi",
      # 条件2：非6、10聚类 + 正常类别 → 正常上皮（Normal Epi）
      !(kmeans_class %in% c(6, 10)) & 
        class == "normal" ~ "Normal Epi",
      # 条件3：6、10聚类（正常聚类）→ 正常上皮（Normal Epi）
      kmeans_class %in% c(6, 10) ~ "Normal Epi",
      # 其他情况：标记为NA（避免错误）
      TRUE ~ NA_character_
    )
  )

# 2. 匹配细胞ID：行名=细胞ID（与Seurat对象的细胞ID一致）
rownames(CNV_score) <- CNV_score$cell

# 3. 验证CNV_score与Seurat对象的细胞ID交集（避免整合时缺失）
common_cells <- intersect(rownames(scRNA_Epi), rownames(CNV_score))

# 4. 将CNV信息（score、Type）整合回Seurat对象
scRNA_Epi <- AddMetaData(
  scRNA_Epi, 
  metadata = CNV_score  # 整合CNV_score、kmeans_class、Type等信息
)
      
# =================================== Infercnv HeatMap ======================================================================
###infercnv.references.txt: 参考细胞的CNV分值矩阵，行名是基因名，列名是细胞名
ref_score  <- read.table('infercnv.references.txt')
ref_score[1:4,1:4]
###infercnv.observations.txt: 受试细胞的CNV分值矩阵，行名是基因名，列名是细胞名
obs_score  <- read.table('infercnv.observations.txt')
obs_score[1:4,1:4]
infercnv_obj <- readRDS("/data/person/g5/wangwy/USOO/scRNA-seq-US/03.CellAnnotation/03.Output//run.final.infercnv_obj")
##取出infercnv评分数据
expr <- infercnv_obj@expr.data
##正常参考细胞的索引
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- c(normal_loc$`5`,normal_loc$`16`)
normal_loc
##待测细胞的索引
test_loc <- infercnv_obj@observation_grouped_cell_indices 
test_loc <- unlist(test_loc)
test_loc <- unname(test_loc)
test_loc[1:4]
# 创建细胞类别注释数据框（整合正常/肿瘤细胞信息）
anno.df <- data.frame(
  cell = c(colnames(expr)[normal_loc], colnames(expr)[test_loc]),  # 所有细胞ID
  class = c(
    rep('normal', length(normal_loc)),         # 正常细胞类别
    anno[colnames(expr)[test_loc], ]           # 肿瘤细胞类别
  )
)
# 查看前6行
head(anno.df)
gn <- rownames(expr)
geneFile <- read.table(
    "/data/person/g5/wangwy/USOO/scRNA-seq-US/03.CellAnnotation/TumorCell/InferCNV/01.Data/geneLocate.txt",
    header = FALSE, sep = "\t",
    stringsAsFactors = FALSE
)  
rownames(geneFile) <- geneFile$V1
common_genes <- intersect(gn, geneFile$V1) 
ub_geneFile <- geneFile[common_genes, ]
expr <- expr[common_genes,]
#k-means聚类
set.seed(1234)
kmeans.result <- kmeans(t(expr),
                        10)
kmeans_df <- data.frame(
    kmeans_class = kmeans.result$cluster,
    cell=names(kmeans.result$cluster)
    ) %>%
    inner_join(anno.df, by="cell")%>%
    arrange(kmeans_class) %>%
    column_to_rownames("cell") %>%
    mutate(kmeans_class = as.factor(kmeans_class))
table(kmeans_df$kmeans_class,kmeans_df$class)
class <- unique(kmeans_df$class)[-11]
class
str(class)
# ========== 1. 重新整理颜色向量 ==========
# 提取kmeans_df的类别信息
unique_classes <- unique(kmeans_df$class)
unique_kmeans <- unique(kmeans_df$kmeans_class)

# 构建class列的颜色向量（命名与unique_classes完全一致）
normal_color <- "blue"
other_colors <- unlist(lapply(c("Paired", "Set3", "Pastel1"), 
                               function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)))[1:(length(unique_classes)-1)]
color_vec <- c(normal_color, other_colors)
names(color_vec) <- unique_classes

# 构建kmeans_class列的颜色向量（命名与unique_kmeans完全一致）
color_v <- brewer.pal(10, "Paired")[1:length(unique_kmeans)]
names(color_v) <- as.character(unique_kmeans)

# ========== 2. 定义左侧注释 ==========
left_anno <- rowAnnotation(
    df = kmeans_df[, c("class", "kmeans_class")],  # 仅保留需要注释的列
    col = list(
        class = color_vec,       # 命名匹配class列的所有类别
        kmeans_class = color_v   # 命名匹配kmeans_class列的所有类别
    ),
    show_annotation_name = TRUE  # 可选：显示注释名称
)
pdf("../03.Output/infercnv_chr.pdf", width = 12, height = 12)

# 关键：按热图基因顺序筛选geneFile
heatmap_genes <- rownames(expr)
sub_geneFile <- geneFile[heatmap_genes, ]

# 构建热图
ht <- Heatmap(
    t(expr)[rownames(kmeans_df), ],  
    col = colorRamp2(
        c(0.8, 1, 1.2), 
        c("#377EB8", "#F0F0F0", "#E41A1C")  
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,  
    show_column_names = FALSE,
    show_row_names = FALSE,  
    column_split = factor(
        sub_geneFile$V2,  # 修正：使用匹配的sub_geneFile$V2
        levels = paste0("chr", 1:22)  
    ),
    column_gap = unit(2, "mm"),  
    heatmap_legend_param = list(
        title = "Modified expression",  
        direction = "vertical",  
        title_position = "leftcenter-rot",  
        at = c(0.8, 1, 1.2),  
        legend_height = unit(3, "cm")  
    ),
    top_annotation = top_anno,  
    left_annotation = left_anno,  
    row_title = NULL,
    column_title = NULL  
)

draw(ht, heatmap_legend_side = "right")
dev.off()  # 注意：这里需要加括号！



















