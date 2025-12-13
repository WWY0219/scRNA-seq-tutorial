# ========================================================= Prepare Environment =====================================================
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
set.seed(1234)
# 查看工作路径下的文件
list.files()

# ========================================================= Load scData for TraceAnalysis ==========================================
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)

# ========================================================= ggplot2 绘制带箭头的图 ===================================================
## Prepare data for draw
umap <- seurat_obj@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  cbind(cell_type = seurat_obj@meta.data$celltype) # 注释后的label信息改为cell_type
head(umap)

## -------------- 计算UMAP坐标最小值-----------------
min_umap1 <- min(umap$umap_1)
min_umap2 <- min(umap$umap_2)
x_end <- min_umap1 + 3
y_end <- min_umap2 + 3

## ---------------- 绘制UMAP图-----------------------
p <- ggplot(umap, aes(x = umap_1, y = umap_2, color = cell_type)) +  
  geom_point(size = 1, alpha = 1) +  
  scale_color_manual(values = cell_colors) +
  theme(
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),   
    panel.border = element_blank(),       
    axis.title = element_blank(),          
    axis.text = element_blank(),           
    axis.ticks = element_blank(),          
    panel.background = element_rect(fill = 'white'),  # 面板背景
    plot.background = element_rect(fill = "white"),   # 绘图区背景
    legend.title = element_blank(),        
    legend.key = element_rect(fill = 'white'), 
    legend.text = element_text(size = 20), 
    legend.key.size = unit(1, 'cm')        
 ) +  
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  annotate(
    "segment",
    x = min_umap1, y = min_umap2,          # 起点
    xend = x_end, yend = min_umap2,        # 终点
    colour = "black", size = 1,
    arrow = arrow(length = unit(0.3, "cm"))
  ) + 
  # 绘制Y轴箭头（用annotate替代geom_segment）
  annotate(
    "segment",
    x = min_umap1, y = min_umap2,          # 起点
    xend = min_umap1, yend = y_end,        # 终点
    colour = "black", size = 1,
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  # 标注X轴文字
  annotate(
    "text", 
    x = min_umap1 + 1.5, y = min_umap2 - 1, 
    label = "UMAP_1", color = "black", size = 3, fontface = "bold"
  ) + 
  # 标注Y轴文字
  annotate(
    "text", 
    x = min_umap1 - 1, y = min_umap2 + 1.5, 
    label = "UMAP_2", color = "black", size = 3, fontface = "bold", angle = 90
  )
print(p)
ggsave("ULM-UMAP.pdf",width=10,height=8,dpi=300)

# ======================================================= ggplot2绘制分面UMAP图 ======================================================
## ------------------------准备数据（合并group+celltype） ------------------------
umap_df <- seurat_obj@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(
    group = seurat_obj@meta.data$group,        
    celltype = seurat_obj@meta.data$celltype   
  )
head(umap_df)
celltype_unique <- unique(as.character(umap_df$celltype))

## -------------------------计算整个UMAP的最小坐标---------------------
min_umap1 <- min(umap_df$umap_1)
min_umap2 <- min(umap_df$umap_2)
x_end <- min_umap1 + 3  # X轴箭头长度
y_end <- min_umap2 + 3  # Y轴箭头长度

## -------------------------绘制分面UMAP-------------------------------
p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype)) +  
  geom_point(size = 1, alpha = 1) +  
  scale_color_manual(values = cell_colors, drop = FALSE) +
  facet_wrap(~group, ncol = 2) +
  theme(
    panel.grid.major = element_blank(),    
    panel.grid.minor = element_blank(),    
    panel.border = element_blank(),        
    axis.title = element_blank(),          
    axis.text = element_blank(),           
    axis.ticks = element_blank(),          
    panel.background = element_rect(fill = 'white'),  
    plot.background = element_rect(fill = "white"),   
    legend.title = element_blank(),        
    legend.key = element_rect(fill = 'white'), 
    legend.text = element_text(size = 12), 
    legend.key.size = unit(0.8, 'cm'),     
    strip.background = element_blank(),    
    strip.text = element_text(size = 16, face = "bold"),
    plot.margin = margin(1, 1, 2, 2, "cm")                 # 预留左下角箭头空间
  ) +  
  guides(color = guide_legend(override.aes = list(size = 5))) +  # 调整图例中点的大小
  # X轴箭头（左下角）
  annotate(
    "segment",
    x = min_umap1, y = min_umap2,          
    xend = x_end, yend = min_umap2,        
    colour = "black", size = 1.2,                             # 箭头加粗更显眼
    arrow = arrow(length = unit(0.4, "cm"), type = "closed")  # 实心箭头
  ) + 
  # Y轴箭头（左下角）
  annotate(
    "segment",
    x = min_umap1, y = min_umap2,          
    xend = min_umap1, yend = y_end,        
    colour = "black", size = 1.2,
    arrow = arrow(length = unit(0.4, "cm"), type = "closed")
  ) +
  # UMAP_1文字标注（箭头下方）
  annotate(
    "text", 
    x = min_umap1 + 1.5, y = min_umap2 - 1.2, 
    label = "UMAP_1", color = "black", size = 4, fontface = "bold"
  ) + 
  # UMAP_2文字标注（箭头左侧，竖排）
  annotate(
    "text", 
    x = min_umap1 - 1.2, y = min_umap2 + 1.5, 
    label = "UMAP_2", color = "black", size = 4, fontface = "bold", angle = 90
  )

print(p)


ggsave("umap_facet_single_arrow.pdf",plot = p,width = 16,  height = 8,
       dpi = 300, device = "pdf", bg = "white")




# ========================================================= ggplot2绘制UMAP图 ========================================================
library(ggunchull)   # 绘制细胞群轮廓线（基于密度的凸包）
library(ggrepel)     

## 定义色板
cols <- c(
  # 蓝色系（3种渐变）
  colorRampPalette(c("#1f78b4", "#a6cee3"))(3),
  # 绿色系（3种渐变）
  colorRampPalette(c("#33a02c", "#b2df8a"))(3),
  # 粉色
  "#df65b0",  
  # 橙色系（3种渐变）
  colorRampPalette(c("#fdbf6f", "#ff7f00"))(3), 
  # 基础单色
  "#696969", "#d8bfd8", "#008b00", "#fb9a99", "#8b5a2b", "#da70d6",
  # 青绿色系（2种渐变）
  colorRampPalette(c("#66cdaa", "#5f9ea0"))(2), 
  # 紫色系（6种渐变）
  colorRampPalette(c("#cab2d6", "#6a3d9a"))(6),
  # 补充单色
  "#1c6597", "#cc6805", "#c6abd4", "#569395", "#dc65d8", "#f99190", "#057605","#575757"
)

## 准备绘图数据（所有单细胞的metadata+UMAP坐标）
obj <- seurat_obj  
meta <- cbind(
  obj@meta.data,                                  # 含celltype的所有细胞metadata
  obj@reductions$umap@cell.embeddings     
)

### 查看UMAP坐标的实际列名（避免列名写错）
cat("UMAP坐标列名：", colnames(obj@reductions$umap@cell.embeddings), "\n")

### 把列名赋值给变量（后续用变量引用，避免硬编码错误）
umap_col1 <- colnames(obj@reductions$umap@cell.embeddings)[1]  # 第一列（x轴）
umap_col2 <- colnames(obj@reductions$umap@cell.embeddings)[2]  # 第二列（y轴）

### 为celltype分配颜色（确保颜色数量匹配celltype水平数）
celltype_levels <- levels(meta$celltype)
cols <- cols[1:length(celltype_levels)]  
names(cols) <- celltype_levels

### 计算每个celltype标签的位置
main_type_med <- meta %>% 
  group_by(celltype) %>% 
  summarise(
    x = median(.data[[umap_col1]], na.rm = TRUE) - 1,  # 加na.rm避免NA值
    y = median(.data[[umap_col2]], na.rm = TRUE) - 1,
    .groups = "drop"
  ) %>%
  filter(!is.na(x) & !is.na(y))                        # 过滤NA行，避免标签绘制报错

### 绘制UMAP图（优化图层顺序：先画凸包，再画点，最后画标签）
umap_plot <- ggplot(
  data = meta,
  aes(x = .data[[umap_col1]], y = .data[[umap_col2]])  
) +
  stat_unchull(
    aes(fill = celltype, color = celltype), 
    alpha = 0.1,        # 提高透明度，降低渲染负载
    size = 0.3,         # 减细轮廓线，避免线条过重
    lty = 1,            # 实线（虚线计算更复杂，提速）
    delta = 1,          # 增大平滑度，降低密度计算复杂度
    k = 5,              # 减少邻居数量，加快渲染
    show.legend = FALSE
  ) +
  geom_point(aes(color = celltype), size = 0.2, show.legend = FALSE) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  geom_text_repel(
    data = main_type_med,
    aes(x = x, y = y, label = celltype, color = cols),  
    inherit.aes = FALSE,                               # 取消继承顶层aes，避免列名冲突
    fontface = "bold", size = 6, show.legend = FALSE,
    box.padding = 0.5, max.overlaps = 30,
    min.segment.length = 0                              # 强制显示所有标签连线（可选）
  ) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # 轴箭头核心设置：arrow()参数控制“小箭头”
    axis.line = element_line(
      color = "black",        # 箭头颜色
      linewidth = 0.4,        # 轴线条粗细（匹配小箭头）
      arrow = grid::arrow(
        length = unit(0.08, "inches"),  # 箭头长度（0.08英寸=约2mm，小巧）
        width = unit(0.05, "inches"),   # 箭头宽度（可选，更小巧）
        type = "closed"                 # 闭合三角形箭头
      )
    ),
    axis.title.x = element_text(hjust = 0.95, size = 12, margin = margin(t = 5)),
    axis.title.y = element_text(hjust = 0.95, size = 12, margin = margin(r = 5)),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )
print(umap_plot)

# 保存PDF（优化保存参数）
ggsave(
  filename = "../03.Output/umap.pdf", 
  plot = umap_plot,
  width = 8,
  height = 8,
  dpi = 300,
  device = "pdf",
  limitsize = FALSE,
  units = "in"  # 明确单位为英寸
)


# ========================================================= 带小箭头的UMAP图 ========================================================
library(tidydr)
color <- c("#96F148","#ff7f00","#e5f5f9","#bebada","#df65b0","#D10000","#0000FF","#fff7fb","#fccde5",
           "#bc80bd","#d9d9d9","#ffed6f","#d6604d","#02818a","#ccecb5","#80b1d3","#fb9a99","#006837",
           "#6a3d9a")

p <- ggplot(seurat_obj@meta.data, aes(x=UMAP1, y=UMAP2, color=cell_state)) +
  geom_point(size=0.1, alpha=0.3,shape = 21, stroke = 0.9) +
  scale_color_manual(values=color) +
  guides( color = guide_legend( title = "", override.aes = list( fill=color,color = "black",stroke = 0.3, size = 4, alpha = 1), ncol = 2 )) +
  theme_dr() +   # 应用带小箭头的坐标轴主题（来自tidydr包）
  theme( panel.grid = element_blank(),
         legend.text = element_text( size = 12, face = "plain",color = "black")
         )   
p

# =============================================== 带标签的UMAP ==================================================================
# 加载必需包
library(ggplot2)
library(ggh4x)

# ===================== 关键修正1：重命名颜色向量（避免与内置函数冲突） =====================
cell_colors <- c("#919ac2","#ffac98","#70a4c8","#a5a9af","#63917d","#dbd1b4","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6")

# ===================== 关键修正2：统一列名/数据引用（避免大小写/数据框错误） =====================
# 确保umap数据框的列名是umap_1/umap_2（而非UMAP_1/UMAP_2）
# 若列名是大写，先标准化：
# colnames(umap)[colnames(umap) == "UMAP_1"] <- "umap_1"
# colnames(umap)[colnames(umap) == "UMAP_2"] <- "umap_2"

# 计算轴范围（基于umap数据框，而非df）
x_lim <- c(min(umap$umap_1), max(umap$umap_1))
y_lim <- c(min(umap$umap_2), max(umap$umap_2))

# ===================== 绘制图形（全半角空格+语法修正） =====================
p <- ggplot(umap, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 0.03, shape = 16, stroke = 0) +
  scale_color_manual(values = cell_colors) +  # 改用重命名的颜色向量
  theme_classic() +  # 使用简洁主题
  theme(
    plot.background = element_blank(),        # 移除背景
    panel.grid.major = element_blank(),       # 移除主要网格线
    panel.grid.minor = element_blank(),       # 移除次要网格线
    plot.margin = margin(5, 5, 5, 5, "mm"),   # 修正margin写法（添加边距值）
    axis.title.x = element_blank(),           # 移除x轴标题
    axis.title.y = element_blank(),           # 移除y轴标题
    axis.text = element_blank(),              # 移除坐标轴刻度标签
    axis.ticks = element_blank(),             # 移除坐标轴刻度线
    # 修正轴箭头：全半角空格+语法规范
    axis.line = element_line(
      colour = "black", 
      linewidth = 0.3,  # 替代size（ggplot2 3.4+推荐linewidth）
      arrow = arrow(length = unit(0.1, "cm"))
    ),  
    strip.background = element_rect(fill = '#e6bac5', color = NA),
    strip.placement = 'outside',
    strip.text = element_text(size = 8),
    legend.position = "none",
    aspect.ratio = 1
  ) +
  scale_x_continuous(limits = x_lim) +  # 使用提前计算的轴范围
  scale_y_continuous(limits = y_lim)

# 打印图形
print(p)




# =============================================== geom_sc_umap.r ==================================================================
source("geom_sc_umap.r")
#读取 Seurat object
pbmc=readRDS("pbmc.rds")
strip.col=c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")
## 带标签的UMAP
ggscplot(object = pbmc, reduction="umap") +
  geom_sc_umap(aes(color = seurat_annotations,
                   cluster = seurat_annotations,
                   cluster_anno = seurat_clusters),
               show.legend = F,
               add_label = T,
               add_legend = T,
               lgd_x = 1.2) +
  theme_sc(r = 0.3)+scale_color_manual(values=strip.col)

## 多列图例
ggscplot(object = pbmc, reduction="tsne") +
  geom_sc_umap(aes(color = seurat_annotations,
                    cluster = seurat_annotations,
                    cluster_anno = seurat_clusters),
                show.legend = F,
                add_label = T,
                add_legend = T,
                lgd_x = 0.5,ncol = 3,lgd_y = -0.25) +
  theme_sc(b = 0.25,x.label = "TSNE 1",y.label = "TSNE 2")

## 
#模拟分组
pbmc$groups = c(rep("Control",round(ncol(pbmc)/2)),rep("Treat",ncol(pbmc)-round(ncol(pbmc)/2)))
ggscplot(object = pbmc, reduction="umap") +
  geom_sc_umap(aes(color = seurat_annotations,
                   cluster = seurat_annotations),
               label.gp = gpar(fontsize = 8,fontface = "bold.italic")) +
  theme_sc(x.label = "UMAP 1",y.label = "UMAP 2") +
  scale_color_manual(values=strip.col)+
  facet_wrap(~groups,scales = "free")
dev.off()
