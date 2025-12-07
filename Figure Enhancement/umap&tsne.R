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
library(harmony)
library(patchwork)
library(RColorBrewer)
library(scRNAtoolVis)
library(CytoTRACE2)
library(monocle3)
set.seed(1234)
# 查看工作路径下的文件
list.files()

# ========================================================= Load scData for TraceAnalysis ==========================================
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)


# ========================================================= ggplot2绘制UMAP图========================================================
library(ggunchull)   # 绘制细胞群轮廓线（基于密度的凸包）
library(ggrepel)     # 避免标签重叠

##定义色板
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

## 2. 准备绘图数据（所有单细胞的metadata+UMAP坐标）
# 提取UMAP_harmony坐标并合并metadata
meta <- cbind(
  seurat_obj@meta.data,                                  # 含celltype的所有细胞metadata
  seurat_obj@reductions$umap@cell.embeddings  
)

# !!!查看UMAP坐标的实际列名（避免列名写错，比如是"UMAP_1"还是"umapharmony_1"）!!!
cat("UMAP坐标列名：", colnames(obj@reductions$umap_harmony@cell.embeddings), "\n")
# 把列名赋值给变量（后续用变量引用，避免硬编码错误）
umap_col1 <- colnames(obj@reductions$umap_harmony@cell.embeddings)[1]  # 第一列（x轴）
umap_col2 <- colnames(obj@reductions$umap_harmony@cell.embeddings)[2]  # 第二列（y轴）

### 为celltype分配颜色
names(cols) <- levels(factor(
  x = levels(meta$celltype), 
  levels = levels(meta$celltype), 
  ordered = TRUE
))

### 4. 计算每个celltype标签的位置（用变量引用列名，避免错误）
main_type_med <- meta %>% 
  group_by(celltype) %>% 
  summarise(
    # 用变量引用UMAP列名，确保和meta中的列名一致
    x = median(.data[[umap_col1]]) - 1,  # 标签x坐标（中位数-偏移）
    y = median(.data[[umap_col2]]) - 1,  # 标签y坐标（中位数-偏移）
    .groups = "drop"
  )

### 5. 绘制UMAP图（全程用变量引用列名，彻底避免列名错误）
umap_plot <- ggplot(
  data = meta,
  aes(x = .data[[umap_col1]], y = .data[[umap_col2]])  # 用变量引用x/y轴列名
) +
  # 2. 绘制单细胞点
  geom_point(aes(color = celltype), size = 0.2, show.legend = FALSE) +
  # 3. 隐藏坐标轴刻度
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  # 4. 应用色板
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  # 5. 添加celltype标签（用main_type_med的x/y列，而非原UMAP列名）
  geom_text_repel(
    data = main_type_med,
    aes(x = x, y = y, label = celltype, color = celltype),  # 用标签数据框的x/y列
    fontface = "bold", size = 6, show.legend = FALSE,
    box.padding = 0.5, max.overlaps = 30
  ) +
  # 6. 设置坐标轴标题
  labs(
    x = "UMAP_harmony 1",
    y = "UMAP_harmony 2"
  ) +
  # 7. 调整主题
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(arrow = arrow(type = "closed")),
    axis.title = element_text(hjust = 0.05, size = 12)
  )
print(umap_plot)
ggsave("../03.Output/umap.pdf", plot = umap_plot,
    width = 8,
    height = 8,
    dpi = 300,
  device = "pdf",
  limitsize = FALSE
  )
p <- p + 
    stat_unchull(
    aes(fill = celltype, color = celltype), 
    alpha = 0.1,        # 提高透明度（降低渲染负载）
    size = 0.3,         # 减细轮廓线（避免线条过重）
    lty = 1,            # 用实线（虚线计算更复杂）
    delta = 1,          # 增大平滑度（关键！降低密度计算复杂度）
    k = 5,              # 减少邻居数量（默认k=10，减小后更快）
    show.legend = FALSE
  )
p
ggsave("../03.Output/umap.pdf", plot = p, width = 8, height = 8, dpi = 300)


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
         )    # 去除所有网格线
p






































  )
