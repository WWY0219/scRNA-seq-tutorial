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
set.seed(1234)
list.files()
dir.create("../03.Output/")


# ============================================ Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")


# ============================================ Preparation Data ============================================
## 筛选所需细胞亚型
meta_data <- seurat_obj@meta.data %>%
  filter(subtype %in% c("SMC", "Fibroblast","T/NK",...)) %>%
  select(orig.ident, subtype, group)

cell_counts <- meta_data %>%
  group_by(orig.ident, group, subtype) %>%  # 分组顺序:orig.ident, group, subtype
  summarise(count = n(), .groups = "drop")  # 等价于ungroup()
cell_counts[1:4,1:4]

## 转换为宽表格
count_wide <- cell_counts %>%
  pivot_wider(names_from = subtype, 
              values_from = count, 
              values_fill = 0) # 缺失值填充为0
count_wide[1:4,]

## 计算相应细胞比例
df_plot <- count_wide %>%
  rename(                
    Naive_CD4_T = `Naive-CD4+T`,            # 替换掉非法字符
    ...
  ) %>%
  mutate(
    Th17_Treg_ratio = ifelse(Treg == 0, NA, Th17_like / Treg),
  ) %>%
  filter(!is.na(Th17_Treg_ratio)) 
df_plot[1:4,]

# ============================================ 统计检验 ============================================
## Kruskal-Wallis检验（整体差异）
kruskal_test <- df_plot %>%
  kruskal_test(Th17_Treg_ratio ~ group)
print(kruskal_test)

## Dunn检验（事后多重比较，校正方法：BH）
dunn_test <- df_plot %>%
  dunn_test(Th17_Treg_ratio ~ group, p.adjust.method = "BH") %>%
  select(group1, group2, p, p.adj, p.adj.signif)
print(dunn_test)


## 定义分组顺序
group_order <- c("NE", "EU", "OMA")

## 生成比较组合
comparisons <- combn(group_order, 2, simplify = FALSE)

## 获取每个比较的校正后p值和显著性标签
sig_labels <- dunn_test %>%
  mutate(
    group1 = factor(group1, levels = group_order),
    group2 = factor(group2, levels = group_order),
    comp = paste(pmin(as.character(group1), as.character(group2)), pmax(as.character(group1), as.character(group2)), sep = " vs ")
  ) %>%
  select(comp, p.adj.signif)

# ============================================ Visulization ========================================
# 定义自定义颜色（三组对应三种颜色，可自行调整）
group_colors <- c("NE" = "#3498db", "EU" = "#e74c3c", "OMA" = "#2ecc71")

# 绘制箱线图
p <- ggplot(df_plot, aes(x = factor(group, levels = group_order), y = Th17_Treg_ratio)) +
  # 箱线图：按group着色，调整透明度和宽度
  geom_boxplot(aes(fill = group), alpha = 0.7, width = 0.6, outlier.shape = NA) +
  # 散点：展示每个样本的具体值，避免重叠
  geom_jitter(aes(color = group), size = 3, alpha = 0.8, position = position_jitter(width = 0.2)) +
  # 标注统计显著性（ggsignif）
  geom_signif(
    comparisons = comparisons, # 比较组合
    data = df_plot,
    map_signif_level = TRUE, # 显示*/**/***（替代p值）
    y_position = seq(max(df_plot$Th17_Treg_ratio) * 1.1, max(df_plot$Th17_Treg_ratio) * 1.3, length.out = length(comparisons)), # 标注位置
    textsize = 4, # 字体大小
    vjust = 0.5 # 垂直调整
  ) +
  # 自定义颜色
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # 坐标轴和标题标签
  labs(
    title = "Distribution of Th17-like/Treg Ratio in NE, EU and OMA Groups",
    x = "Group",
    y = "Th17-like / Treg Ratio",
    fill = "Group",
    color = "Group"
  ) +
  # 主题美化
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # 标题居中、加粗
    axis.title = element_text(size = 12), # 轴标题大小
    axis.text = element_text(size = 10), # 轴标签大小
    legend.position = "none", # 隐藏图例（x轴已显示group）
    panel.grid = element_blank() # 隐藏网格线（可选）
  )

# 显示图形
print(p)



## ============================================ ggpubr-Visulization ========================================
### 定义分组顺序（确保绘图顺序为NE、EU、OMA）
group_order <- c("NE", "EU", "OMA")
df_plot$group <- factor(df_plot$group, levels = group_order) # 强制分组顺序

# 定义自定义颜色（三组对应三种颜色）
group_colors <- c("NE" = "#3498db", "EU" = "#e74c3c", "OMA" = "#2ecc71")
# ------------- 步骤1：定义组间比较组合（与原代码一致）-------------
comparisons <- combn(group_order, 2, simplify = FALSE) # 生成所有两两比较组合：NE vs EU, NE vs OMA, EU vs OMA

# ------------- 步骤2：绘制箱线图并添加统计检验 -------------
p <- ggboxplot(
  df_plot, 
  x = "group",               # x轴：分组
  y = "Th17_Treg_ratio",     # y轴：比值
  color = "group",           # 箱线图边框颜色按group分组
  fill = "group",            # 箱线图填充颜色按group分组
  palette = group_colors,    # 自定义颜色（与group_colors匹配）
  alpha = 0.7,               # 填充透明度
  width = 0.6,               # 箱线图宽度
  outlier.shape = NA,        # 隐藏箱线图的异常值点（后续用散点展示）
  add = "jitter",            # 添加散点（展示每个样本的具体值）
  add.params = list(         # 散点参数设置
    size = 3, alpha = 0.8, position = position_jitter(width = 0.2)
  ),
  legend = "right",          
  title = "Distribution of Th17-like/Treg Ratio in NE, EU and OMA Groups", # 标题
  xlab = "Group",            # x轴标签
  ylab = "Th17-like / Treg Ratio",  # y轴标签
    legend.title = "Group"
) +
   font("legend.title", color="black", face = "bold",size = 15)+
  font("legend.text", color = "black")+
  # ------------- 统计检验1：添加整体检验（如Kruskal-Wallis/ANOVA）-------------
  stat_compare_means(
    method = "anova", # 整体检验方法（非参数，适合单细胞数据）；若符合正态分布，用"anova"
    label.x = 1.5,           # 整体检验标签的x轴位置（可根据分组数调整）
    label.y = max(df_plot$Th17_Treg_ratio, na.rm = TRUE) * 1.3, # y轴位置（在最上方）
    size = 4,                # 字体大小
    face = "bold"            # 字体加粗
  ) +
  # ------------- 统计检验2：添加组间两两比较（多重比较）-------------
  stat_compare_means(
    comparisons = comparisons, # 两两比较组合
    method = "wilcox.test",      # kruskal.test
    p.adjust.method = "BH",    # 多重检验校正方法（Benjamini-Hochberg）
    map_signif_level = TRUE,   # 显示*/**/***（替代p值）
    y.position = seq(          # 每个比较的标注位置（从上到下）
      max(df_plot$Th17_Treg_ratio, na.rm = TRUE) * 1.25,
      max(df_plot$Th17_Treg_ratio, na.rm = TRUE) * 1.1,
      length.out = length(comparisons)
    ),
    size = 4                   # 显著性标签字体大小
  ) +
  #theme_pubr() + # ggpubr的默认美观主题
  border("black") +
  theme(
       axis.line = element_line(color = "black", size = 0.5), # 边线颜色+粗细
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # 标题居中、加粗
    axis.title = element_text(size = 12), # 轴标题大小
    axis.text = element_text(size = 10), # 轴标签大小
    panel.grid = element_blank() # 隐藏网格线（可选）
  )

# 显示图形
print(p)






















