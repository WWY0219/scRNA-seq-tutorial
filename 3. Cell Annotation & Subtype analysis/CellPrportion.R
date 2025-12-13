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
  pivot_wider(names_from = subtype, values_from = count, values_fill = 0) # 缺失值填充为0
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

# 提取显著性标注的位置和标签（用于绘图）
# 定义分组顺序（NE、EU、OMA）
group_order <- c("NE", "EU", "OMA")
# 生成比较组合
comparisons <- combn(group_order, 2, simplify = FALSE)
# 获取每个比较的校正后p值和显著性标签
sig_labels <- dunn_test %>%
  mutate(
    # 统一group1和group2的顺序，与comparisons匹配
    group1 = factor(group1, levels = group_order),
    group2 = factor(group2, levels = group_order),
    comp = paste(pmin(as.character(group1), as.character(group2)), pmax(as.character(group1), as.character(group2)), sep = " vs ")
  ) %>%
  select(comp, p.adj.signif)

