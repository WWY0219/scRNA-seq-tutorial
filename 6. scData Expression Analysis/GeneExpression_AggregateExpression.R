# ============================================ Prepare Environment ===================================================
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
library(clustree)
library(cowplot)
library(stringr)
library(ggsci)
library(pheatmap)
library(ggrepel)
library(tidyr)
library(rstatix)
set.seed(1234)
list.files()
dir.create("../03.Output/")


# ================================================ Load Data =========================================================
seurat_obj <- qread("seurat_obj.qs")
# ============================================== AggregateExpression =============================================
## ------------------------------data处理----------------------------------
av <-AggregateExpression(seurat_obj,
                         group.by = c("orig.ident","group","celltype"),
                         assays = "RNA",
                         slot = "data",                              # counts适用于DEseq2；data适用于可视化、相对表达量比较                       
                         return.seurat = FALSE) 
## av$RNA的列名格式为：orig.ident_group_celltype（由group.by的顺序决定，如"Sample1_Control_MSC"）
av=as.data.frame(av[[1]])

foxp1_agg <- as.data.frame(av["FOXP1", , drop = FALSE])              # 提取FOXP1行
colnames(foxp1_agg) <- gsub("\\.", "_", colnames(foxp1_agg))         # 若列名有小数点，替换为下划线（可选）

## 将列名拆分为orig.ident、group、celltype（关键：匹配group.by的顺序）
foxp1_df <- foxp1_agg %>%
  t() %>%                                                            # 转置：行=orig.ident_group_celltype，列=FOXP1
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group_info") %>%                 # 列名转为列
  separate(
    col = group_info, 
    into = c("orig.ident", "group", "celltype"),                     # 按分隔符拆分（默认分隔符：下划线/连字符/点，需匹配实际列名）
    sep = "_",                                                       # 若列名分隔符是其他（如"-"），改为sep = "-"
    remove = TRUE
  ) %>%
  rename(FOXP1_expression = FOXP1) %>%                               # 重命名表达量列
  mutate(
    # 转换为因子（便于后续分组和可视化）
    orig.ident = factor(orig.ident),
    group = factor(group),
    celltype = factor(celltype),
    # 转换表达量为数值型（避免字符型）
    FOXP1_expression = as.numeric(FOXP1_expression)
  )

## 检查数据结构（确保分组信息正确）
str(foxp1_df)
table(foxp1_df$celltype, foxp1_df$group)                              # 查看每个细胞群的group分布

## ----------------------------------统计检验---------------------------------------
### 按celltype分组，做每个细胞群内的group差异检验
### 结果存储为列表，每个元素是一个细胞群的检验结果
foxp1_diff_list <- foxp1_df %>%
  group_by(celltype) %>%
  group_split() %>%                                    # 按celltype拆分为列表
  lapply(function(df) {
    # 跳过样本数不足的细胞群
    if (n_distinct(df$group) < 2) {
      message(paste("细胞群", unique(df$celltype), "的group数量<2，跳过检验"))
      return(NULL)
    }
    
    # 整体检验（两组用wilcox.test，多组用kruskal.test/anova）
    group_n <- n_distinct(df$group)
    if (group_n == 2) {
      # 两组：非参数Wilcoxon检验（推荐，适用于小样本/非正态分布）
      global_test <- df %>% wilcox_test(FOXP1_expression ~ group)
    } else {
      # 多组：非参数Kruskal检验（替代ANOVA，无分布假设）
      global_test <- df %>% kruskal_test(FOXP1_expression ~ group)
      # 多组：若数据符合正态分布，可用ANOVA：df %>% anova_test(FOXP1_expression ~ group)
    }
    
    # 两两比较（多组时，做group间的两两Wilcoxon检验，带多重检验校正）
    if (group_n > 2) {
      pairwise_test <- df %>%
        wilcox_test(FOXP1_expression ~ group, p.adjust.method = "BH") %>%  # BH校正减少假阳性
        add_significance("p.adj")                                          # 添加显著性标记（*/* */***）
    } else {
      pairwise_test <- global_test %>% add_significance("p")               # 两组直接用整体检验结果
    }
    
    # 整理结果：添加细胞群名称
    result <- list(
      celltype = unique(df$celltype),
      global_test = global_test,
      pairwise_test = pairwise_test
    )
    return(result)
  })

### 过滤空结果，合并为数据框（便于查看和保存）
foxp1_diff_df <- bind_rows(
  lapply(foxp1_diff_list, function(x) {
    if (!is.null(x)) {
      x$pairwise_test %>% mutate(celltype = x$celltype)
    }
  })
)

### 查看差异检验结果（按细胞群和group分组）
print(foxp1_diff_df)
### 保存结果到CSV
write.csv(foxp1_diff_df, "FOXP1_diff_by_celltype_and_group.csv", row.names = FALSE)

## ----------------------------------小提琴图---------------------------------------
### 先计算每个细胞群内FOXP1的最大表达量（用于统一Y轴或标签位置）
foxp1_max_by_celltype <- foxp1_df %>%
  group_by(celltype) %>%
  summarize(max_expr = max(FOXP1_expression, na.rm = TRUE)) %>%
  ungroup()

### 合并最大表达量到主数据框
foxp1_df <- foxp1_df %>%
  left_join(foxp1_max_by_celltype, by = "celltype")

### 绘制分面箱线图（x=group，y=FOXP1，facet=celltype）
p_facet <- ggboxplot(
  data = foxp1_df,
  x = "group",
  y = "FOXP1_expression",
  fill = "group",
  palette = "jco",
  add = "jitter",
  facet.by = "celltype",
  scales = "free_y",
  add.params = list(alpha = 0.5)
) +
  # 关键修改：移除自定义函数，直接传入字符串（或让函数自动选择）
  stat_compare_means(
    aes(group = group),
    method = ifelse(n_distinct(foxp1_df$group) == 2, "wilcox.test", "kruskal.test"),
    label = "p.signif",
    size = 3,
    color = "black"
  ) +
  labs(
    x = "Group",
    y = "FOXP1 Expression (Aggregated)",
    title = "FOXP1 Expression in Each Cell Type by Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_facet)
ggsave("FOXP1_facet_by_celltype.pdf", plot = p_facet, width = 12, height = 8, dpi = 300)
ggsave("FOXP1_facet_by_celltype.pdf", plot = p_facet, width = 12, height = 8, dpi = 300)
