# scTenifoldKnk Tutorial
## Python
### 导入模块
```py
from scTenifold.data import get_test_df
from scTenifold import scTenifoldKnk
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
```
### 准备转录组数据
数据形式应该整理为（基因 x 细胞）的形式<br>
```py
df = get_test_df(n_cells=200)
df
```
### 执行虚拟敲除（注意这步比较耗时，耐心等待）
```py
sc = scTenifoldKnk(data=df,
                   ko_method="default",
                   ko_genes=["NG-568","NG-1"],  
                   qc_kws={"min_lib_size": 10, "min_percent": 0.001},
                   )
result = sc.build()
```
> data：传入转录组数据<br>
> ko_method: 设定敲除基因的方法, 这里不用改，使用“default”默认即可<br>
> ko_genes：指定要虚拟敲除的基因名。这里是列表的形式，所以可以传入多个基因，每个基因用逗号隔开。本例我同时敲除了两个基因：NG-568和NG-1<br>
> min_lib_size: 指单个基因在所有细胞或样本中的最小总表达量。如果某个基因的总表达量小于这个阈值，会被过滤掉。这里 10 表示总表达量低于 10 的基因不参与网络构建<br>
> min_percent：指基因在细胞中至少有多少比例是非零表达才保留。例如 0.001 表示至少在 0.1% 的细胞中表达，否则该基因会被过滤掉。这个参数用于去除几乎不表达的基因，避免噪声影响网络<br>
### 查看结果
```py
result
```
* 列名含义解释<br>
> Gene：基因名<br>
> Distance：这个基因在虚拟敲除（KO）网络与原始网络中的扰动程度（网络距离）。数值越大，说明该基因受KO影响越大<br>
> boxcox-transformed distance：为了满足统计假设，对原始 Distance 做 Box-Cox 变换后的值，通常用于后续Z统计量计算，使数据更接近正态分布<br>
> Z：标准化后的扰动量：衡量该基因的距离相对于其他基因的离群程度。Z 越大（正或负）表示扰动越明显<br>
> FC：反映虚拟敲除后基因表达或扰动的变化倍数<br>
> p-value: 原始P值<br>
> adjusted p-value: 多重假设校正后的P值<br>

### 绘图（绘制上调及下调幅度排名前10的基因）
```py
sns.set(style="whitegrid", context="talk")

ko_genes = ["NG-1", "NG-568"] #这里替换成你所敲除的基因名
filtered_result = result[~result["Gene"].isin(ko_genes)]

z_pos_top = filtered_result[filtered_result["Z"] > 0].sort_values("Z", ascending=False).head(10)
z_neg_top = filtered_result[filtered_result["Z"] < 0].sort_values("Z").head(10)

top_genes = pd.concat([z_pos_top, z_neg_top])
z_values = top_genes["Z"].values
genes = top_genes["Gene"].values
colors = ['#d73027' if z > 0 else '#4575b4' for z in z_values]

plt.figure(figsize=(10,7))
plt.barh(genes, z_values, color=colors, height=0.6)
plt.xlabel("Z-score", fontsize=14)
plt.ylabel("")
plt.title("Top 10 Up- and Down-regulated Genes (excluding KO genes)", fontsize=16, fontweight='bold')
plt.gca().invert_yaxis()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(axis='x', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
```
## R
