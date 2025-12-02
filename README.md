# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!
## 0. Cellranger
In this field, I don not know about this app.

## 1. Quality Control
In this part, I will process quality control for my scRNA-seq data. 
首先我们需要明确文件存放结构
>./01.QC/
>├── data/                 # 原始数据/预处理数据目录<br>
 │   ├── raw/              # 原始测序数据（如fastq、count矩阵）<br>
 │   └── processed/        # 清洗/归一化后的数据（如RDS、h5ad）<br>
 ├── scripts/              # 分析脚本目录<br>
 │   ├── 01_preprocess.R   # 数据预处理脚本<br>
 │   ├── 02_cluster.R      # 细胞聚类脚本<br>
 │   └── 03_trajectory.R   # 轨迹分析脚本<br>
 ├── results/              # 结果输出目录<br>
 │   ├── figures/          # 可视化图表（PDF/PNG）<br>
 │   ├── tables/           # 统计表格（CSV/Excel）<br>
 │   └── rds/              # 保存的Seurat/Monocle对象<br>
 ├── docs/                 # 文档目录（如实验方案、说明文档）<br>
 ├── README.md             # 项目说明文档<br>
 └── .gitignore            # Git忽略文件配置<br>
### 01.Run sc_run.q
对单样本进行单独质控，包括去除线粒体相关基因、红细胞相关基因和核糖体相关基因。样本命名应为不添加"_"的字符(*eg.Seuratobj*)，输出文件后应更名为*seuratobj_qc*.<br>
根据每个样本的相应表达进行质控并筛选细胞<br>
输出文件为*seuratobj_qc*,进行下一步操作<br>

### 02.Run sc_doublefinder.R
同样进行单样本去除双细胞，通过*Doublefinder*去除质量低下的双细胞
### 03.Run sc_decontx.R
使用*Decontx*去除环境RNA污染（可选）

## 2. Data integration and FindCluster
在这一部分，我将对质控完成的数据进行多样本整合，并进行聚类
### 04.Harmony
### 05.Annotation

## 3. Cell Annotation

## 4. Subtype Annotation

## 5. Infercnv Analysis

## 6. Cellchat Analysis

## 7. Cytotrace and Monocle Analysis

## 8. 
