# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!
## 0. Cellranger
In this field, I don not know about this app.

## 1. Quality Control
In this part, I will process quality control for my scRNA-seq data. 

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
质控后的单细胞数据整合后进行通过*sc_harmony.R*去除批次效应
### 05.Finder Best Resolution

## 3. Cell Annotation & Subtype Annotation
在这一部分将对最适Resolution的各个cluster进行细胞注释相关操作
### 06. Major CellAnnotation
通过各谱系marker进行注释
### 07. Subtype Annotation
* T/NK cell
* Monoctye
* Fibroblast 

## 4. Infercnv Analysis

## 5. Cellchat Analysis

## 6. Cytotrace and Monocle Analysis
### CytoTRCAE2 
### Monocle2 & Monocle3
### slingplot
