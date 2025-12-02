# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!
## 0. Cellranger
In this field, I don not know about this app.

## 1. Quality Control
In this part, I will process quality control for my scRNA-seq data. 
首先我们需要明确文件存放结构
>./01.QC/
>>../01.Rawdata
>>../02.Code
>>../03.Output
>>>../Seurat_obj/
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
