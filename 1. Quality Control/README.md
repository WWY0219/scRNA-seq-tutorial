# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!

## 1. Quality Control
In this part, I will process quality control for my scRNA-seq data. 
首先我们需要明确文件存放结构
>./01.QC/<br>
--01.Rawdata<br>
--02.Code<br>
--03.Output<br>
-----Seurat_obj_1<br>
-----Seurat_obj_1<br>
-----...

### 01.Run sc_run.q
对单样本进行单独质控，包括去除线粒体相关基因、红细胞相关基因和核糖体相关基因。样本命名应为不添加"_"的字符(*eg.Seuratobj*)，输出文件后应更名为*seuratobj_qc*.<br>
根据每个样本的相应表达进行质控并筛选细胞<br>
输出文件为*seuratobj_qc*,进行下一步操作<br>

### 02.Run sc_doublefinder.R
同样进行单样本去除双细胞，通过*Doublefinder*去除质量低下的双细胞
### 03.Run sc_decontx.R
使用*Decontx*去除环境RNA污染（可选）

