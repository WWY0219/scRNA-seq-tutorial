# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!

## 3. Cell Annotation & Subtype Analysis
In this part, I will process quality control for my scRNA-seq data.<br>
First, we need to clarify the file storage structure.
>./02.DataIntegeration&FindCluster/<br>
--01.Data<br>
-----Seurat_obj_merge_qc.qs<br>
-----Seurat_obj_merge_qc.rds<br>
-----/Infercnv/genev37.txt<br>
-----...<br>
--02.Code<br>
-----31. CellAnnotation.R<br>
-----32. CellAnnotation_TumorClassification.R<br>
-----33. CellAnnotation_T_NK.R<br>
-----34. CellAnnotation_Mono_Macro.R<br>
--03.Output<br>
-----../CellMarkers/<br>
-----../InferCNV/<br>
-----../CopyKAT/<br>
-----../T_NK.qs<br>
-----../Monocyte.qs<br>
-----../Tumor.qs<br>

### 01.Major Cell Annotation & Subtypecell Subset
首先我们将对第二section 找到的cluster进行大类注释。随后进行亚群提取。

### 02.TumorClassification
在这一部分我们将使用inferCNV及cnv来鉴定肿瘤细胞亚群.

### 03.T&NK cell subtype  
We should find a best resolution for scData.<br>

### 04. SMC&Fibr subtype 
After this, we will recieve scData with cluster_res. <br>

### Final Step
整合所有的细胞类型并输出为一个qs文件
