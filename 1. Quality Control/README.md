# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!

## 1. Quality Control
In this part, I will process quality control for my scRNA-seq data.<br>
First, we need to clarify the file storage structure.
>./01.QC/<br>
--01.Rawdata<br>
-----Seurat_obj_1.h5<br>
-----Seurat_obj_2.h5<br>
-----...<br>
--02.Code<br>
-----QC.R<br>
--03.Output<br>
-----../Seurat_obj_1/<br>
-----../Seurat_obj_1<br>
-----../...<br>

### 01.Run sc_run.q
Perform separate QC for each individual sample, including the removal of Mit-related genes, RBC-related genes, and ribo-related genes. Name the sample (e.g., *`Seurat_obj_1`*) and name the output file *`seurat_obj_1_qc`*.<br>
Perform QC and filtering based on the corresponding expression profile of each sample, and name the resulting object *`seurat_obj_1_filtered`*.<br>

### 02.Run sc_doublefinder.R
Remove doublets from individual samples (inferior - quality doublets are eliminated using <mark>***Doublefinder***</mark>), and name the output file *`seurat_obj_1_db`*.

### 03.Run sc_decontx.R
使用<mark>***Decontx***<mark>去除环境RNA污染（可选）<br>

### Merge Files
Save *`seurat_obj_merge`* as *`seurat_obj_merge_qc`* with format of .qs or .rds. <br>

