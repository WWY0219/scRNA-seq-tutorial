# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!

## 2. Data Integeration and FindCluster
In this part, I will process quality control for my scRNA-seq data.<br>
First, we need to clarify the file storage structure.
>./02.DataIntegeration&FindCluster/<br>
--01.Data<br>
-----Seurat_obj_merge_qc.qs<br>
-----Seurat_obj_merge_qc.rds<br>
-----...<br>
--02.Code<br>
-----harmony.R<br>
--03.Output<br>
-----../Seurat_obj_resolution_res1/<br>
-----../Seurat_obj_resolution_res2/<br>
-----../...<br>

### 01.Run Harmony
Perform separate QC for each individual sample, including the removal of Mit-related genes, RBC-related genes, and ribo-related genes. Name the sample (e.g., *`Seurat_obj_1`*) and name the output file *`seurat_obj_1_qc`*.<br>
Perform QC and filtering based on the corresponding expression profile of each sample, and name the resulting object *`seurat_obj_1_filtered`*.<br>

### 02.Run sc_resolutionfinder.R
Remove doublets from individual samples (inferior - quality doublets are eliminated using <mark>***Doublefinder***</mark>), and name the output file *`seurat_obj_1_db`*.

### 03.Run sc_decontx.R
使用<mark>***Decontx***</mark>去除环境RNA污染（可选）<br>

### Merge Files
Save *`seurat_obj_merge`* as *`seurat_obj_merge_qc`* with format of .qs or .rds. <br>
