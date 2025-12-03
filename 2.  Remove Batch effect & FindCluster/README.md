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
-----RemoveBatchEffect_harmony.R<br>
--03.Output<br>
-----../Seurat_obj_resolution_res1/<br>
-----../Seurat_obj_resolution_res2/<br>
-----../Seurat_obj_cellamrkers_fig/<br>

### 01.Run Harmony
First, we perform *`NormalizeData`*, *`FindVariableFeatures`*, *`CellCycleScoring`*, *`ScaleData`* and *`RunPCA`*.<br>
Next, we should remove *batch effect* of our multiple samples. We use *`harmony`* to operate.<br>

### 02.Run sc_resolutionfinder.R
In this part, we will find the best resolution with *Function* *`sc_resolutionfinder.R`* .

### 03.Determined Best-Resolution & Cellmarekrs
Draw figures of cellmarkers.<br>

### 04. 
Save *`seurat_obj_merge`* as *`seurat_obj_merge_qc`* with format of .qs or .rds. <br>
