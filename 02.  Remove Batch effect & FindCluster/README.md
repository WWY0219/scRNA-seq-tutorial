# 2. Data Integeration and FindCluster

## Environment
> Conda::R432<br>

## Description of files and folders
In this part, ....<br>
First, we need to clarify the file storage structure.<br>
### 01.Data
### 02.Code
### 03.Output
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
-----../Seurat_obj_harmony.qs<br>

## Tutorial
### 01.Run Harmony
First, we perform *`NormalizeData`*, *`FindVariableFeatures`*, *`CellCycleScoring`*, *`ScaleData`* and *`RunPCA`*.<br>
Next, we should remove *batch effect* of our multiple samples. We use *`harmony`* to operate.<br>

### 02.Run sc_resolutionfinder.R
In this part, we will find the best resolution with *Function* *`sc_resolutionfinder.R`* .

### 03.Determined Best-Resolution 
We should find a best resolution for scData.<br>

### Output Profile 
After this, we will recieve scData with cluster_res. <br>
