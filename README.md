# scRNA-seq-tutorial
The relevant codes used by WWY for analyzing **scRNA-seq**. Please be careful!!!
## Cellranger
In this field, I don not know about this app.
## Quality Control
In this part, I will process quality control for my scRNA-seq data. 
### 01.Run sc_run.q
对单样本进行单独质控，包括去除线粒体相关基因、红细胞相关基因和核糖体相关基因。样本命名应为不添加"_"的字符，输出文件后应更名为  
思华
### 02.Run sc_doublefinder.R
### 03.Run sc_decontx.R
## Data integration and FindCluster

