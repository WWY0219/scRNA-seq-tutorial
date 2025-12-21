# ====================================批量替换meta.data信息==================================
name_mapping <- c(
  "Normal" = "MED12-Positive",
  "Leiomyoma" = "MED12-Negative"
)
smc$group <- recode(smc$group, !!!name_mapping)

# ====================================批量替换meta.data信息==================================
seurat_obj$group[seurat_obj$orig.ident == "ne1"] <- "NE"
seurat_obj$celltype <- stringr::str_replace(seurat_obj$celltype, "^Epthelial$", "Epithelial")

seurat_obj@meta.data$orig.ident <- factor(seurat_obj@meta.data$orig.ident, 
                                          levels=c("ne1","ne2","ne3", ... ))


# ============================ Add groupdata into meta.data ===================================
## identify groupnames
healthy_samples <- c("Healthy_01", "Healthy_02")  
tumor_samples   <- c("Tumor_01", "Tumor_02")  

seurat_obj@meta.data$group <- ifelse(
  seurat_obj@meta.data$orig.ident %in% healthy_samples,         # 对照样本ID
  "Healthy",                                             
  "Tumor"                                            
)
