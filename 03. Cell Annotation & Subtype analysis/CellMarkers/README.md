# Cellmarker
各种细胞的Marker<br>
## Tumor



## Ovarian Aging
### Human
```R
H_Ovary_markers <- list( Oocytes =c("ZP3","TUBB8B","TUBB8","GDF9","FIGLA","OOSP2"),
                         GCs = c("GSTA1","AMH","HSD17B1","INH8B","TNNI3"),
                         TCs = c("DCN","STAR","CYP17A1","HHIP","LHCGR"),                                                                     
                         T_cell = c("CD3E", "CD3D", "CD4", "CD8A"),
                         B_cell = c("MS4A1", "CD19", "CD79A"),
                         Monocyte = c("CD14", "CD68", "CSF1R", "FCGR3A", "LYZ","ITGAM"),
                         Dcs = c("CD1C", "CLEC9A", "HLA-DRA", "CD83"),
                         NK = c("NKG7", "GNLY", "KLRD1", "NCAM1"),                      
                         Endothelial_cell = c("PECAM1", "VWF", "CDH5", "TM4SF1"),     
                         SMC = c("ACTA2", "MUSTN1"),                            
                         Epithelial_cell = c("EPCAM", "KRT8", "KRT18", "CDH1")
                    )
```

### Mouse
```R
M_Ovary_markers <- list(Oocytes =c("Zp3"),                                                                  
                     T_cell = c("CD3E", "CD3D", "CD4", "CD8A"),
                     B_cell = c("Ms4a1", "Cd19", "Cd79a"),
                     Monocyte = c("CD14", "CD68", "CSF1R", "FCGR3A", "LYZ","ITGAM"),
                     Dcs = c("CD1C", "CLEC9A", "HLA-DRA", "CD83"),
                     NK = c("NKG7", "GNLY", "KLRD1", "NCAM1"),                      
                     TCs = c("Srd5a1"),
                     GCs = c("Amh"),      
                     Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),     
                     SMC = c("ACTA2", "MYH11", "TAGLN"),                            
                     luteal_cell =c("Ptgfr"),
                     Epithelial = c("EPCAM", "KRT8", "KRT18", "CDH1")
                    )
```
