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

## Immune Cell
### T/NK
#### T/NK Alltype
```R
TNK_markers <- list(T_cell = c("CD3D","CD3E","CD3G"),
                    CD4_T  = c("CD4","IL7R"),
                    CD8_T  = c("CD8A","CD8B"),
                    γδ_T   = c("TRDC","TRGC1","TRGC2"),
                    NK_T   = c("ZBTB16","KLRB1"),
                    NK_cell= c("NCAM1","KLRD1","NKG7"),
                    CD56dim_NK_Cytotox = c("FCGR3A","PRF1","GZMB"),
                    CD56bright_NK_reg  = c("XCL1","XCL2"))
```
#### CD4/8 Tcell Subtype
```R
CD4_T_sub_markers <- list(
                 CD4_Naive = c("CCR7","LST1","IL7R"),
                 CD4_Tcm   = c("SELL"),
                 CD4_Tem   = c("CD45RO","GZMB","GZMK"),
                 Th1       = c("TBX21","IFNG","CXCR3"),
                 Th2       = c("GATA3","IL4","IL5","CCR4"),
                 Th17      = c("IL17A","RORC"),
                 Treg      = c("FOXP3","IL2RA","CTLA4"))
```

```R
CD8_T_sub_markers <- list(CD8_Naive = c("CCR7","SELL","TCF7"),
                          CD8_Tcm   = c("IL7R","CCR7"),
                          CD8_Tem   = c("GZMB","GZMK","PRF1"),
                          CD8_T_Effector = c("IFNG"),
                          CD8_T_resident = c("CD69,"ITGAE","CXCR6"),
                          CD8_Tex   = c("PDCD1","LAG3","TIGIT","TOX"))
                 




















