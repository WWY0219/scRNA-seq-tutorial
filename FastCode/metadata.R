# ====================================批量替换meta.data信息==================================
name_mapping <- c(
  "Normal" = "MED12-Positive",
  "Leiomyoma" = "MED12-Negative"
)
smc$group <- recode(smc$group, !!!name_mapping)
