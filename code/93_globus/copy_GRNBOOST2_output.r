modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")
dim(modules) # 671322      5
colnames(modules)[1] = "predictor"
write.csv(modules, "processed-data/93_globus/MBv_GRNBOOST2_F-test-padj-05_edge-list.csv", row.names=F)


regulons = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_adj_with-logcounts-corr.csv")
dim(regulons) # 615640      5
write.csv(regulons, "processed-data/93_globus/MBv_GRNBOOST2_human-TF_edge-list.csv", row.names=F)


mod_aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
colnames(mod_aucell) = gsub("Regulon\\.for\\.", "", colnames(mod_aucell))
colnames(mod_aucell)[grep("seurat_label", colnames(mod_aucell))] = "domain_ct"
colnames(mod_aucell)[grep("smoothed", colnames(mod_aucell))] = "domain_sp"
dim(mod_aucell) # 513200     21
write.csv(mod_aucell, "processed-data/93_globus/MBv_n119_DEG-modules_AUCell.csv", row.names=T)


reg_aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(reg_aucell)[1:(grep("seurat_label", colnames(reg_aucell))-1)] = sapply(colnames(reg_aucell)[1:(grep("seurat_label", colnames(reg_aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))
colnames(reg_aucell)[grep("seurat_label", colnames(reg_aucell))] = "domain_ct"
colnames(reg_aucell)[grep("smoothed", colnames(reg_aucell))] = "domain_sp"
dim(reg_aucell) # 469396     30
## filter to regulons with at least 10 components
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
reg_subset = regulons$TF[regulons$set_size>=10]
reg_aucell = reg_aucell[, c(reg_subset, colnames(reg_aucell)[25:30])]
dim(reg_aucell) # 469396     24
write.csv(reg_aucell, "processed-data/93_globus/MBv_n109_TF-regulons_AUCell.csv", row.names=T)
