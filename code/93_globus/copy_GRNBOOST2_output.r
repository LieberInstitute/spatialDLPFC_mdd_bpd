modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr.csv")
dim(modules) # 671322      5
colnames(modules)[1] = "predictor"
write.csv(modules, "processed-data/93_globus/MBv_GRNBOOST2_F-test-padj-05_edge-list.csv", row.names=F)


regulons = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_adj_with-logcounts-corr.csv")
dim(regulons) # 615640      5
write.csv(regulons, "processed-data/93_globus/MBv_GRNBOOST2_human-TF_edge-list.csv", row.names=F)
