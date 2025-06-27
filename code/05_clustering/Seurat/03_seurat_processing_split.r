setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
})
set.seed(123)

load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control.Rdata")

#subset to only features that passed MBv pseudobulk gene filters
cat("\nSubset to only genes that passed MBv pseudobulk gene filters...\n")
var = read.csv("processed-data/05_clustering/Seurat/var_SZBDMulti-Seq_filtered.csv")
stopifnot(identical(seu_con[["RNA"]]@meta.data$featureid, var$featureid))
seu_con[["RNA"]]@meta.data$featurekey = var$featurekey
cat("\nOriginal number of SZBDMulti-seq genes:",nrow(var),"\n")

rdata = read.csv("processed-data/06_pseudobulk/pseudobulk-sample-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)
cat("\nOriginal number of MBv pseudobulk genes:",nrow(rdata),"\n")

seu_con = seu_con[seu_con[["RNA"]]@meta.data$featureid %in% rownames(rdata),]
cat("\nNumber of MBv-filtered SZBDMulti-seq genes:",nrow(seu_con),"\n")
seu_con

#now split to glia (with a handful of neurons for transfer purposes) and neurons
seu_con$azimuth_broad = factor(seu_con$azimuth, levels=c("Astro","Oligo","OPC",
                                                         "Micro","PC","SMC","Immune","VLMC","Endo",
                                                         "Pvalb","Lamp5","Vip","Sst","Pax6","Chandelier","Sncg","Lamp5 Lhx6","Sst Chodl",
                                                         "L2/3 IT","L4 IT",
                                                         "L5 IT","L5 ET","L5/6 NP",
                                                         "L6 IT","L6 IT Car3","L6b","L6 CT"),
                               labels=c("Astro","Oligo","OPC",
                                        "Vasc/Immune","Vasc/Immune","Vasc/Immune","Vasc/Immune","Vasc/Immune",
                                        "Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb","Inhb",
                                        "L2/3","L4",
                                        "L5","L5","L5",
                                        "L6","L6","L6","L6"))

cat("\nSubset to glia and inhb. neurons...\n")
seu_glia = seu_con[,seu_con$azimuth_broad %in% c("Astro","Oligo","OPC","Vasc/Immune","Inhb")]
dim(seu_glia)

cat("\nSubset to excitatory neurons...\n")
seu_nrn = seu_con[,seu_con$azimuth_broad %in% c("L2/3","L4","L5","L6")]
dim(seu_nrn)

#cat("\nSubset to glia (with random 5k neurons)...\n")
#sub_nrn = sample(colnames(seu_nrn), size = 5000, replace=F)
#sub_glia = seu_con[,c(colnames(seu_glia), sub_nrn)]
#dim(sub_glia) #11472 33524



cat("\n\nSCT (neurons)...\n")
format(Sys.time())
seu_nrn <- SCTransform(seu_nrn, ncells=5000)

save(seu_nrn, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_exc-nrn.Rdata")
cat("\nSaved checkpoint file to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_exc-nrn.Rdata\n\n")



cat("\n\nSCT (glia)...\n")
format(Sys.time())
seu_glia <- SCTransform(seu_glia, ncells=5000)

save(seu_glia, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_glia-inhb.Rdata")
cat("\nSaved checkpoint file to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_glia-inhb.Rdata\n\n")



cat("\n\nPCA and UMAP (neurons)...\n")
format(Sys.time())
seu_nrn <- RunPCA(seu_nrn)
seu_nrn <- RunUMAP(seu_nrn, dims=1:20)

#seu_con <- FindNeighbors(seu_con, dims = 1:20, k.param=50) 
#seu_con <- FindClusters(seu_con, resolution=.2) 

save(seu_nrn, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_exc-nrn.Rdata")
cat("\nSaved to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_exc-nrn.Rdata\n")



cat("\n\nPCA and UMAP (glia)...\n")
format(Sys.time())
seu_glia <- RunPCA(seu_glia)
seu_glia <- RunUMAP(seu_glia, dims=1:20)

#seu_con <- FindNeighbors(seu_con, dims = 1:20, k.param=50) 
#seu_con <- FindClusters(seu_con, resolution=.2) 

save(seu_glia, file="processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_glia-inhb.Rdata")
cat("\nSaved to: processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_glia-inhb.Rdata\n")


print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
