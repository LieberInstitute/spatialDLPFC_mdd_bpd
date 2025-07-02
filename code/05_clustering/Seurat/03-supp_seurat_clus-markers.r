setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(Seurat)
})
set.seed(123)

#load seu_con MBv-filtered to extract annotations
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
mdata = seu_con@meta.data
rm(seu_con)

#load seu_con not MBv-filtered to do FindMarkers on genes (not limited by SRT expr)
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_processed-SCT.Rdata")
stopifnot(identical(rownames(seu_con@meta.data), rownames(mdata)))

seu_con$seurat_annotated <- mdata$seurat_annotated
seu_con$seurat_low.res <- mdata$seurat_low.res

Idents(seu_con) <- "seurat_low.res"

#log2 norm RNA??
#spe <- computeLibraryFactors(spe)
#spe <- logNormCounts(spe)

markers = FindMarkers(seu_con, ident.1="L4", min.pct=.1, test.use="wilcox")
saveRDS(markers, "processed-data/05_clustering/Seurat/results_SZBDMulti-seq_control_FindMarkers_L4_SCT-wilcox.rda")
cat("\n\nFindMarkers results saved to: processed-data/05_clustering/Seurat/results_SZBDMulti-seq_control_FindMarkers_L4_SCT-wilcox.rda\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
