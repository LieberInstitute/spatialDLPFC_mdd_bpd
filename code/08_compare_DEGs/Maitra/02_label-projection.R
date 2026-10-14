setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(Seurat)
})
set.seed(123)

#load objects
#cat("\nProcessing MDD males...\n")
#load("processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_male_MBv-filtered_processed-SCT.Rdata")
#dim(seu_m)
cat("\nProcessing MDD females...\n")
load("processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_female_MBv-filtered_processed-SCT.Rdata")
dim(seu_f)

## load reference
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")

#integrate
cat("\nStart integration...\n")
format(Sys.time())
anchors <- FindTransferAnchors(reference = seu_con,
       features=geneList$qual_genes,
       query = seu_f, normalization.method = "SCT",
       k.anchor=50, k.score=50, max.features=500, npcs=30, dims=1:30, mapping.score.k=T)

save(anchors, file="processed-data/08_compare_DEGs/Maitra/anchors_Maitra2023_MBv-filtered_ref-SZBD-control_query-female_qual-genes-kanchor-50-pc30.Rdata")
cat("\nAnchors saved to: processed-data/08_compare_DEGs/Maitra/anchors_Maitra2023_MBv-filtered_ref-SZBD-control_query-female_qual-genes-kanchor-50-pc30.Rdata\n")

#transfer labels
cat("\nTransfer labels...\n")
format(Sys.time())

mdd_query <- TransferData(anchorset = anchors, refdata = seu_con$seurat_low.res,
        weight.reduction = seu_f[["pca"]], dims = 1:30, k.weight=50)

write.csv(mdd_query, "processed-data/08_compare_DEGs/Maitra/results_label-transfer_MBv-filtered_ref-SZBD-control_query-female_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=T)
cat("\nSaved label transfer to: processed-data/08_compare_DEGs/Maitra/results_label-transfer_MBv-filtered_ref-SZBD-control_query-female_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv\n")


## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

