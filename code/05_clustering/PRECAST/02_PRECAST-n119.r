setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(dplyr)
})
set.seed(123)

load("processed-data/05_clustering/PRECAST/srt-list_spe_n1663_counts.Rdata")
#load("processed-data/05_clustering/PRECAST/srt-list_spe_H-M-markers_counts.Rdata")

#use customGeneList to further filter out low spcov and batch effect genes before normalization
geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1)
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% geneList$final_svgs]
cglist = rownames(avg.expr)[avg.expr$gene_name %in% geneList$qual_genes]
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% setdiff(geneList$qual_genes, geneList$all_batch_effect)]
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% setdiff(geneList$qual_genes, geneList$top.decile_low.spcov)]
length(cglist)

#prep precast
preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=cglist,
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
cat("\n\nStart model:",format(Sys.time(), tz="EST"),"\n")
PRECASTObj <- PRECAST(PRECASTObj, K=9)
#save precast object to try and get feature loadings
save(PRECASTObj, file=paste0("processed-data/05_clustering/PRECAST/precastObj_precast_k-9_n",length(cglist),".Rdata"))

# pick model (necessary but only changes things if more than 1 K) and integrate
#PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
#seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

#save(seuInt, file=paste0("processed-data/05_clustering/PRECAST/srt_precast_k-9_n",length(cglist),".Rdata"))
#cat("\n\nSaved to:",paste0("processed-data/05_clustering/PRECAST/srt_precast_k-9_n",length(cglist),".Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
