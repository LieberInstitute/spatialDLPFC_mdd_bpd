setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(dplyr)
})
set.seed(123)

load("processed-data/05_clustering/PRECAST/srt-list_spe_n1663_counts.Rdata")
#load("processed-data/05_clustering/PRECAST/srt-list_spe_H-M-markers_counts.Rdata")

#remove problem area spots and 3MAD outlier spots that are in low UMI cluster
flags = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=1)
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
flag.cluster = cbind(flags[rownames(cdata),], cdata[,c("precast_k9_1663","precast_k9_1663_f")])
flag.cluster$seurat_key = paste(flag.cluster$sample_id, rownames(flag.cluster), sep="_")

flag.cluster$discard_revised = flag.cluster$problem_area_flag
cat("\nSpots excluded b/c classified as problem area:\n")
table(flag.cluster$discard_revised)
#FALSE   TRUE 
#526409   8839

flag.cluster$discard_revised = ifelse(flag.cluster$genes_3MAD.outlier_binary==T & flag.cluster$precast_k9_1663_f=="low UMI",
                                      T, flag.cluster$discard_revised)
cat("\nSpots excluded b/c classified as problem area OR 3MAD outlier that was annotated to low UMI cluster:\n") 
table(flag.cluster$discard_revised)
#FALSE   TRUE 
#523352  11896 

keep.spots = flag.cluster$seurat_key[!flag.cluster$discard_revised]
cat("\nDouble check number of spots to keep:\n")
length(keep.spots)

#filter srt.sets
cat("\nNumber of spots before removing problem areas and 3MAD outliers:\n")
rowSums(sapply(srt.sets, dim))[2]
srt.sets = lapply(srt.sets, function(x) x[,colnames(x) %in% keep.spots])
cat("\nNumber of spots after removing problem areas and 3MAD outliers:\n")
rowSums(sapply(srt.sets, dim))[2]

#use customGeneList to further filter out low spcov and batch effect genes before normalization
#geneList = readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
#avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1)
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% geneList$final_svgs]
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% geneList$qual_genes]
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% setdiff(geneList$qual_genes, geneList$all_batch_effect)]
#cglist = rownames(avg.expr)[avg.expr$gene_name %in% setdiff(geneList$qual_genes, geneList$top.decile_low.spcov)]
#length(cglist)

#prep precast
preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=rownames(srt.sets[[1]]),#cglist,
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
cat("\n\nStart model:",format(Sys.time(), tz="EST"),"\n")
PRECASTObj <- PRECAST(PRECASTObj, K=9)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

#save(seuInt, file=paste0("processed-data/05_clustering/PRECAST/srt_precast_k-12_n",length(cglist),".Rdata"))
#cat("\n\nSaved to:",paste0("processed-data/05_clustering/PRECAST/srt_precast_k-12_n",length(cglist),".Rdata"))
save(seuInt, file="processed-data/05_clustering/PRECAST/srt_no-problem-areas-outliers_precast_k-9_n1663.Rdata")
cat("\n\nSaved to: processed-data/05_clustering/PRECAST/srt_no-problem-areas-outliers_precast_k-9_n1663.Rdata\n")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
