setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scuttle)
})
set.seed(123)
setAutoBlockSize(1e9)

#load spe
cat("\nFull spe...\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
dim(spe) 

#load domains
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))
cdata$smoothed_k9_1663 = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                                labels=c("L1","L2","L3.4","L5","L6","WM","low.UMI","Vasc","GABA"))

#load seurat labels
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res = res[rownames(cdata),]
stopifnot(identical(rownames(res), rownames(cdata)))
cdata$seurat_label = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                            labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))

#remove low UMI cluster
cat("\nRemoving low UMI domain...\n")
cdata = cdata[cdata$smoothed_k9_1663!="low.UMI",]
cdata$smoothed_k9_1663 = droplevels(cdata$smoothed_k9_1663)
dim(cdata)

#create custom clusters
cdata$custom_cluster = as.character(cdata$seurat_label)
cdata[cdata$custom_cluster=="Astro" & cdata$smoothed_k9_1663 %in% c("L1","Vasc","WM"), "custom_cluster"] = "Astro.L1"
cdata[cdata$custom_cluster=="Astro" & !cdata$smoothed_k9_1663 %in% c("L1","Vasc","WM"), "custom_cluster"] = "Astro.Nrn"
cdata[cdata$custom_cluster=="L2.3" & cdata$smoothed_k9_1663=="L2", "custom_cluster"] = "L2"
cdata[cdata$custom_cluster=="L2.3" & cdata$smoothed_k9_1663!="L2", "custom_cluster"] = "L3"

cdata$custom_cluster = factor(cdata$custom_cluster, levels=c("Micro.Vasc","Astro.L1","Astro.Nrn","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("Micro.Vasc","Astro.L1","Astro.Nrn","L2","L3","L4","Inhb","L5","L6","WM"))

cat("\nCustom clusters to spe...\n")
spe = spe[,rownames(cdata)]
stopifnot(identical(rownames(cdata), colnames(spe)))
spe$custom_cluster = cdata$custom_cluster
table(spe$custom_cluster, useNA="ifany")

#pseudobulk raw counts
spe
cat("\nPseudobulk spe by: sample_id, custom_cluster\n")
spe_pseudo <- aggregateAcrossCells(spe, ids=colData(spe)[,c("sample_id","custom_cluster")], statistics="sum", store.number="nspots", use.assay.type="counts")
dim(spe_pseudo)

#quick save checkpoints
save(spe_pseudo, file="processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster.Rdata")


#remove repeated colData column for sample_id and cluster
g1 = grep("custom", colnames(colData(spe_pseudo)))
if(length(g1)>1) colData(spe_pseudo)[,g1[[2]]] <- NULL
g2 = grep("sample_id", colnames(colData(spe_pseudo)))
if(length(g2)>1) colData(spe_pseudo)[,g2[[2]]] <- NULL
#change name of ncells to nspots
colnames(colData(spe_pseudo))[grep("ncells", colnames(colData(spe_pseudo)))] = "nspots"

#remove reduced dims
if(length(reducedDimNames(spe_pseudo))>0) {
	for(i in reducedDimNames(spe_pseudo)) {
		reducedDim(spe_pseudo, i) <- NULL
	}
}
imgData(spe_pseudo) <- NULL
spatialCoords(spe_pseudo) <-NULL

#keep only sample level coldata
colData(spe_pseudo) = colData(spe_pseudo)[,c("sample_id","brnum","age","sex","condition","PMI","RIN",
	"slide","array","MBv_sample","seq","round",
	"custom_cluster","nspots")]
colData(spe_pseudo)$condition = factor(spe_pseudo$condition, levels=c("NTC","MDD","BPD"))
spe_pseudo

Sys.time()
save(spe_pseudo, file="processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster.Rdata")

#update spe tracker
#write(c(paste("******* Created pseudobulked spe on",format(Sys.time()),"EST"),
#        "******* Old file location: processed-data/04_feature_selection/spe_n119_postQC-conservative_norm_",
#        "******* New file location: processed-data/06_pseudobulk/Seurat/spe_n119_conservative_pseudo_sample-seurat-pc20-filtered.Rdata",
#        "******* Source code: code/06_pseudobulk/Seurat/01_create_pseudobulk.r",
#        "*******","*******","*******"), "spe_tracker_current.txt", append=TRUE)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
