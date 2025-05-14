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
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata = cdata[rownames(colData(spe)),]
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots that couldn't be saved
cat("\nHow many spots are dropped with low UMI cluster:\n")
table(cdata$precast_k9_1663_f=="low UMI")

cdata2 = cdata[cdata$precast_k9_1663_f!="low UMI",]
spe = spe[,rownames(cdata2)]
spe$precast_k9_1663 = factor(cdata2$precast_k9_1663_f, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))
cat("\nTransferred combined cluster to spe:\n")
table(spe$precast_k9_1663, useNA="ifany")

#pseudobulk raw counts
#assays(spe, "logcounts") <- NULL
spe
cat("\nPseudobulk spe by: sample_id, precast_k9_1663 (low UMI removed)\n")
Sys.time()
spe_pseudo <- aggregateAcrossCells(spe, ids=colData(spe)[,c("sample_id","precast_k9_1663")], statistics="sum", store.number="nspots", use.assay.type="counts")
colnames(colData(spe_pseudo))[grep("ncells", colnames(colData(spe_pseudo)))] = "nspots"
#remove reduced dims
#reducedDim(spe_pseudo, "PCA_1663") <- NULL
#reducedDim(spe_pseudo, "PCA_1104") <- NULL
#reducedDim(spe_pseudo, "PCA_1079") <- NULL
#imgData(spe_pseudo) <- NULL
spatialCoords(spe_pseudo) <-NULL

#tmp save
save(spe_pseudo, file="processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9.Rdata")

#remove empty coldata
na.sums = colSums(is.na(colData(spe_pseudo)))
names(na.sums) = colnames(colData(spe_pseudo))
keep.cols = na.sums[na.sums!=dim(spe_pseudo)[[2]]]
colData(spe_pseudo) = colData(spe_pseudo)[,names(keep.cols)]
colData(spe_pseudo)$condition = factor(spe_pseudo$condition, levels=c("NTC","MDD","BPD"))

spe_pseudo

Sys.time()
save(spe_pseudo, file="processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9.Rdata")

#update spe tracker
#write(c(paste("******* Created pseudobulked spe on",format(Sys.time()),"EST"),
#        "******* Old file location: processed-data/04_feature_selection/spe_n120_postQC_norm_",
#        "******* New file location: processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus.Rdata",
#        "******* Source code: code/06_pseudobulk/01_create_pseudobulk.r",
#        "*******","*******","*******"), "spe_tracker_current.txt", append=TRUE)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
