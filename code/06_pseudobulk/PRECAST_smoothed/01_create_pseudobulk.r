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
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_")

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

#drop low UMI spots that couldn't be saved
cat("\nHow many spots are dropped after smoothing clusters:\n")
#table(cdata[cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),"smoothed_k9_1663_f"])
table(cdata[cdata$smoothed_k7_1626=="low UMI", "smoothed_k7_1626"])

#cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),]
cdata2 = cdata[cdata$smoothed_k7_1626!="low UMI",]
spe = spe[,rownames(cdata2)]
#spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM"),
#	labels=c("L1","L2","L3.4","L5","L6","WM"))
#cat("\nTransferred smoothed PRECAST k=9 n1663 to spe:\n")
#table(spe$smoothed_k9_1663, useNA="ifany")

spe$smoothed_k7_1626 = factor(cdata2$smoothed_k7_1626, levels=c("L1","L2","L3/4","L5","L6","WM"),
                                labels=c("L1","L2","L3.4","L5","L6","WM"))
cat("\nTransferred smoothed PRECAST k=7 n1626 to conservative spe:\n")
table(spe$smoothed_k7_1626, useNA="ifany")

# drop spots but keep low UMI
#cat("\nHow many spots are dropped after smoothing clusters (keep low UMI):\n")
#table(cdata[cdata$smoothed_k9_1663_f %in% c("GABA","Vasc"),"smoothed_k9_1663_f"])
#cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("GABA","Vasc"),]
#spe = spe[,rownames(cdata2)]
#spe$smoothed_k9_1663 = factor(cdata2$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
#        labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))
#cat("\nTransferred smoothed PRECAST k=9 n1663 to spe:\n")
#table(spe$smoothed_k9_1663, useNA="ifany")

## conservative keep low UMI
#spe$smoothed_k7_1626 = factor(cdata$smoothed_k7_1626, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
#                                labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))
#cat("\nTransferred smoothed PRECAST k=7 n1626 to conservative spe:\n")
#table(spe$smoothed_k7_1626, useNA="ifany")

#pseudobulk raw counts
spe
cat("\nPseudobulk spe by: sample_id, smoothed_k7_1626\n")
spe_pseudo <- aggregateAcrossCells(spe, ids=colData(spe)[,c("sample_id","smoothed_k7_1626")], 
	statistics="sum", store.number="nspots", use.assay.type="counts")
dim(spe_pseudo)

#remove repeated colData column for sample_id and cluster
g1 = grep("smoothed_k7_1626", colnames(colData(spe_pseudo)))
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
colData(spe_pseudo) = colData(spe_pseudo)[,c("sample_id","brnum","age","sex","condition","PMI","RIN","slide","array",
	"MBv_sample","seq","round","smoothed_k7_1626","nspots")]
colData(spe_pseudo)$condition = factor(spe_pseudo$condition, levels=c("NTC","MDD","BPD"))
spe_pseudo

Sys.time()
save(spe_pseudo, file="processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7.Rdata")

update spe tracker
write(c(paste("******* Created pseudobulked spe on",format(Sys.time()),"EST"),
        "******* Old file location: processed-data/04_feature_selection/spe_n119_postQC-conservative_norm_",
        "******* New file location: processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7.Rdata",
        "******* Source code: code/06_pseudobulk/PRECAST_smoothed/01_create_pseudobulk.r",
        "*******","*******","*******"), "spe_tracker_current.txt", append=TRUE)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
