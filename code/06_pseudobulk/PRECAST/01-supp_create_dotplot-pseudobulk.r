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
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))

##drop low UMI spots that couldn't be saved
#cat("\nHow many spots are dropped after smoothing clusters:\n")
#table(cdata[cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),"smoothed_k9_1663_f"])

#cdata2 = cdata[!cdata$smoothed_k9_1663_f %in% c("low UMI","GABA","Vasc"),]
#spe = spe[,rownames(cdata2)]
spe$precast_k9_1663 = factor(cdata$smoothed_k9_1663_f, levels=c("Vasc","L1","L2","L3/4","GABA","L5","L6","WM","low UMI"),
	labels=c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI"))
cat("\nTransferred PRECAST k=9 n1663 to spe:\n")
table(spe$precast_k9_1663, useNA="ifany")

#i keep getting a weird error when trying with SpatialExperiment that fdoesn't happen with the Seurat SingleCellExperiment so try to convert to sce first
sce <- SingleCellExperiment(assays = list(logcounts = logcounts(spe)), colData=colData(spe))
rowData(sce) = rowData(spe)

#pseudobulk raw counts
cat("\nPseudobulk spe by: precast_k9_1663\n")

spe_summ = scuttle::aggregateAcrossCells(sce, ids=colData(spe)[,c("precast_k9_1663")], 
                            statistics=c("mean","prop.detected"),
                            use.assay.type="logcounts")

dim(spe_summ)

#quick save checkpoints
save(spe_summ, file="processed-data/06_pseudobulk/PRECAST/spe_n119_pseudo-dotplot_precast-n1663-k9.Rdata")

#remove repeated colData column for sample_id and cluster
g1 = grep("precast", colnames(colData(spe_summ)))
if(length(g1)>1) colData(spe_summ)[,g1[[2]]] <- NULL
#g2 = grep("condition", colnames(colData(spe_summ)))
#if(length(g2)>1) colData(spe_summ)[,g2[[2]]] <- NULL
#g3 = grep("sex", colnames(colData(spe_summ)))
#if(length(g3)>1) colData(spe_summ)[,g3[[2]]] <- NULL

#change name of ncells to nspots
colnames(colData(spe_summ))[grep("ncells", colnames(colData(spe_summ)))] = "nspots"

#remove reduced dims
if(length(reducedDimNames(spe_summ))>0) {
	for(i in reducedDimNames(spe_summ)) {
		reducedDim(spe_summ, i) <- NULL
	}
}
#imgData(spe_summ) <- NULL
#spatialCoords(spe_summ) <-NULL

#keep only sample level coldata
colData(spe_summ) = colData(spe_summ)[,c(#"sample_id","brnum","age",
	#"sex","condition",
	#"PMI","RIN","slide","array","MBv_sample","seq","round",
	"precast_k9_1663",
	"nspots")]
#colData(spe_summ)$condition = factor(spe_summ$condition, levels=c("NTC","MDD","BPD"))
#colData(spe_summ)$sex = factor(spe_summ$sex, levels=c("F","M"))

Sys.time()
save(spe_summ, file="processed-data/06_pseudobulk/PRECAST/spe_n119_pseudo-dotplot_precast-n1663-k9.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/PRECAST/spe_n119_pseudo-dotplot_precast-n1663-k9.Rdata")

#update spe tracker
#write(c(paste("******* Created pseudobulked spe on",format(Sys.time()),"EST"),
#        "******* Old file location: processed-data/04_feature_selection/spe_n119_postQC_norm_",
#        "******* New file location: processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20.Rdata",
#        "******* Source code: code/06_pseudobulk/Seurat/01_create_pseudobulk.r",
#        "*******","*******","*******"), "spe_tracker_current.txt", append=TRUE)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
