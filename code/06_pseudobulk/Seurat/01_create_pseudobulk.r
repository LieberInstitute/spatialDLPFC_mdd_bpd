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
res2 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc20_red-precast-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(res2)))
#res4 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
#stopifnot(identical(rownames(colData(spe)), rownames(res4)))

#transfer label IDs and combine L2/3
spe$seurat_qual.genes_pc20.kweight50 = factor(res2$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
########### unmerge L2 and L3 for pc20
	labels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
	#labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))
###########

#spe$seurat_qual.genes_pc30.kweight50 = factor(res4$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
#        labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))
cat("\nTransferred label transfer from qual genes results with PC=20 and k.weights=50 to spe:\n")
table(spe$seurat_qual.genes_pc20.kweight50, useNA="ifany")
#table(spe$seurat_qual.genes_pc30.kweight50, useNA="ifany")

#pseudobulk raw counts
spe
cat("\nPseudobulk spe by: sample_id, seurat_qual.genes_pc20.kweight50\n")
spe_pseudo <- aggregateAcrossCells(spe, ids=colData(spe)[,c("sample_id","seurat_qual.genes_pc20.kweight50")], statistics="sum", store.number="nspots", use.assay.type="counts")
dim(spe_pseudo)

#quick save checkpoints
save(spe_pseudo, file="processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20.Rdata")

#remove repeated colData column for sample_id and cluster
g1 = grep("seurat", colnames(colData(spe_pseudo)))
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
	"seurat_qual.genes_pc20.kweight50",
	#"seurat_qual.genes_pc30.kweight50",
	"nspots")]
colData(spe_pseudo)$condition = factor(spe_pseudo$condition, levels=c("NTC","MDD","BPD"))
spe_pseudo

Sys.time()
save(spe_pseudo, file="processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20.Rdata")

#update spe tracker
write(c(paste("******* Created pseudobulked spe on",format(Sys.time()),"EST"),
        "******* Old file location: processed-data/04_feature_selection/spe_n119_postQC_norm_",
        "******* New file location: processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20.Rdata",
        "******* Source code: code/06_pseudobulk/Seurat/01_create_pseudobulk.r",
        "*******","*******","*******"), "spe_tracker_current.txt", append=TRUE)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
