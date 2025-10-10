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
dim(spe)
#load clusters
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc20_red-precast-kweight-50-low-res.csv", row.names=1)
#stopifnot(identical(rownames(colData(spe)), rownames(res)))
#res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-MBv-conservative_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)

spe = spe[,rownames(res)]
cat("\nDim spe after filtering out spots without PRECAST embeddings:\n")
dim(spe)
stopifnot(identical(rownames(colData(spe)), rownames(res)))

#transfer label IDs and combine L2/3
spe$seurat_pc20 = factor(res$predicted.id, #levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
	levels=c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))
########### unmerge L2 and L3 for pc20
	#labels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
	#labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))
###########

#spe$seurat_qual.genes_pc30.kweight50 = factor(res4$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
#        labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))
cat("\nTransferred label transfer from qual genes results with PC=20 and k.weights=50 to spe:\n")
table(spe$seurat_pc20, useNA="ifany")
#table(spe$seurat_qual.genes_pc30.kweight50, useNA="ifany")


##remove low UMI cluster
#cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
#stopifnot(identical(rownames(colData(spe)), rownames(cdata)))
#
#spe$lowUMI_cluster = cdata$precast_k9_1663_f=="low UMI"
#
#spe = spe[,spe$lowUMI_cluster==F]
#cat("\n\nRemoved low UMI cluster spots...\n")
#dim(spe)
#table(spe$seurat_qual.genes_pc30.kweight50, useNA="ifany")


#realize counts
cat("\nRealize logcounts matrix...\n")
format(Sys.time())
regular_matrix_counts <- as.matrix(assays(spe)[["logcounts"]])
sparse_matrix_counts <- as(regular_matrix_counts, "dgCMatrix")
rownames(sparse_matrix_counts) <- rownames(logcounts(spe))
colnames(sparse_matrix_counts) <- colnames(logcounts(spe))
#logcounts(spe) <- sparse_matrix_counts

#i keep getting a weird error when trying with SpatialExperiment that fdoesn't happen with the Seurat SingleCellExperiment so try to convert to sce first
sce <- SingleCellExperiment(assays = list(logcounts = sparse_matrix_counts), colData=colData(spe))
rowData(sce) = rowData(spe)

#pseudobulk logcounts
cat("\nPseudobulk spe by: condition, sex, seurat_pc20\n")
spe_summ = aggregateAcrossCells(sce, ids=colData(sce)[,c("condition","sex","seurat_pc20")], 
                            statistics=c("mean","prop.detected"),
                            use.assay.type="logcounts")
dim(spe_summ)

#quick save checkpoints
save(spe_summ, file="processed-data/06_pseudobulk/Seurat/spe_n119_conservative_pseudo-dotplot_dx-sex-seurat-pc20.Rdata")

#remove repeated colData column for sample_id and cluster
g1 = grep("seurat", colnames(colData(spe_summ)))
if(length(g1)>1) colData(spe_summ)[,g1[[2]]] <- NULL
g2 = grep("condition", colnames(colData(spe_summ)))
if(length(g2)>1) colData(spe_summ)[,g2[[2]]] <- NULL
g3 = grep("sex", colnames(colData(spe_summ)))
if(length(g3)>1) colData(spe_summ)[,g3[[2]]] <- NULL

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
	"sex","condition",
	#"PMI","RIN","slide","array","MBv_sample","seq","round",
	#"seurat_qual.genes_pc20.kweight50",
	"seurat_pc20",
	"nspots")]
colData(spe_summ)$condition = factor(spe_summ$condition, levels=c("NTC","MDD","BPD"))
colData(spe_summ)$seurat_label = colData(spe_summ)$seurat_pc20


Sys.time()
save(spe_summ, file="processed-data/06_pseudobulk/Seurat/spe_n119_conservative_pseudo-dotplot_dx-sex-seurat-pc20.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/Seurat/spe_n119_conservative_pseudo-dotplot_dx-sex-seurat-pc20.Rdata")

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
