setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(Seurat)
	library(scuttle)
})
set.seed(123)

#load seu_con MBv-filtered with metadata
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
mdata= seu_con@meta.data
varfeat = VariableFeatures(seu_con)
rm(seu_con)

#load seu_con with original genes
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_processed-SCT.Rdata")

##transfer cluster labels
#seu_con$seurat_low.res <- factor(as.character(mdata$seurat_low.res), 
#	levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
#	labels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
#seu_con$azimuth_broad <- factor(as.character(seu_con$azimuth),
#                                levels=c("Endo","PC","SMC","VLMC",
#                                         "Immune","Micro",
#                                         "Astro",
#                                         "L2/3 IT", "L4 IT",
#                                         "L5 IT", "L5 ET", "L5/6 NP",
#                                         "L6 IT","L6 IT Car3","L6 CT","L6b",
#                                         "Oligo","OPC",
#                                         "Sncg","Pax6","Vip",
#                                         "Lamp5","Lamp5 Lhx6",
#                                         "Sst","Sst Chodl","Pvalb",
#                                         "Chandelier"),
#                                labels=c("Vasc","Vasc","Vasc","Vasc",
#                                         "Micro","Micro",
#                                         "Astro",
#                                         "L2/3","L4",
#                                         "L5","L5","L5",
#                                         "L6","L6","L6","L6",
#                                         "Oligo","Oligo",
#                                         "CGE VIP","CGE VIP","CGE VIP",
#                                         "CGE LAMP5","CGE LAMP5",
#                                         "MGE SST", "MGE SST", "MGE PVALB",
#                                         "Chandelier")
#)
seu_con$azimuth_super.broad <- factor(as.character(seu_con$azimuth),
                                levels=c("Endo","PC","SMC","VLMC",
                                         "Immune","Micro",
                                         "Astro",
                                         "L2/3 IT", "L4 IT",
                                         "L5 IT", "L5 ET", "L5/6 NP",
                                         "L6 IT","L6 IT Car3","L6 CT","L6b",
                                         "Oligo","OPC",
                                         "Sncg","Pax6","Vip",
                                         "Lamp5","Lamp5 Lhx6",
                                         "Sst","Sst Chodl","Pvalb",
                                         "Chandelier"),
                                labels=c("Vasc","Vasc","Vasc","Vasc",
                                         "Micro","Micro",
                                         "Astro",
                                         "ExcN","ExcN",
                                         "ExcN","ExcN","ExcN",
                                         "ExcN","ExcN","ExcN","ExcN",
                                         "Oligo","Oligo",
                                         "InhN","InhN","InhN",
                                         "InhN","InhN",
                                         "InhN", "InhN", "InhN",
                                         "InhN")
)
#to sce
### feature/row data
fdata = seu_con[["RNA"]]@meta.data
fdata$gene_name = rownames(seu_con[["RNA"]])
colnames(fdata)[1] = "gene_id"
rownames(fdata) = fdata$gene_id
fdata$VarFeat_MBv.filtered = fdata$gene_name %in% varfeat
### pull count matrix
mtx = seu_con[["RNA"]]$counts
rownames(mtx) = fdata$gene_id
colnames(mtx) = colnames(seu_con[["RNA"]])

### metadata/coldata
cdata = seu_con@meta.data[,c("Channel","demux_type","assignment",#"anno","subclass","azimuth",
	"individualID","Cohort","Biological_Sex","Age_death","Disorder","azimuth_super.broad")]
	#"seurat_low.res")]

rm(seu_con)
sce_con <- SingleCellExperiment(assays = list(counts = mtx), colData=cdata)
rowData(sce_con) = fdata
sce_con
#table(colData(sce_con)$seurat_low.res, useNA="ifany")
table(colData(sce_con)$azimuth_super.broad, useNA="ifany")

#normalization
format(Sys.time())
cat("\nCompute library factors and normalize counts")
sce_con <- computeLibraryFactors(sce_con)
sce_con <- logNormCounts(sce_con)

#pseudobulk raw counts
#cat("\nPseudobulk sce by: Disorder, seurat_low.res\n")
cat("\nPseudobulk sce by: Disorder, azimuth_super.broad\n")

sce_summ = aggregateAcrossCells(sce_con, ids=colData(sce_con)[,c("Disorder","azimuth_super.broad")],
		#"seurat_low.res")], 
                            statistics=c("mean","prop.detected"),
                            use.assay.type="logcounts")

dim(sce_summ)

#quick save checkpoints
save(sce_summ, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")

#remove repeated colData column for sample_id and cluster
#g1 = grep("seurat", colnames(colData(sce_summ)))
g1 = grep("azimuth", colnames(colData(sce_summ)))
if(length(g1)>1) colData(sce_summ)[,g1[[2]]] <- NULL
g2 = grep("Disorder", colnames(colData(sce_summ)))
if(length(g2)>1) colData(sce_summ)[,g2[[2]]] <- NULL

##change name of ncells to nspots
#colnames(colData(spe_summ))[grep("ncells", colnames(colData(spe_summ)))] = "nspots"

##remove reduced dims
#if(length(reducedDimNames(spe_summ))>0) {
#	for(i in reducedDimNames(spe_summ)) {
#		reducedDim(spe_summ, i) <- NULL
#	}
#}
#imgData(spe_summ) <- NULL
#spatialCoords(spe_summ) <-NULL

#keep only sample level coldata
#colData(sce_summ) = colData(sce_summ)[,c("Disorder","seurat_low.res","ncells")]
colData(sce_summ) = colData(sce_summ)[,c("Disorder","azimuth_super.broad","ncells")]


Sys.time()
save(sce_summ, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")
cat("\nPseudobulk sce saved to: processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_azimuth-super-broad.Rdata")

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
