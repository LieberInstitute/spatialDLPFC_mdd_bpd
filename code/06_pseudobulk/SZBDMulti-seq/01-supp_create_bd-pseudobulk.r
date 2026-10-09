setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(SingleCellExperiment)
	library(scuttle)
})
set.seed(123)

#load seu_bd with original genes
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar.Rdata")
#seu_bd = readRDS("processed-data/00_external_snRNAseq/seurat_SZBDMultiseq.rda")
# this is actually every donor but didn't want to rename all the objects in this script
seu_bd

#head(seu_bd@meta.data)
#write.csv(seu_bd@meta.data, "processed-data/00_external_snRNAseq/check_metadata.csv")

## retrieve channel
#pullfirst = unlist(lapply(strsplit(as.character(seu_bd$Channel), "_"), function(x) x[[1]]))
#seu_bd@meta.data$Channel = pullfirst

## make azimuth broad
#seu_bd$azimuth_broad <- factor(as.character(seu_bd$azimuth),
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
#					 "Inhb","Inhb","Inhb",
#					 "Inhb","Inhb",
#					 "Inhb","Inhb","Inhb",
#					 "Inhb"))

##transfer results
res1 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(colnames(seu_bd), rownames(res1)))
seu_bd@meta.data$seurat_low.res = res1$predicted.id

#to sce
### feature/row data
fdata = seu_bd[["RNA"]]@meta.data
fdata$gene_name = rownames(seu_bd[["RNA"]])
colnames(fdata)[1] = "gene_id"
rownames(fdata) = fdata$gene_id

### metadata/coldata
cdata = seu_bd@meta.data[,c("Channel","demux_type","assignment",#"anno","subclass","azimuth",
	"individualID","Cohort","Biological_Sex","Age_death","Disorder",
#	"azimuth_broad")]
	"seurat_low.res")]


# pull count matrix
mtx = seu_bd[["RNA"]]$counts
rownames(mtx) = fdata$gene_id
colnames(mtx) = colnames(seu_bd[["RNA"]])

rm(seu_bd)
sce_bd <- SingleCellExperiment(assays = list(counts = mtx), colData=cdata)
rowData(sce_bd) = fdata


# load con too
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")
mdata= seu_con@meta.data
varfeat = VariableFeatures(seu_con)
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control.Rdata")
fdata = seu_con[["RNA"]]@meta.data
fdata$gene_name = rownames(seu_con[["RNA"]])
colnames(fdata)[1] = "gene_id"
rownames(fdata) = fdata$gene_id
fdata$VarFeat_MBv.filtered = fdata$gene_name %in% varfeat

stopifnot(identical(colnames(seu_con), rownames(mdata)))
seu_con$seurat_low.res <- as.character(mdata$seurat_low.res)
cdata = seu_con@meta.data[,c("Channel","demux_type","assignment",#"anno","subclass","azimuth",
        "individualID","Cohort","Biological_Sex","Age_death","Disorder",
#       "azimuth_broad")]
        "seurat_low.res")]

mtx = seu_con[["RNA"]]$counts
rownames(mtx) = fdata$gene_id
colnames(mtx) = colnames(seu_con[["RNA"]])

rm(seu_con)
sce_con <- SingleCellExperiment(assays = list(counts = mtx), colData=cdata)
rowData(sce_con) = fdata

sce_both = cbind(sce_con, sce_bd)
sce_both$Channel = unlist(lapply(strsplit(as.character(sce_both$Channel), "_"), function(x) x[[1]]))

#pseudobulk raw counts
cat("\nPseudobulk sce by: individualID, seurat_low.res, channel\n")
sce_pseudo <- aggregateAcrossCells(sce_both, ids=colData(sce_both)[,c("individualID",
	"seurat_low.res")],#"Channel")],
#	"azimuth_broad")],#,"Channel")],
	statistics="sum", use.assay.type="counts")
dim(sce_pseudo)

#remove repeated colData columns for pseudobulking levels
g1 = grep("seurat", colnames(colData(sce_pseudo)))
#g1 = grep("azimuth", colnames(colData(sce_pseudo)))
if(length(g1)>1) colData(sce_pseudo)[,g1[[2]]] <- NULL
g2 = grep("individualID", colnames(colData(sce_pseudo)))
if(length(g2)>1) colData(sce_pseudo)[,g2[[2]]] <- NULL
#g3 = grep("Channel", colnames(colData(sce_pseudo)))
#if(length(g3)>1) colData(sce_pseudo)[,g3[[2]]] <- NULL

Sys.time()
save(sce_pseudo, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_pseudo_indivID-low-res.Rdata")
cat("\nPseudobulk sce saved to: processed-data/06_pseudobulk/SZBDMulti-seq/sce_pseudo_indivID-low-res.Rdata")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

