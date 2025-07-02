setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(SingleCellExperiment)
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

#transfer cluster labels
seu_con$seurat_low.res <- factor(as.character(mdata$seurat_low.res), 
	levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
	labels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))

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
	"individualID","Cohort","Biological_Sex","Age_death","Disorder","seurat_low.res")]

rm(seu_con)
sce_con <- SingleCellExperiment(assays = list(counts = mtx), colData=cdata)
rowData(sce_con) = fdata
sce_con
table(colData(sce_con)$seurat_low.res, useNA="ifany")


#pseudobulk raw counts
cat("\nPseudobulk sce by: individualID, seurat_low.res\n")
sce_pseudo <- aggregateAcrossCells(sce_con, ids=colData(sce_con)[,c("individualID","seurat_low.res")], statistics="sum", use.assay.type="counts")
dim(sce_pseudo)

#remove repeated colData columns for pseudobulking levels
g1 = grep("seurat", colnames(colData(sce_pseudo)))
if(length(g1)>1) colData(sce_pseudo)[,g1[[2]]] <- NULL
g2 = grep("individualID", colnames(colData(sce_pseudo)))
if(length(g2)>1) colData(sce_pseudo)[,g2[[2]]] <- NULL

Sys.time()
save(sce_pseudo, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res.Rdata")
cat("\nPseudobulk sce saved to: processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res.Rdata")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
