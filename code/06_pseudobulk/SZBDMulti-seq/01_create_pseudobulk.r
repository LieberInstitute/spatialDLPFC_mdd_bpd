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
umap = seu_con@reductions$umap@cell.embeddings
colnames(umap) = c("UMAP1","UMAP2")
rm(seu_con)

#load seu_con with original genes
load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_processed-SCT.Rdata")

#transfer cluster labels
seu_con$seurat_low.res <- factor(as.character(mdata$seurat_low.res), 
	levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
	labels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))
seu_con$azimuth_broad <- factor(as.character(seu_con$azimuth),
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
                                         "L2/3","L4",
                                         "L5","L5","L5",
                                         "L6","L6","L6","L6",
                                         "Oligo","Oligo",
                                         "CGE CNR1","CGE CNR1","CGE CNR1",
                                         "CGE LAMP5","CGE LAMP5",
                                         "MGE SST", "MGE SST", "MGE PV",
                                         "MGE PV")
)

seu_con = seu_con[,seu_con$azimuth_broad %in% c("CGE CNR1","CGE LAMP5","MGE SST","MGE PV")]

cat("\nFiltering to inhibitory neurons...\n")
dim(seu_con)
seu_con$inhn_type = droplevels(seu_con$azimuth_broad)
seu_con$azimuth_broad <- NULL

umap = umap[colnames(seu_con),]

#to sce
### feature/row data
fdata = seu_con[["RNA"]]@meta.data
fdata$gene_name = rownames(seu_con[["RNA"]])
colnames(fdata)[1] = "gene_id"
rownames(fdata) = fdata$gene_id
#fdata$VarFeat_MBv.filtered = fdata$gene_name %in% varfeat

### pull count matrix
mtx = seu_con[["RNA"]]$counts
rownames(mtx) = fdata$gene_id
colnames(mtx) = colnames(seu_con[["RNA"]])

### metadata/coldata
cdata = seu_con@meta.data[,c("Channel","demux_type","assignment","azimuth",#"anno","subclass","azimuth",
	"individualID","Cohort","Biological_Sex","Age_death","Disorder",#"azimuth_broad",
	"inhn_type")]

rm(seu_con)
sce_con <- SingleCellExperiment(assays = list(counts = mtx), colData=cdata)
rowData(sce_con) = fdata
reducedDim(sce_con, "UMAP") <- umap
sce_con
#table(colData(sce_con)$seurat_low.res, useNA="ifany")
table(colData(sce_con)$inhn_type, useNA="ifany")

#normalize
sce_con <- computeLibraryFactors(sce_con)
sce_con <- logNormCounts(sce_con)

save(sce_con, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_InhN-only.Rdata")
cat("\nSaved full size sce object to: processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_InhN-only.Rdata\n")
#stop("Re-run only to save full size object")

#pseudobulk raw counts
cat("\nPseudobulk sce by: individualID, seurat_low.res\n")
sce_con$azimuth <- NULL
sce_pseudo <- aggregateAcrossCells(sce_con, ids=colData(sce_con)[,c("individualID","inhn_type")],
	#"seurat_low.res")],
	statistics="sum", use.assay.type="counts")
dim(sce_pseudo)

#remove repeated colData columns for pseudobulking levels
#g1 = grep("seurat", colnames(colData(sce_pseudo)))
g1 = grep("inhn_type", colnames(colData(sce_pseudo)))
if(length(g1)>1) colData(sce_pseudo)[,g1[[2]]] <- NULL
g2 = grep("individualID", colnames(colData(sce_pseudo)))
if(length(g2)>1) colData(sce_pseudo)[,g2[[2]]] <- NULL

Sys.time()
save(sce_pseudo, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-InhN_indivID-inhn-type.Rdata")
cat("\nPseudobulk sce saved to: processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-InhN_indivID-inhn-type.Rdata")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
