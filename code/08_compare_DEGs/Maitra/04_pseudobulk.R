setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(SingleCellExperiment)
	library(scuttle)
})
set.seed(123)

# load in full dataset
seu_all = readRDS("processed-data/00_external_snRNAseq/seurat_Maitra2023_qc-filtered.rda")
counts = seu_all[["RNA"]]@layers$counts
colnames(counts) = colnames(seu_all)
rownames(counts) = rownames(seu_all)

rdata = seu_all[["RNA"]]@meta.data
rownames(rdata) <- rdata$gene_id

# load in post-transfer processed objects
load("processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_male_MBv-filtered_processed-SCT_post.Rdata")
load("processed-data/08_compare_DEGs/Maitra/seurat_Maitra2023_female_MBv-filtered_processed-SCT_post.Rdata")

cdata = rbind(seu_m@meta.data[,c(4:34,39)], seu_f@meta.data[,c(4:34,39)])

# remove individuals flagged in sn QC
cat("\nRemove donors flagged in sn QC...\n")
dim(cdata)
cdata = cdata[-grep("M24", cdata$Sample),]
dim(cdata)
cdata = cdata[cdata$Sample!="F20",]
dim(cdata)
cdata =	cdata[cdata$Sample!="F18",]
dim(cdata)

# subset counts
counts = counts[,rownames(cdata)]

# build SCE and pseudobulk
sce <- SingleCellExperiment(assays = list(counts = counts), colData=cdata)
rowData(sce) = rdata

dim(sce)
table(colData(sce)[,c("Sex","Condition")])

cat("\nPseudobulk...\n")
sce_pseudo <- aggregateAcrossCells(sce, ids=colData(sce)[,c("Sample","seurat_low.res")],
	statistics="sum", use.assay.type="counts")
dim(sce_pseudo)

# remove repeated colData columns for pseudobulking levels
g1 = grep("seurat", colnames(colData(sce_pseudo)))
if(length(g1)>1) colData(sce_pseudo)[,g1[[2]]] <- NULL
g2 = grep("Sample", colnames(colData(sce_pseudo)))
if(length(g2)>1) colData(sce_pseudo)[,g2[[2]]] <- NULL

# rename factor levels
sce_pseudo$condition = factor(sce_pseudo$Condition, levels=c("Control","Case"), labels=c("NTC","MDD"))
sce_pseudo$sex = factor(sce_pseudo$Sex, levels=c("Female","Male"), labels=c("F","M"))
sce_pseudo$seurat_low.res = factor(sce_pseudo$seurat_low.res, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"))

format(Sys.time())
save(sce_pseudo, file="processed-data/08_compare_DEGs/Maitra/sce_Maitra_pseudo_sample-low-res.Rdata")
cat("\nPseudobulk sce saved to: processed-data/08_compare_DEGs/Maitra/sce_Maitra_pseudo_sample-low-res.Rdata\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
