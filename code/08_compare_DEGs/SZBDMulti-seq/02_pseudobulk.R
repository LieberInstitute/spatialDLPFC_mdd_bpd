setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(Seurat)
        library(SingleCellExperiment)
        library(scuttle)
})
set.seed(123)

# load in full dataset
seu_all = readRDS("processed-data/00_external_snRNAseq/seurat_SZBDMultiseq_qc-filtered.rda")
counts = seu_all[["RNA"]]@layers$counts
colnames(counts) = colnames(seu_all)
rownames(counts) = rownames(seu_all)

rdata = seu_all[["RNA"]]@meta.data
rownames(rdata) <- rdata$featureid
rownames(counts) <- rdata$featureid

cdata = seu_all@meta.data[,c(1,7:10, 14:22)]
# load in post-transfer processed objects
load("processed-data/08_compare_DEGs/SZBDMulti-seq/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT_post.Rdata")
# "transfer_low.res"
load("processed-data/08_compare_DEGs/SZBDMulti-seq/seurat_SZBDMulti-seq_bipolar_MBv-filtered_processed-SCT_post.Rdata")
# seurat_low.res
# remove BD4
seu_bd = seu_bd[,seu_bd$individualID!="BD4"]


cdata = cdata[c(colnames(seu_con), colnames(seu_bd)),]
cdata$seurat_low.res = NA
cdata[colnames(seu_con), "seurat_low.res"] = as.character(seu_con$transfer_low.res)
cdata[colnames(seu_bd), "seurat_low.res"] = as.character(seu_bd$seurat_low.res)
table(cdata$seurat_low.res)

# subset counts
counts = counts[,rownames(cdata)]

# build SCE and pseudobulk
sce <- SingleCellExperiment(assays = list(counts = counts), colData=cdata)
rowData(sce) = rdata

dim(sce)
table(colData(sce)[,c("Biological_Sex","Disorder")])

cat("\nPseudobulk...\n")
sce_pseudo <- aggregateAcrossCells(sce, ids=colData(sce)[,c("individualID","seurat_low.res")],
        statistics="sum", use.assay.type="counts")
dim(sce_pseudo)

# remove repeated colData columns for pseudobulking levels
g1 = grep("seurat", colnames(colData(sce_pseudo)))
if(length(g1)>1) colData(sce_pseudo)[,g1[[2]]] <- NULL
g2 = grep("individualID", colnames(colData(sce_pseudo)))
if(length(g2)>1) colData(sce_pseudo)[,g2[[2]]] <- NULL

# rename factor levels
sce_pseudo$condition = factor(sce_pseudo$Disorder, levels=c("control","Bipolar Disorder"), labels=c("NTC","BD"))
sce_pseudo$sex = factor(sce_pseudo$Biological_Sex, levels=c("female","male"), labels=c("F","M"))
sce_pseudo$seurat_low.res = factor(sce_pseudo$seurat_low.res, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
        labels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg"))

format(Sys.time())
save(sce_pseudo, file="processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res.Rdata")
cat("\nPseudobulk sce saved to: processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-low-res.Rdata\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
