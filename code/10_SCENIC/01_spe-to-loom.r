#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(loomR)
})
set.seed(123)
setAutoBlockSize(1e9)

#load spe
cat("\nFull spe...\n")
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
dim(spe) # 28965 535248

#load clusters
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
stopifnot(identical(rownames(colData(spe)), rownames(cdata)))
cdata$smoothed_k9_1663 = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"),
                                labels=c("L1","L2","L3.4","L5","L6","WM","low.UMI","Vasc","GABA"))
cdata$key = colData(spe)$key
cdata$nUMI = colData(spe)$sum_umi
cdata$nGene = colData(spe)$sum_gene

# add seurat label
res = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
res = res[rownames(cdata),]
stopifnot(identical(rownames(res), rownames(cdata)))
cdata$seurat_label = factor(res$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
	labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))

#for rowdata want to import pseudobulk edgeR filterByExpr results
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata")
cat("\nPseudobulk spe for gene filtering...\n")
dim(spe_pseudo) # 21080   690
rowData(spe_pseudo)$high_expr_group_sample_id2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)

rdata = as.data.frame(rowData(spe_pseudo)[,c("gene_id","gene_name","gene_type")])
rdata$DE_input = rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T
table(rdata$DE_input)

cat("\nFiltered gene set to DE input...\n")
rdata = rdata[rdata$DE_input,]
dim(rdata)

spe = spe[rownames(rdata),]


# remove low UMI cluster
cat("\nRemoving low UMI cluster...\n")
cdata = cdata[cdata$smoothed_k9_1663!="low.UMI",]
cdata$smoothed_k9_1663 = droplevels(cdata$smoothed_k9_1663)

spe = spe[,rownames(cdata)]
dim(spe)

# DE input version
rowData(spe)[c("ENSG00000187522"),"gene_name"] = "MSTANTD7"

# make new rdata df with updated info and swap rownames to gene_name
rdata2 = rowData(spe)[,c("gene_id","gene_name","gene_type")]
rdata2$DE_input = rdata[rownames(rdata2),"DE_input"]
rownames(rdata2) = rdata2$gene_name
cat("\nNumber of genes after revising to unique gene names:\n")
nrow(rdata2)


# remove 10 npas4 samples
cat("\n\nFull  dataset:", dim(spe), "\n")
npas4.outliers = c("Br5666","Br5594","Br5599","Br5448","Br6316","Br6021","Br5993","Br6192","Br5454","Br8073")
spe = spe[,!spe$brnum %in% npas4.outliers]
cat("\nRemoved 10 Npas4 outlier samples:", dim(spe),"\n")
(fn3 = paste0("processed-data/10_SCENIC/expr_loom/spe-n109_", nrow(rdata2), "-genes_no-lowUMI.loom"))

#cat("\n\nFull dataset...\n")
#(fn3 = paste0("processed-data/10_SCENIC/expr_loom/spe-n119_", nrow(rdata2), "-genes_no-lowUMI_logcounts.loom"))

# extract counts matrix and change rownames
mtx3 = counts(spe)
#mtx3 = logcounts(spe)
stopifnot(identical(rownames(mtx3), rdata2$gene_id))
rownames(mtx3) = rdata2$gene_name
dim(mtx3)

lfile = create(fn3, mtx3, #gene.attrs=as.list(rdata),
               cell.attrs=as.list(cdata[colnames(spe),c("key","sample_id","brnum","condition","sex","smoothed_k9_1663","seurat_label","nUMI","nGene")]),
               calc.numi = F, chunk.dims=c(100,500))
stopifnot(identical(rownames(rdata2), lfile[["row_attrs/Gene"]][]))
lfile$add.row.attribute(list("gene_id" = rdata2$gene_id,
                             "gene_name" = rownames(rdata2),
                             "gene_type" = rdata2$gene_type,
                             "DE_input" = rdata2$DE_input), 
                        overwrite = TRUE)
lfile$close_all()


## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
