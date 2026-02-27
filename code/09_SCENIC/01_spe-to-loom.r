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

#for rowdata want to import pseudobulk edgeR filterByExpr results
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
cat("\nPseudobulk spe for gene filtering...\n")
dim(spe_pseudo) # 21080   690
rowData(spe_pseudo)$high_expr_group_sample_id2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- edgeR::filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)

rdata = as.data.frame(rowData(spe_pseudo)[,c("gene_id","gene_name","gene_type")])
rdata$DE_input = rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T
table(rdata$DE_input)
#FALSE  TRUE 
# 7236 13844


## single sample
#cat("\n\nSingle sample test subset...\n")
#(ssample = "V13B23-329_A1")
#(fn1 = "processed-data/09_SCENIC/single-sample_V13B23-329-A1_17300-genes.loom")
#spe_sub1 = spe[rownames(rdata),spe$sample_id==ssample]
##spe_sub1
#t1 = rowSums(counts(spe_sub1)>0)
##table(t1>3)
#spe_sub1 = spe_sub1[t1>3,]
#dim(spe_sub1)
##identical(rdata[rownames(spe_sub1),"gene_id"], rownames(spe_sub1))
##identical(cdata[colnames(spe_sub1),"key"], colnames(spe_sub1))
#
#lfile = create(fn1, counts(spe_sub1), gene.attrs=as.list(rdata[rownames(spe_sub1),]),
#       cell.attrs=as.list(cdata[colnames(spe_sub1),c("key","sample_id","brnum","condition","sex","smoothed_k9_1663","nUMI","nGene")]),
#       calc.numi = F, chunk.dims=c(100,500))
##lfile
##lfile[["matrix"]]
##lfile[["col_attrs"]]
##lfile[["row_attrs"]]
## for whatever reason the gene attributes didn't get added so add them now 
#stopifnot(identical(rownames(spe_sub1), lfile[["row_attrs/Gene"]][]))
#lfile$add.row.attribute(list("gene_id" = rownames(spe_sub1),
#                             "gene_name" = rdata[rownames(spe_sub1),"gene_name"],
#                             "gene_type" = rdata[rownames(spe_sub1),"gene_type"],
#                             "DE_input" = rdata[rownames(spe_sub1),"DE_input"]), 
#                        overwrite = TRUE)
#
#lfile$close_all()

## single slide
#cat("\n\nSingle slide test subset...\n")
#(sslide = "V13B23-329")
#(fn2 = "processed-data/09_SCENIC/single-slide_V13B23-329_17300-genes.loom")
#spe_sub2 = spe[rownames(rdata), spe$slide==sslide]
#spe_sub2 = spe_sub2[t1>3,]
#dim(spe_sub2)
#
#lfile = create(fn2, counts(spe_sub2), gene.attrs=as.list(rdata[rownames(spe_sub2),]),
#               cell.attrs=as.list(cdata[colnames(spe_sub2),c("key","sample_id","brnum","condition","sex","smoothed_k9_1663","nUMI","nGene")]),
#               calc.numi = F, chunk.dims=c(100,500))
#stopifnot(identical(rownames(spe_sub2), lfile[["row_attrs/Gene"]][]))
#lfile$add.row.attribute(list("gene_id" = rownames(spe_sub2),
#                             "gene_name" = rdata[rownames(spe_sub2),"gene_name"],
#                             "gene_type" = rdata[rownames(spe_sub2),"gene_type"],
#                             "DE_input" = rdata[rownames(spe_sub2),"DE_input"]), 
#                        overwrite = TRUE)
#lfile$close_all()


# all data
cat("\n\nFull dataset...\n")
spe = spe[rownames(rdata),]
dim(spe)
(fn3 = "processed-data/09_SCENIC/spe-n119_21080-genes.loom")

lfile = create(fn3, counts(spe), gene.attrs=as.list(rdata),
               cell.attrs=as.list(cdata[colnames(spe),c("key","sample_id","brnum","condition","sex","smoothed_k9_1663","nUMI","nGene")]),
               calc.numi = F, chunk.dims=c(100,500))
stopifnot(identical(rownames(spe), lfile[["row_attrs/Gene"]][]))
lfile$add.row.attribute(list("gene_id" = rownames(spe),
                             "gene_name" = rdata[rownames(spe),"gene_name"],
                             "gene_type" = rdata[rownames(spe),"gene_type"],
                             "DE_input" = rdata[rownames(spe),"DE_input"]), 
                        overwrite = TRUE)
lfile$close_all()


## Reproducibility information
cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
