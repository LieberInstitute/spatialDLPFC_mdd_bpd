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

cat("\nFiltered gene set to DE input...\n")
rdata = rdata[rdata$DE_input,]
dim(rdata)

spe = spe[rownames(rdata),]

#need to switch rownames to gene_name to match with SCENIC_aux files
##check for repeated gene names (multiple gene ids)
#t1 = table(rdata$gene_name)
#t1[t1==2]
## there are 8 gene_names that are duplicated
## for the DE input there are only 2
#avg.expr = read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1663-k9_filtered-genes_avg-logcounts.csv", row.names=1)
#fixme = avg.expr[avg.expr$gene_name %in% names(t1)[t1==2],]
#fixme = fixme[order(fixme$gene_name),]

# i googled the gene ids to decide when to name differently and when to combine
# name differently
#rowData(spe)[c("ENSG00000187522","ENSG00000271858","ENSG00000261186","ENSG00000285053","ENSG00000269226"),"gene_name"] = c("MSTANTD7","LOC101928965","LINC03100","GGPS1-TBCE","TMSB15C")

# DE input version
rowData(spe)[c("ENSG00000187522","ENSG00000271858"),"gene_name"] = c("MSTANTD7","LOC101928965")

#rowData(spe)["ENSG00000285053","gene_type"] = "lncRNA"
# combine
#counts(spe)["ENSG00000271147",] = colSums(counts(spe)[c("ENSG00000271147","ENSG00000286237"),])
#counts(spe)["ENSG00000188626",] = colSums(counts(spe)[c("ENSG00000188626","ENSG00000261480"),])
#counts(spe)["ENSG00000015479",] = colSums(counts(spe)[c("ENSG00000015479","ENSG00000280987"),])
# remove redundant combinations
#spe = spe[setdiff(rownames(spe),c("ENSG00000286237","ENSG00000261480","ENSG00000280987")),]

# make new rdata df with updated info and swap rownames to gene_name
rdata2 = rowData(spe)[,c("gene_id","gene_name","gene_type")]
rdata2$DE_input = rdata[rownames(rdata2),"DE_input"]
rownames(rdata2) = rdata2$gene_name
cat("\nNumber of genes after revising to unique gene names:\n")
nrow(rdata2)

## single sample
#cat("\n\nSingle sample test subset...\n")
#(ssample = "V13B23-329_A1")
##(fn1 = "processed-data/09_SCENIC/single-sample_V13B23-329-A1_17300-genes.loom")
#spe_sub1 = spe[,spe$sample_id==ssample]
#
## extract counts matrix and change rownames
#mtx1 = counts(spe_sub1)
#rownames(mtx1) = rdata2$gene_name
#
## reduce to genes expressed in at least three spots
#t1 = rowSums(mtx1>0)
#mtx1 = mtx1[t1>3,]
#dim(mtx1)
#
## update rdata for sample and filename
#(fn1 = paste0("processed-data/09_SCENIC/single-sample_V13B23-329-A1_", nrow(mtx1), "-genes.loom"))
#rdata2_ss = rdata2[rownames(mtx1),]
#
#lfile = create(fn1, mtx1, #gene.attrs=as.list(rdata2[rownames(spe_sub1),]),
#       cell.attrs=as.list(cdata[colnames(spe_sub1),c("key","sample_id","brnum","condition","sex","smoothed_k9_1663","nUMI","nGene")]),
#       calc.numi = F, chunk.dims=c(100,500))
##lfile
##lfile[["matrix"]]
##lfile[["col_attrs"]]
##lfile[["row_attrs"]]
## for whatever reason the gene attributes didn't get added so add them now 
#stopifnot(identical(rownames(rdata2_ss), lfile[["row_attrs/Gene"]][]))
#lfile$add.row.attribute(list("gene_id" = rdata2_ss$gene_id,
#                             "gene_name" = rownames(rdata2_ss),
#                             "gene_type" = rdata2_ss$gene_type,
#                             "DE_input" = rdata2_ss$DE_input), 
#                        overwrite = TRUE)
#
#lfile$close_all()
#
## single slide
#cat("\n\nSingle slide test subset...\n")
#(sslide = "V13B23-329")
#(fn2 = paste0("processed-data/09_SCENIC/single-slide_V13B23-329_", nrow(mtx1), "-genes.loom"))
#spe_sub2 = spe[,spe$slide==sslide]
#
## extract counts matrix and change rownames
#mtx2 = counts(spe_sub2)
#rownames(mtx2) = rdata2$gene_name
#
## reduce to genes expressed in at least three spots
#mtx2 = mtx2[t1>3,]
#dim(mtx2)
#
#lfile = create(fn2, mtx2, #gene.attrs=as.list(rdata[rownames(spe_sub2),]),
#               cell.attrs=as.list(cdata[colnames(spe_sub2),c("key","sample_id","brnum","condition","sex","smoothed_k9_1663","nUMI","nGene")]),
#               calc.numi = F, chunk.dims=c(100,500))
#stopifnot(identical(rownames(rdata2_ss), lfile[["row_attrs/Gene"]][]))
#lfile$add.row.attribute(list("gene_id" = rdata2_ss$gene_id,
#                             "gene_name" = rownames(rdata2_ss),
#                             "gene_type" = rdata2_ss$gene_type,
#                             "DE_input" = rdata2_ss$DE_input), 
#                        overwrite = TRUE)
#lfile$close_all()


# all data
cat("\n\nFull dataset...\n")
(fn3 = paste0("processed-data/09_SCENIC/spe-n119_", nrow(rdata2), "-genes.loom"))

# extract counts matrix and change rownames
mtx3 = counts(spe)
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
