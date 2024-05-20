library(SpatialExperiment)
library(scater)
library(SpotSweeper)
library(here)

load(here("processed-data","02_build_spe","spe_raw.Rdata"))
colData(spe) = cbind(colData(spe),
do.call(rbind, lapply(strsplit(spe$sample_id, split="_"), function(x) cbind.data.frame("slide"=x[[1]],"position"=x[[2]]))))

colData(spe) = colData(spe)[, c("key","sample_id","slide","position","in_tissue","array_row","array_col","sum_umi","sum_gene","expr_chrM","expr_chrM_ratio")]

spe <- spe[,spe$in_tissue==TRUE]
colData(spe)$lg10.sum = log10(colData(spe)$sum_umi)
colData(spe)$lg10.genes = log10(colData(spe)$sum_gene)

#add outlier info to coldata
colData(spe)$sum_3MAD.outlier_slide = isOutlier(colData(spe)$lg10.sum, batch=colData(spe)$slide, type="lower", nmads=3)

colData(spe)$genes_3MAD.outlier_slide = isOutlier(colData(spe)$lg10.genes, batch=colData(spe)$slide, type="lower", nmads=3)

spe <- localOutliers(spe, metric = "sum_umi", direction = "lower", log = TRUE)
colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))] = "sum_local.outlier"

spe <- localOutliers(spe, metric = "sum_genes", direction = "lower", log = TRUE)
colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))] = "genes_local.outlier"

spe <- localOutliers(spe, metric = "expr_chrM_ratio", direction = "higher", log = FALSE)
colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))] = "chrM.ratio_local.outlier"

save(spe, here("processed-data","03_QC","spe_batch1.Rdata"))
