setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scater)
	library(SpotSweeper)
	library(here)
})

cat(paste(format(Sys.time(), tz="UTC"),"Load raw spe and specify metadata columns...","\n"))
load(here("processed-data","03_QC","spe_demo-filt.Rdata"))
colData(spe)$lg10.sum = log10(colData(spe)$sum_umi)
colData(spe)$lg10.genes = log10(colData(spe)$sum_gene)

cat(paste(format(Sys.time(), tz="UTC"), "Calculate 3MAD outliers...","\n"))
colData(spe)$sum_3MAD.outlier_slide = isOutlier(colData(spe)$lg10.sum, batch=colData(spe)$slide, type="lower", nmads=3)
colData(spe)$sum_3MAD.outlier_sample = isOutlier(colData(spe)$lg10.sum, batch=colData(spe)$sample_id, type="lower", nmads=3)

colData(spe)$genes_3MAD.outlier_slide = isOutlier(colData(spe)$lg10.genes, batch=colData(spe)$slide, type="lower", nmads=3)
colData(spe)$genes_3MAD.outlier_sample = isOutlier(colData(spe)$lg10.genes, batch=colData(spe)$sample_id, type="lower", nmads=3)

colData(spe)$chrM.ratio_3MAD.outlier_slide = isOutlier(colData(spe)$expr_chrM_ratio, batch=colData(spe)$slide, type="higher", nmads=3)
colData(spe)$chrM.ratio_3MAD.outlier_sample = isOutlier(colData(spe)$expr_chrM_ratio, batch=colData(spe)$sample_id, type="higher", nmads=3)

cat(paste(format(Sys.time(), tz="UTC"), "Calculate local outliers (UMI counts)...","\n"))
spe <- localOutliers(spe, metric = "sum_umi", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "sum_local.outlier"

cat(paste(format(Sys.time(), tz="UTC"), "Calculate local outliers (n genes)...","\n"))
spe <- localOutliers(spe, metric = "sum_gene", direction = "lower", log = TRUE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "genes_local.outlier"

cat(paste(format(Sys.time(), tz="UTC"), "Calculate local outliers (mito %)...","\n"))
spe <- localOutliers(spe, metric = "expr_chrM_ratio", direction = "higher", log = FALSE)
#colData(spe)[,c(ncol(colData(spe)),ncol(colData(spe))-2)] = NULL
colnames(colData(spe))[ncol(colData(spe))-1] = "chrM.ratio_local.outlier"

cat(paste(format(Sys.time(), tz="UTC"), "Save spe object to", here("processed-data","03_QC","spe_demo-filt.Rdata"),"\n"))
save(spe, file=here("processed-data","03_QC","spe_demo-filt.Rdata"))

write(c(paste("Modified spe_demo-filt on",format(Sys.time(), tz="UTC"),"UTC"),
        paste("File location:",here("processed-data","03_QC","spe_demo-filt.Rdata")),
        paste("Source code:",here("code","03_QC","02_outliers.r")),
        "*","*","*"), here("spe_tracker_current.txt"), append=TRUE)
## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
