setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(scuttle)
})
set.seed(123)
setAutoBlockSize(1e9)

##load spe
#spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")

##i keep getting a weird error when trying with SpatialExperiment that fdoesn't happen with the Seurat SingleCellExperiment so try to convert to sce first
#sce <- SingleCellExperiment(assays = list(logcounts = logcounts(spe)), colData=colData(spe))
#rowData(sce) = rowData(spe)

##pseudobulk raw counts
#cat("\nPseudobulk spe by: sample_id\n")

#spe_summ = scuttle::aggregateAcrossCells(sce, ids=colData(spe)[,c("sample_id")],
#                                         statistics=c("mean","prop.detected"),
#                                         use.assay.type="logcounts")

#save(spe_summ, file="processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")

cat("\nLoad saved pseudobulk...\n")
load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")
#change name of ncells to nspots
colnames(colData(spe_summ))[grep("ncells", colnames(colData(spe_summ)))] = "nspots"

#keep only sample level coldata
colData(spe_summ) = colData(spe_summ)[,c("sample_id","brnum","age",
  "sex","condition",
  "PMI","RIN","slide","array","MBv_sample","seq","round",
  "nspots")]

colData(spe_summ)$condition = factor(spe_summ$condition, levels=c("NTC","MDD","BPD"))
colData(spe_summ)$sex = factor(spe_summ$sex, levels=c("F","M"))

Sys.time()
save(spe_summ, file="processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")
cat("\nPseudobulk spe saved to: processed-data/06_pseudobulk/spe_n119_pseudo-dotplot_sample-id.Rdata")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
