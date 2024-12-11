#args = commandArgs(TRUE) 
#print(args[[1]]) 
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(nnSVG)
})
set.seed(123)

sample = "V13B23-342.Rdata"
load(file=paste0("processed-data/04_feature_selection/per-slide_spe/",sample))
cat("\nCalculating nnSVG... ",format(Sys.time(),tz="EST"),"\n")
dim(tmp)
results <- nnSVG(tmp, n_threads=12)
svg = rowData(results)
write.csv(svg, paste0("processed-data/04_feature_selection/per-slide_svgs/",gsub("\\.Rdata","_nnSVG-results",sample),".csv"), row.names=T)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
