setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
})

# modules first
refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")
write.csv(select(refined.modules, source=TF, target, importance, rho), "processed-data/publication/supp_tables/modules_filtered-refined.csv", row.names=F)

# select only regulons greater than 10 components
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons <- filter(regulons, set_size>10)
write.csv(select(regulons, TF, regulon_size=set_size, regulon_targets=set_str), "processed-data/publication/supp_tables/regulons_n109-min-10.csv", row.names=F)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
