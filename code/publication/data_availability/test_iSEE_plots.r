setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')

suppressPackageStartupMessages({
  library(SpatialExperiment)
})

set.seed(1234)

# read in custom functions
source("code/publication/data_availability/iSEE_custom-plot_utils.r")

# create test CORT plot for each pseudobulk dataset
domain1 = c("domain-SP","domain-CT")
plist = lapply(domain1, function(x) {
	spe <- readRDS(paste0("processed-data/publication/iSEE_pseudobulk-spe_donor-", x, ".rds"))	

	suppressMessages(iSEEplots(spe, "CORT", annot_name=domain1))
})

ggsave(file="plots/publication/test_iSEE_CORT.pdf", 
	gridExtra::arrangeGrob(grobs=plist, ncol=1, top=NULL),
	width=8, height=7)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
