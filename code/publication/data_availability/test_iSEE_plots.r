setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')

suppressPackageStartupMessages({
  library(SpatialExperiment)
})

#test_gene = "CORT"
test_gene = "COL14A1"

set.seed(1234)

# read in custom functions
source("code/publication/data_availability/iSEE_custom-plot_utils.r")

# read in spe
spe <- readRDS("processed-data/publication/iSEE_pseudobulk-spe_both-annotations.rds")

# swap rownames, update metadata rownames
rownames(spe) <- rowData(spe)$gene_name
rownames(metadata(spe)[[1]]) <- metadata(spe)[[1]]$gene_name
rownames(metadata(spe)[[2]]) <-	metadata(spe)[[2]]$gene_name


p1 <- CUSTOM_VIOLIN(spe, test_gene, mode="whole-tissue")
p2 <- CUSTOM_VIOLIN(spe, test_gene, mode="domain-SP")
p3 <- CUSTOM_VIOLIN(spe, test_gene, mode="domain-CT")

ggsave(file=paste0("plots/publication/test_iSEE-revised_", test_gene, ".pdf"), 
	gridExtra::grid.arrange(p1, p2, p3, ncol=1),
	width=7, height=10)

# NOT RUN, template for use within iSEE
#library(iSEE)
#GENERATOR <- createCustomPlot(CUSTOM_VIOLIN)
#custom.p1 <- GENERATOR(mode="whole-tissue")
#custom.p2 <- GENERATOR(mode="domain-CT")
#custom.p3 <- GENERATOR(mode="domain-SP")
#app <- iSEE(spe, initial=list(custom.p1, custom.p2, custom.p3))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
