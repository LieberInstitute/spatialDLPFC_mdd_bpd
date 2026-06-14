setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
})

spe <- readRDS("processed-data/publication/spe_n24_example-samples.rds")

# aucell for modules
aucell = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_AUCell.csv", row.names=1)
colnames(aucell) = gsub("Regulon\\.for\\.", "", colnames(aucell))

aucell_mod = aucell[rownames(aucell) %in% colnames(spe),1:(grep("seurat_label", colnames(aucell))-1)]

## no spots missing from aucell_mod
#missing.modules = setdiff(colnames(spe), rownames(aucell_mod))
#length(missing.modules)

reducedDim(spe, "AUCell_module") <- aucell_mod[colnames(spe),]

# aucell for regulons
aucell = read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))

## filter to regulons with at least 10 components
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
reg_subset = regulons$TF[regulons$set_size>=10]

aucell_reg = aucell[rownames(aucell) %in% colnames(spe), reg_subset]

## no spots missing from aucell_reg
#missing.regulons = setdiff(colnames(spe), rownames(aucell_reg))
#length(missing.regulons)

reducedDim(spe, "AUCell_regulon") <- aucell_reg[colnames(spe),]

saveRDS(spe, "processed-data/publication/spe_n24_example-samples.rds")


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
