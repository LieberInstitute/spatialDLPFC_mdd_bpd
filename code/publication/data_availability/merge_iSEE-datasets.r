setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
})

spe_sp <- readRDS("processed-data/publication/iSEE_pseudobulk-spe_donor-domain-SP.rds")
spe_sp$annotation= "domain-SP"
rowData(spe_sp) <- rowData(spe_sp)[,c("gene_id","gene_name","gene_type")]

spe_ct <- readRDS("processed-data/publication/iSEE_pseudobulk-spe_donor-domain-CT.rds")
spe_ct$annotation= "domain-CT"
rowData(spe_ct) <- rowData(spe_ct)[,c("gene_id","gene_name","gene_type")]

each.genes = union(rownames(metadata(spe_sp)[[1]]), rownames(metadata(spe_ct)[[1]]))
missing.sp = setdiff(each.genes, rownames(spe_sp))
both.genes = setdiff(each.genes, missing.sp)


attempt1 <- cbind(spe_sp[both.genes,], spe_ct[both.genes,])
attempt1

dim(attempt1)

saveRDS(attempt1, "processed-data/publication/iSEE_pseudobulk-spe_both-annotations.rds")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
