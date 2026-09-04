setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
})

set.seed(123)

#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9.Rdata")
#spe_sp <- spe_pseudo

#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30.Rdata")
#spe_ct <- spe_pseudo

spe_sp <- readRDS("processed-data/publication/iSEE_pseudobulk-spe_donor-domain-SP.rds")
spe_sp$annotation= "domain-SP"
rowData(spe_sp) <- rowData(spe_sp)[,c("gene_id","gene_name","gene_type")]

spe_ct <- readRDS("processed-data/publication/iSEE_pseudobulk-spe_donor-domain-CT.rds")
spe_ct$annotation= "domain-CT"
rowData(spe_ct) <- rowData(spe_ct)[,c("gene_id","gene_name","gene_type")]

both.genes = intersect(rownames(spe_sp), rownames(spe_ct))

spe <- cbind(spe_sp[both.genes,], spe_ct[both.genes,])

rownames(spe) <- rowData(spe)$gene_name
rownames(metadata(spe)[[1]]) <- metadata(spe)[[1]]$gene_name
rownames(metadata(spe)[[2]]) <-	metadata(spe)[[2]]$gene_name

saveRDS(spe, "processed-data/publication/MBv-shiny_pseudobulk-spe_both-annotations.rds")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

