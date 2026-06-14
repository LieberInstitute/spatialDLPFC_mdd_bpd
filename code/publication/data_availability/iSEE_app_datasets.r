setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
})


for(results_set in c("smoothed-k9-1663","seurat-pc30")) {
  
  if(results_set=="smoothed-k9-1663") {
    load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
    annot_name = "domain-SP"
  }
  if(results_set=="seurat-pc30") {
    load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
    annot_name = "domain-CT"
  }
  
  colnames(colData(spe_pseudo))[grep("seurat_label|smoothed_k9_1663", colnames(colData(spe_pseudo)))] = "domain"

  # subset to genelist with DE
  f.df = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set, "_rev-gene-input_F-test.csv"))
  spe_pseudo = spe_pseudo[f.df$gene_id,]
  
  # extract rowData
  rdata = rowData(spe_pseudo)[,c("gene_id","gene_name","gene_type")]
  
  
  # add whole-tissue model stats to rdata
  tmp = f.df[,c("F","adj.P.Val","gene_id","gene_name")]
  colnames(tmp)[1:2] = c("whole.tissue_F_stat","whole.tissue_F_adj.P.Val")
  rdata2 = merge(rdata, tmp, sort=F)
  
  t.df = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set, "_rev-gene-input_moderated-t-test.csv"))
  tmp = select(t.df, adj.P.Val, gene_id, gene_name, coef) %>%
    tidyr::pivot_wider(names_from="coef", values_from="adj.P.Val", 
                       names_prefix = "whole.tissue_t_adj.P.Val_")
  rdata2 = merge(rdata2, tmp, sort=F)
  
  
  # add domain-restricted model stats to rdata
  f.df2 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set, "_rev-gene-input_F-test_all.csv"))
  tmp = f.df2[,c("F","adj.P.Val","gene_id","gene_name")]
  colnames(tmp)[1:2] = c("domain.restricted_F_stat","domain.restricted_F_adj.P.Val")
  rdata2 = merge(rdata2, tmp, sort=F)
  
  
  t.df2 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set, "_rev-gene-input_moderated-t-test.csv"))
  tmp = select(t.df2, adj.P.Val, gene_id, gene_name, coef) %>%
    tidyr::pivot_wider(names_from="coef", values_from="adj.P.Val", 
                       names_prefix = "domain.restricted_t_adj.P.Val_")
  rdata2 = merge(rdata2, tmp, sort=F)
  
  
  # add ammended rdata back to spe
  rownames(rdata2) = rdata2$gene_id
  stopifnot(identical(rownames(rowData(spe_pseudo)), rownames(rdata2)))
  rowData(spe_pseudo) <- rdata2
  
  # remove counts assay to save space
  counts(spe_pseudo) <- NULL
  
  print(head(colData(spe_pseudo)))
  
  saveRDS(spe_pseudo, paste0("processed-data/publication/iSEE_pseudobulk-spe_donor-", annot_name, ".rds"))
  cat("\nSaved rds to:",paste0("processed-data/publication/iSEE_pseudobulk-spe_donor-", annot_name, ".rds"),"\n")
}



cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

