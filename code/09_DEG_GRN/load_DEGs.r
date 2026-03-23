library(dplyr)

cpList <- readRDS("plots/colorPalettes.rds")
cpList$transfer.bright = cpList$transfer.bright[c(1:4,8,5:7)]

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

sigList <- lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  if(results_set=="smoothed-k9-1663") comp_names= names(cpList$smoothed.bright)
  if(results_set=="seurat-pc30") comp_names= names(cpList$transfer.bright)
  
  la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  
  lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  
  lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
           sex.group=factor(coef, levels=comparisons2),
           cluster="L-A") %>%
    filter(gene_id %in% la.degs$gene_id, adj.P.Val<.05)
  
  lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
           cluster=factor(cluster, levels=comp_names),
           sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2)) %>%
    filter(gene_id %in% lr.degs$gene_id, adj.P.Val<.05)
  
  return(list("la"=lat, "lr"=lrt))
})

names(sigList) <- c("sm","se")
sigList = do.call(c, sigList)
