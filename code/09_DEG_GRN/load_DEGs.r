library(dplyr)

cpList <- readRDS("plots/colorPalettes.rds")
cpList$transfer.bright = cpList$transfer.bright[c(1:4,8,5:7)]

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

deList <- lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  if(results_set=="smoothed-k9-1663") {
	comp_names= names(cpList$smoothed.bright)
	set_name="sm"
  }
  if(results_set=="seurat-pc30") {
	comp_names= names(cpList$transfer.bright)
	set_name="se"
  }

  la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  
  lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                            "_dx-sex_degs-F-test-t-test.csv"))
  
  lat = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
	   # artificially increase the pval of any gene that doesn't have a sig interaction term
	   adj.P.Val2=ifelse(gene_name %in% la.degs$gene_name, adj.P.Val, .5),
           sex.group=factor(coef, levels=comparisons2),
           cluster="L-A", source=set_name)
  
  lrt = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_",
                        results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(dir=factor(sign(logFC), levels=c(1,-1), labels=c("Inc.","Dec.")), 
           # artificially increase the pval of any gene that doesn't have a sig interaction term 
	   adj.P.Val2=ifelse(gene_name %in% lr.degs$gene_name, adj.P.Val, .5),
           cluster=factor(cluster, levels=comp_names), source=set_name,
           sex.group = factor(paste(sex, group, sep="_"), levels=comparisons2))
  
  return(list("la"=lat, "lr"=lrt))
})

names(deList) <- c("sm","se")
deList = do.call(c, deList)

de.df = do.call(rbind, deList) %>%
  mutate(cluster_source = paste(cluster, source))

sigList <- lapply(deList, function(x) {
	filter(x, adj.P.Val2<.05)
})

sig.df = do.call(rbind, sigList) %>%
  mutate(cluster_source = paste(cluster, source))

cpList <- readRDS("plots/colorPalettes.rds")
