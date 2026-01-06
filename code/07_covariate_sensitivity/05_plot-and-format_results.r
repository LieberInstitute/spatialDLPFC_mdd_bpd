setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
})

set.seed(123)
source("code/07_dx_DE/custom_functions.r")
cpList = readRDS("plots/colorPalettes.rds")

covars = c("age","BMI","Smoking","RIN","PMI")
names(covars) <- covars

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
results_set = "smoothed-k9-1663"
comp_names = c("L1","L2","L3dot4","L5","L6","WM")
names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
#results_set = "seurat-pc30"
#comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
#names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")

cat("\nAnnotation results set:", results_set, "\n\n")

# layer agnostic
cat("\n********* Layer-agnostic model ********\n")
resList <- lapply(covars, function(x) 
	readRDS(paste0("processed-data/07_covariate_sensitivity/lmFit-voom_layer-agnostic_", results_set, 
		"_condition-sex_rev-gene-input_covars-pc3-", x, ".rda")))

ebList <- lapply(covars, function(x) eBayes(groupContrasts(resList[[x]], comparisons, 
	add_covar=c("pc3",x)), trend=T))

## pc3 sig
cat("\nMain effect of pc3 covariate in each model:\n")
do.call(rbind, lapply(ebList, function(x) {
  pc3.sig = topTable(x, coef="pc3", p.value=.05, n=Inf)
  c("n_sig"=nrow(pc3.sig), "min"=round(min(pc3.sig$logFC),3), "max"=round(max(pc3.sig$logFC),3))
  }))


## covar sig
cat("\nMain effect of tested covariate in each model:\n")
do.call(rbind, lapply(covars, function(x) {
  var.sig = topTable(ebList[[x]], coef=x, p.value=.05, n=Inf)
  c("n_sig"=nrow(var.sig), "min"=round(min(var.sig$logFC),3), "max"=round(max(var.sig$logFC),3))
}))

## dx*sex effect
sex.resList = lapply(covars, function(x) sexTopTable(ebList[[x]], add_covar=c("pc3",x), phist=T))

pdf(file=paste0("plots/07_covariate_sensitivity/layer-agnostic_", results_set, "_rev-gene-input_p-val-histogram.pdf"))
sex.resList[["age"]]$phist+labs(title="Age", subtitle=paste0("Layer-agnostic (", results_set, ")"))
sex.resList[["BMI"]]$phist+labs(title="BMI", subtitle=paste0("Layer-agnostic (", results_set, ")"))
sex.resList[["Smoking"]]$phist+labs(title="Smoking", subtitle=paste0("Layer-agnostic (", results_set, ")"))
sex.resList[["RIN"]]$phist+labs(title="RIN", subtitle=paste0("Layer-agnostic (", results_set, ")"))
sex.resList[["PMI"]]$phist+labs(title="PMI", subtitle=paste0("Layer-agnostic (", results_set, ")"))
dev.off()

cat("\ndx*sex p value histogram saved to:", paste0("plots/07_covariate_sensitivity/layer-agnostic_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"),"\n")

## format and save results
res.df = do.call(rbind, lapply(covars, function(x) mutate(sex.resList[[x]]$results, covar=x)))

cat("\ndx*sex DEGs (padj<.05, abs(logFC)>.2):\n")
filter(res.df, adj.P.Val<.05, abs(logFC)>.2) %>% 
	group_by(coef, covar) %>% tally() %>%
	tidyr::pivot_wider(names_from="covar", values_from="n")

write.csv(res.df, paste0("processed-data/07_covariate_sensitivity/layer-agnostic_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), row.names=F)
cat("\n\nSaved results dframe to:", paste0("processed-data/07_covariate_sensitivity/layer-agnostic_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), "\n\n")


# layer adjusted
cat("\n********* Layer-adjusted model ********\n")
resList <- lapply(covars, function(x) 
	readRDS(paste0("processed-data/07_covariate_sensitivity/lmFit-voom_layer-adjusted_", results_set, 
		"_condition-sex_rev-gene-input_covars-pc3-", x, ".rda")))

ebList <- lapply(covars, function(x) eBayes(groupContrasts(resList[[x]], comparisons, 
	add_covar=c("pc3",x)), trend=T))

## pc3 sig
cat("\nMain effect of pc3 covariate in each model:\n")
do.call(rbind, lapply(ebList, function(x) {
  pc3.sig = topTable(x, coef="pc3", p.value=.05, n=Inf)
  c("n_sig"=nrow(pc3.sig), "min"=round(min(pc3.sig$logFC),3), "max"=round(max(pc3.sig$logFC),3))
  }))


## covar sig
cat("\nMain effect of tested covariate in each model:\n")
do.call(rbind, lapply(covars, function(x) {
  var.sig = topTable(ebList[[x]], coef=x, p.value=.05, n=Inf)
  c("n_sig"=nrow(var.sig), "min"=round(min(var.sig$logFC),3), "max"=round(max(var.sig$logFC),3))
}))

## dx*sex effect
sex.resList = lapply(covars, function(x) sexTopTable(ebList[[x]], add_covar=c("pc3",x), phist=T))

pdf(file=paste0("plots/07_covariate_sensitivity/layer-adjusted_", results_set, "_rev-gene-input_p-val-histogram.pdf"))
sex.resList[["age"]]$phist+labs(title="Age", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["BMI"]]$phist+labs(title="BMI", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["Smoking"]]$phist+labs(title="Smoking", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["RIN"]]$phist+labs(title="RIN", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["PMI"]]$phist+labs(title="PMI", subtitle=paste0("Layer-adjusted (", results_set, ")"))
dev.off()

cat("\ndx*sex p value histogram saved to:", paste0("plots/07_covariate_sensitivity/layer-adjusted_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"),"\n")

## format and save results
res.df = do.call(rbind, lapply(covars, function(x) mutate(sex.resList[[x]]$results, covar=x)))

cat("\ndx*sex DEGs (padj<.05, abs(logFC)>.2):\n")
filter(res.df, adj.P.Val<.05, abs(logFC)>.2) %>% 
	group_by(coef, covar) %>% tally() %>%
	tidyr::pivot_wider(names_from="covar", values_from="n")

write.csv(res.df, paste0("processed-data/07_covariate_sensitivity/layer-adjusted_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), row.names=F)
cat("\n\nSaved results dframe to:", paste0("processed-data/07_covariate_sensitivity/layer-adjusted_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), "\n\n")

# layer restricted
cat("\n********* Layer-restricted model ********\n")
resList <- lapply(covars, function(x) 
	readRDS(paste0("processed-data/07_covariate_sensitivity/lmFit-voom_layer-restricted_", results_set, 
		"_condition-sex_rev-gene-input_covars-pc3-", x, ".rda")))

#need to have "dot" instead of "." in comparisons for cluster name 
comparisons = unlist(lapply(comp_names, 
                            function(x) paste(x, comparisons, sep="_")))

ebList <- lapply(covars, function(x) eBayes(groupContrasts(resList[[x]], comparisons, 
	add_covar=c("pc3",x)), trend=T))

## pc3 sig
cat("\nMain effect of pc3 covariate in each model:\n")
do.call(rbind, lapply(ebList, function(x) {
  pc3.sig = topTable(x, coef="pc3", p.value=.05, n=Inf)
  c("n_sig"=nrow(pc3.sig), "min"=round(min(pc3.sig$logFC),3), "max"=round(max(pc3.sig$logFC),3))
  }))


## covar sig
cat("\nMain effect of tested covariate in each model:\n")
do.call(rbind, lapply(covars, function(x) {
  var.sig = topTable(ebList[[x]], coef=x, p.value=.05, n=Inf)
  c("n_sig"=nrow(var.sig), "min"=round(min(var.sig$logFC),3), "max"=round(max(var.sig$logFC),3))
}))

## dx*sex effect
### get top table for all comparisons
sex.resList <- lapply(covars, function(y) {
  lapply(names(comp_names), function(x) sexTopTable(ebList[[y]], phist=T, add_covar=c("pc3",y), cluster=x))
})

sex.res_dfList <- lapply(covars, function(x) {
  sex.res_df = do.call(rbind, lapply(sex.resList[[x]], function(y) y$results)) %>%
    mutate(cluster=factor(cluster, levels=names(comp_names)),
           covar=x)
  return(sex.res_df)
})

### p value histogram
plist <- lapply(covars, function(x) {
  tmp = mutate(sex.res_dfList[[x]], sex.group=factor(paste(group, sex), levels=c("NTC.MDD F","NTC.MDD M","NTC.BPD F","NTC.BPD M","MDD.BPD F","MDD.BPD M")))
  ggplot(tmp, aes(x=P.Value))+
    geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
    facet_grid(rows=vars(cluster), cols=vars(sex.group), scales="free_y")+
    labs(title=x, subtitle=paste0("Layer-restricted (", results_set, ")"))+
    theme_minimal()+theme(panel.grid.minor=element_blank())
})
pdf(file=paste0("plots/07_covariate_sensitivity/layer-restricted_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"))
plist[["age"]]
plist[["BMI"]]
plist[["Smoking"]]
plist[["RIN"]]
plist[["PMI"]]
dev.off()

cat("\ndx*sex p value histogram saved to:", paste0("plots/07_covariate_sensitivity/layer-restricted_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"),"\n")

### format and save results
res.df = do.call(rbind, sex.res_dfList) %>%
  mutate(covar= factor(covar, levels=covars))

cat("\ndx*sex DEGs (padj<.05, abs(logFC)>.2):\n")
for(i in levels(res.df$cluster)) {
  cat("\n",i,"\n")
  print(filter(res.df, cluster==i, adj.P.Val<.05, abs(logFC)>.2) %>% 
    group_by(sex, group, covar, .drop=F) %>% tally() %>%
    tidyr::pivot_wider(names_from="covar", values_from="n")
  )
}

write.csv(res.df, paste0("processed-data/07_covariate_sensitivity/layer-restricted_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), row.names=F)
cat("\n\nSaved results dframe to:", paste0("processed-data/07_covariate_sensitivity/layer-restricted_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), "\n\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
