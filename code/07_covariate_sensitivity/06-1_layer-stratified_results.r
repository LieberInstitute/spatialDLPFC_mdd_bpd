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

#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#results_set = "smoothed-k9-1663"
#comp_names = c("L1","L2","L3dot4","L5","L6","WM")
#names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
results_set = "seurat-pc30"
comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")

# layer stratified
cat("\n********* Layer-stratified model ********\n")
resList <- lapply(covars, function(x) readRDS(paste0("processed-data/07_covariate_sensitivity/lmFit-voom_layer-stratified_", results_set, 
		"_condition-sex_rev-gene-input_covars-pc3-", x, ".rda")))


ebList <- lapply(covars, function(x) {
  lapply(resList[[x]], function(y) eBayes(groupContrasts(y, comparisons, 
                                                         add_covar=c("pc3",x)), trend=T)
         )
})


## dx*sex effect
### get top table for all comparisons
sex.resList <- lapply(covars, function(x) {
	do.call(rbind, lapply(names(comp_names), function(y) 
		mutate(sexTopTable(ebList[[x]][[y]], add_covar=c("pc3",x), phist=F), cluster=y)
		)
	) %>% mutate(cluster=factor(cluster, levels=names(comp_names)), covar=x)
})


### p value histogram
plist <- lapply(covars, function(x) {
  tmp = mutate(sex.resList[[x]], coef=factor(coef, levels=comparisons))
  ggplot(tmp, aes(x=P.Value))+
    geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
    facet_grid(rows=vars(cluster), cols=vars(coef), scales="free_y")+
    labs(title=x, subtitle=paste0("Layer-stratified (", results_set, ")"))+
    theme_minimal()+theme(panel.grid.minor=element_blank())
})
pdf(file=paste0("plots/07_covariate_sensitivity/layer-stratified_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"))
plist[["age"]]
plist[["BMI"]]
plist[["Smoking"]]
plist[["RIN"]]
plist[["PMI"]]
dev.off()

cat("\ndx*sex p value histogram saved to:", paste0("plots/07_covariate_sensitivity/layer-stratified_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"),"\n")

### format and save results
res.df = do.call(rbind, sex.resList) %>%
  mutate(covar= factor(covar, levels=covars))

cat("\ndx*sex DEGs (padj<.05, abs(logFC)>.2):\n")
for(i in levels(res.df$cluster)) {
  cat("\n",i,"\n")
  print(filter(res.df, cluster==i, adj.P.Val<.05, abs(logFC)>.2) %>% 
    group_by(sex, group, covar, .drop=F) %>% tally() %>%
    tidyr::pivot_wider(names_from="covar", values_from="n")
  )
}

write.csv(res.df, paste0("processed-data/07_covariate_sensitivity/layer-stratified_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), row.names=F)
cat("\n\nSaved results dframe to:", paste0("processed-data/07_covariate_sensitivity/layer-stratified_", results_set, 
	"_rev-gene-input_all-covar-model-results.csv"), "\n\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
