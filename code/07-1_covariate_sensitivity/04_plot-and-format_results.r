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

covars = c("detected","nspots","age","BMI","Smoking","RIN","PMI")
names(covars) <- covars

#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
#results_set = "smoothed-k9-1663"
#comp_names = c("L1","L2","L3dot4","L5","L6","WM")
#names(comp_names) = c("L1","L2","L3.4","L5","L6","WM")

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
results_set = "seurat-pc30"
comp_names = c("MicrodotVasc","Astro","L2dot3","L4","Inhb","L5","L6","Oligo")
names(comp_names) = c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")

cat("\nAnnotation results set:", results_set, "\n\n")

# layer adjusted
cat("\n********* Layer-adjusted model ********\n")
resList <- lapply(covars, function(x) 
	readRDS(paste0("processed-data/07-1_covariate_sensitivity/lmFit-voom_layer-adjusted_", results_set, 
		"_condition-sex_rev-gene-input_covars-pc3-", x, ".rda")))

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")

ebList <- lapply(covars, function(x) eBayes(groupContrasts(resList[[x]], comparisons, 
	add_covar=c("pc3",x)), trend=T, robust=T))

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
### F statistics
ebList_f <- lapply(covars, function(x) {
  getTopTable(ebList[[x]], .coef=comparisons) %>%
	mutate("F_design"="~ 0 + dx*sex",
		covar=x)
})
names(ebList_f) <- covars

ebList_f_sig = lapply(ebList_f, function(x) filter(x, adj.P.Val<.05)$gene_id)

eb_f_df = do.call(rbind, ebList_f)

### p hist
test1 = do.call(rbind, ebList_f) %>%
  mutate(covar=factor(covar, levels=covars))
ph = ggplot(test1, aes(x=P.Value))+
  geom_histogram(bins=30)+
  facet_wrap(vars(covar))+
  scale_x_continuous(breaks=c(0,.5,1))+
  labs(subtitle=paste0("Layer-adjusted (", results_set, ")"),
    title="~ 0 + dx*sex", x="F test p value")+
  theme_minimal()+theme(panel.grid.minor=element_blank())

## moderated t statistics
sex.resList = lapply(covars, function(x) {
  tmp = sexTopTable(ebList[[x]], add_covar=c("pc3",x), phist=T)
  ##### have to add covars because they are part of the contrast matrix but it doesn't change that the comparisons are just dx*sex

  ### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
  ### details: https://github.com/cran/limma/blob/master/R/decidetests.R
  tmp$results$adj.P.Val = p.adjust(tmp$results$P.Value, method="BH")
  tmp$results$covar = x
  return(tmp)
})

cat("\nLayer-adjusted analysis with adj p<.05 and F adj p<.05:\n\n")
do.call(rbind, lapply(covars, function(x) filter(sex.resList[[x]]$results, adj.P.Val<.05, gene_id %in% ebList_f_sig[[x]]) %>% 
	group_by(sex, group, covar, .drop=F) %>% tally()
)) %>% mutate(covar= factor(covar, levels=covars)) %>%
	tidyr::pivot_wider(names_from="covar", values_from="n", values_fill=0)


pdf(file=paste0("plots/07-1_covariate_sensitivity/layer-adjusted_", results_set, "_rev-gene-input_p-val-histogram.pdf"))
ph
sex.resList[["detected"]]$phist+labs(title="detected", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["nspots"]]$phist+labs(title="nspots", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["age"]]$phist+labs(title="Age", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["BMI"]]$phist+labs(title="BMI", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["Smoking"]]$phist+labs(title="Smoking", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["RIN"]]$phist+labs(title="RIN", subtitle=paste0("Layer-adjusted (", results_set, ")"))
sex.resList[["PMI"]]$phist+labs(title="PMI", subtitle=paste0("Layer-adjusted (", results_set, ")"))
dev.off()

cat("\nLayer-adjusted dx*sex p value histogram saved to:", paste0("plots/07-1_covariate_sensitivity/layer-adjusted_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"),"\n")

write.csv(eb_f_df, paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", results_set,
	"_rev-gene-input_F-test_all-covar-models.csv"), row.names=F)
cat("\n\n\nSaved F test results dframe to:", paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", results_set,
        "_rev-gene-input_F-test_all-covar-models.csv"), "\n")


res.df = do.call(rbind, lapply(covars, function(x) mutate(sex.resList[[x]]$results, covar=x)))

write.csv(res.df, paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", results_set, 
	"_rev-gene-input_moderated-t-test_all-covar-models.csv"), row.names=F)
cat("\nSaved moderated t test results dframe to:", paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", results_set, 
	"_rev-gene-input_moderated-t-test_all-covar-models.csv"), "\n\n")


# layer restricted
cat("\n********* Layer-restricted model ********\n")
resList <- lapply(covars, function(x) 
	readRDS(paste0("processed-data/07-1_covariate_sensitivity/lmFit-voom_layer-restricted_", results_set, 
		"_condition-sex_rev-gene-input_covars-pc3-", x, ".rda")))

#need to have "dot" instead of "." in comparisons for cluster name 
comparisons2 = unlist(lapply(comp_names, 
                            function(x) paste(x, comparisons, sep="_")))

ebList <- lapply(covars, function(x) eBayes(groupContrasts(resList[[x]], comparisons2, 
	add_covar=c("pc3",x)), trend=T, robust=T))

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
### F statistics
comparisons2_names = unlist(lapply(names(comp_names),
                            function(x) paste(x, comparisons, sep="_")))
eb_fList_all = lapply(covars, function(x) 
  getTopTable(ebList[[x]], .coef=comparisons2_names) %>%
	mutate(cluster="all", F_design="~ 0 + dx*sex*cluster",
		covar=x)
)

f_sig_list_all = lapply(eb_fList_all, function(x)
  filter(x, adj.P.Val<.05)$gene_id
)

eb_fList_all.df = do.call(rbind, eb_fList_all) %>%
  mutate(covar=factor(covar, levels=covars))
ph2 <- ggplot(eb_fList_all.df, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  facet_wrap(vars(covar), scales="free_y")+
  labs(subtitle=paste0("Layer-restricted (", results_set, ")"),
    title="~ 0 + dx*sex*cluster", x="F test p value")+
  theme_minimal()+theme(panel.grid.minor=element_blank())

#### nested approach
eb_fList_each <- lapply(covars, function(z) {
  tmp = lapply(names(comp_names), function(x) {
    tmp = getTopTable(ebList[[z]], .coef=grep(x, colnames(coef(ebList[[z]])), value=T))
    colnames(tmp)[1:6] = sapply(strsplit(colnames(tmp)[1:6], "_"), function(y) paste(y[[2]], y[[3]], sep="_"))
    tmp$cluster = x
    tmp$covar = z
    tmp$F_design = paste0(x, ": ~ 0 + dx*sex")
    return(tmp)
  })
  ### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
  ### details: https://github.com/cran/limma/blob/master/R/decidetests.R
  tmp = do.call(rbind, tmp) %>% 
    mutate(cluster=factor(cluster, levels=names(comp_names)))
  tmp$adj.P.Val_global = p.adjust(tmp$P.Value, method="BH")
  return(tmp)
})

#f_sig_list_each = lapply(eb_fList_each, function(z) {
#  tmp = lapply(names(comp_names), function(x)
#	filter(z, cluster==x, adj.P.Val<.05)$gene_id
#  )
#  names(tmp) <- names(comp_names)
#  return(tmp)
#})

eb_fList_each.df = do.call(rbind, eb_fList_each) %>%
  mutate(covar=factor(covar, levels=covars))

ph3 <- ggplot(eb_fList_each.df, aes(x=P.Value))+
  geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
  facet_grid(rows=vars(cluster), cols=vars(covar), scales="free_y")+
  labs(subtitle=paste0("Layer-restricted (", results_set, ")"),
    title="cluster: ~ 0 + dx*sex", x="F test p value")+
  theme_minimal()+theme(panel.grid.minor=element_blank())


## moderated t statistics
sex.resList <- lapply(covars, function(y) {
  tmp = do.call(rbind, lapply(names(comp_names), function(x) sexTopTable(ebList[[y]], phist=F, add_covar=c("pc3",y), cluster=x) %>%
	mutate(covar=y))) %>%
  mutate(cluster= factor(cluster, levels=names(comp_names)))
  ### implement global padj method (without this, different t cutoffs were sig because of different group sizes)
  ### details: https://github.com/cran/limma/blob/master/R/decidetests.R
  tmp$adj.P.Val = p.adjust(tmp$P.Value, method="BH")
  return(tmp)
})


### p value histogram
plist <- lapply(covars, function(x) {
  tmp = mutate(sex.resList[[x]], sex.group=factor(paste(group, sex), levels=c("NTC.MDD F","NTC.MDD M","NTC.BPD F","NTC.BPD M","MDD.BPD F","MDD.BPD M")))
  ggplot(tmp, aes(x=P.Value))+
    geom_histogram(bins=30)+scale_x_continuous(breaks=c(0,.5,1))+
    facet_grid(rows=vars(cluster), cols=vars(sex.group), scales="free_y")+
    labs(title=x, subtitle=paste0("Layer-restricted (", results_set, ")"))+
    theme_minimal()+theme(panel.grid.minor=element_blank())
})
pdf(file=paste0("plots/07-1_covariate_sensitivity/layer-restricted_", results_set, 
	"_rev-gene-input_p-val-histogram.pdf"), width=8, height=11)
gridExtra::grid.arrange(ph2, ph3, layout_matrix=matrix(c(1,2,2,2)))
plist[["detected"]]
plist[["nspots"]]
plist[["age"]]
plist[["BMI"]]
plist[["Smoking"]]
plist[["RIN"]]
plist[["PMI"]]
dev.off()

cat("\n\nLayer-adjusted dx*sex p value histograms saved to:", paste0("plots/07-1_covariate_sensitivity/layer-restricted_", results_set,
        "_rev-gene-input_p-val-histogram.pdf"), "\n")


### format and save results
cat("\nLayer-restricted dx*sex DEGs (padj<.05, F padj<.05):\n")
for(i in covars) {
  cat("\n",i,"\n")
  print(filter(sex.resList[[i]], adj.P.Val<.05, gene_id %in% f_sig_list_all[[i]]) %>% 
    group_by(sex, group, cluster) %>% tally() %>%
    tidyr::pivot_wider(names_from="cluster", values_from="n", values_fill=0))
}


write.csv(eb_fList_all.df, paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", results_set,
        "_rev-gene-input_F-test_all_all-covar-models.csv"), row.names=F)
write.csv(eb_fList_each.df, paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", results_set,
        "_rev-gene-input_F-test_each_all-covar-models.csv"), row.names=F)
cat("\n\n\nSaved F test results dframe to:", paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", results_set,
        "_rev-gene-input_F-test_all-covar-models.csv"), "\n")


res.df = do.call(rbind, sex.resList)

write.csv(res.df, paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", results_set,
        "_rev-gene-input_moderated-t-test_all-covar-models.csv"), row.names=F)
cat("\nSaved moderated t test results dframe to:", paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", results_set,
        "_rev-gene-input_moderated-t-test_all-covar-models.csv"), "\n\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
