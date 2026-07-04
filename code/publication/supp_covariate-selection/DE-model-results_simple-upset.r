setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(UpSetR)
})
set.seed(123)

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons
comparisons2 = comparisons[c(1,3,5,2,4,6)]

covars = c("none","pc3 only","nspots","chrM_ratio","age","BMI","Smoking","RIN")
names(covars) <- covars

# layer adjusted
la.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  ## F test
  f1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                       results_set, "_rev-gene-input_F-test.csv")) %>%
    mutate(covar="none")
  f2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                       results_set, "_rev-gene-input_F-test.csv")) %>%
    mutate(covar="pc3 only")
  f3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                       results_set, "_rev-gene-input_F-test_all-covar-models.csv"))
  
  la.all = bind_rows(f1, f2, f3) %>% 
    filter(covar %in% covars) %>%
    mutate(covar=factor(covar, levels=covars), annotation=results_set)
  return(la.all)
})) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                 labels=c("domain-SP","domain-CT")))

plot.df = filter(la.all, adj.P.Val<.05) %>% group_by(covar, annotation) %>% tally()

sp.upset = lapply(covars, function(x) {
	filter(la.all, adj.P.Val<.05, covar==x, annotation=="domain-SP")$gene_name
})

ct.upset = lapply(covars, function(x) {
        filter(la.all, adj.P.Val<.05, covar==x, annotation=="domain-CT")$gene_name 
})

## t test
lat.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  t1 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-no-covars_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="none")
  
  t2 = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="pc3 only")
  
  t3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-adjusted_", 
                       results_set, "_rev-gene-input_moderated-t-test_all-covar-models.csv"))
  
  la.t.all = bind_rows(t1, t2, t3) %>% 
    filter(covar %in% covars) %>%
    mutate(coef= factor(coef, levels=comparisons),
           covar=factor(covar, levels=covars),
           annotation= results_set)
  return(la.t.all)
})) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                labels=c("domain-SP","domain-CT")))

### compare number of DEGs
la.degs = do.call(rbind, lapply(covars, function(y) {
  sig.genes = filter(lat.all, covar==y, annotation=="domain-SP", adj.P.Val<.05)$gene_id
  f1 = filter(la.all, covar==y, adj.P.Val<.05, annotation=="domain-SP", gene_id %in% sig.genes)
  sig.genes = filter(lat.all, covar==y, annotation=="domain-CT", adj.P.Val<.05)$gene_id
  f2 = filter(la.all, covar==y, adj.P.Val<.05, annotation=="domain-CT", gene_id %in% sig.genes)
  return(rbind(f1, f2))
})) 

sp.upset = lapply(covars, function(x) {
        filter(la.degs, covar==x, annotation=="domain-SP")$gene_name
})
names(sp.upset)[2] = "pc3_only"

cat("\nNumber of DEGs (whole-tissue, domain-SP):\n")
print(sapply(sp.upset, length))


getIntersections <- function(upsetList) {

	totally.unique = lapply(names(upsetList), function(x) {
		ref1 = upsetList[[x]]
		q1 = unlist(upsetList[setdiff(names(upsetList), x)])
		return(setdiff(ref1, q1))
	})
	names(totally.unique) <- names(upsetList)

	all.overlaps = table(unlist(upsetList))
	all.overlaps = names(all.overlaps)[all.overlaps==length(upsetList)]

	pc3.overlaps = table(unlist(upsetList[-1]))
	pc3.overlaps = setdiff(names(pc3.overlaps)[pc3.overlaps==(length(upsetList)-1)], all.overlaps)

	data.frame(group=c(names(upsetList), "all_models", "pc3_models"),
		n_DEGs= c(sapply(totally.unique, length), length(all.overlaps), length(pc3.overlaps)))
		

}

df1 = getIntersections(sp.upset) %>% mutate(annotation="domain-SP", model="whole-tissue")

ct.upset = lapply(covars, function(x) {
        filter(la.degs, covar==x, annotation=="domain-CT")$gene_name
})
names(ct.upset)[2] = "pc3_only"

cat("\nNumber of DEGs (whole-tissue, domain-CT):\n")
print(sapply(ct.upset, length))

df2 = getIntersections(ct.upset) %>% mutate(annotation="domain-CT", model="whole-tissue")


# layer-restricted
lr.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  ## F test
  f1.all = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                           results_set, "_rev-gene-input_F-test_all.csv")) %>%
    mutate(covar="none")
  
  f2.all = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                           results_set, "_rev-gene-input_F-test_all.csv")) %>%
    mutate(covar="pc3 only")
  
  f3.all = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                           results_set, "_rev-gene-input_F-test_all_all-covar-models.csv"))
  
  lr.all = bind_rows(f1.all, f2.all, f3.all) %>%
    filter(covar %in% covars) %>%
    mutate(covar=factor(covar, levels=covars), annotation= results_set)
  return(lr.all[,c("AveExpr","F","P.Value","adj.P.Val","df.prior","df.total","gene_id","gene_name",
                   "cluster","F_design","covar","annotation")])
})) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                labels=c("domain-SP","domain-CT")))



#sp.upset = lapply(covars, function(x) {
#        filter(lr.all, adj.P.Val<.05, covar==x, annotation=="domain-SP")$gene_name
#})

#ct.upset = lapply(covars, function(x) {
#        filter(lr.all, adj.P.Val<.05, covar==x, annotation=="domain-CT")$gene_name
#})

## t test
lrt.all = do.call(rbind, lapply(c("smoothed-k9-1663","seurat-pc30"), function(results_set) {
  if(results_set == "smoothed-k9-1663") clus_levels=c("L1","L2","L3.4","L5","L6","WM")
  if(results_set == "seurat-pc30") clus_levels=c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
  
  t1 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-no-covars_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="none")
  
  t2 = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3_", 
                       results_set, "_rev-gene-input_moderated-t-test.csv")) %>%
    mutate(covar="pc3 only")
  
  t3 = read.csv(paste0("processed-data/07-1_covariate_sensitivity/layer-restricted_", 
                       results_set, "_rev-gene-input_moderated-t-test_all-covar-models.csv"))
  
  lr.t.all = bind_rows(t1, t2, t3) %>% 
    filter(covar %in% covars) %>%
    mutate(sex.group= factor(paste(sex, group, sep="_"), levels=comparisons),
           covar=factor(covar, levels=covars),
           cluster=factor(cluster, levels=clus_levels),
           annotation= results_set)
  return(lr.t.all)
  })) %>% mutate(annotation=factor(annotation, levels=c("smoothed-k9-1663", "seurat-pc30"), 
                                   labels=c("domain-SP","domain-CT")))


### compare number of DEGs
lr.degs = do.call(rbind, lapply(covars, function(y) {
  sig.genes = filter(lrt.all, covar==y, annotation=="domain-SP", adj.P.Val<.05)$gene_id
  f1 = filter(lr.all, covar==y, adj.P.Val<.05, annotation=="domain-SP", gene_id %in% sig.genes)
  sig.genes = filter(lrt.all, covar==y, annotation=="domain-CT", adj.P.Val<.05)$gene_id
  f2 = filter(lr.all, covar==y, adj.P.Val<.05, annotation=="domain-CT", gene_id %in% sig.genes)
  return(rbind(f1, f2))
})) 


sp.upset = lapply(covars, function(x) {
        filter(lr.degs, covar==x, annotation=="domain-SP")$gene_name
})
names(sp.upset)[2] = "pc3_only"

cat("\nNumber of DEGs (domain-restricted, domain-SP):\n")
print(sapply(sp.upset, length))

df3 = getIntersections(sp.upset) %>% mutate(annotation="domain-SP", model="domain-restricted")


ct.upset = lapply(covars, function(x) {
        filter(lr.degs, covar==x, annotation=="domain-CT")$gene_name
})
names(ct.upset)[2] = "pc3_only"

cat("\nNumber of DEGs (domain-restricted, domain-CT):\n")
print(sapply(ct.upset, length))


df4 = getIntersections(ct.upset) %>% mutate(annotation="domain-CT", model="domain-restricted")

df.all = bind_rows(df1, df2, df3, df4) %>%
	mutate(group= factor(group, levels=c("all_models", "pc3_models", "none", "pc3_only", "age", "nspots", "BMI", "chrM_ratio", "Smoking", "RIN")),
		annotation = factor(annotation, levels=c("domain-SP",  "domain-CT")),
		model = factor(model, levels=c("whole-tissue", "domain-restricted"))
	)

p1 <- ggplot(df.all, aes(x=group, y=n_DEGs))+
	geom_bar(stat="identity", width=.7)+
	geom_text(aes(label=n_DEGs), hjust=.5, vjust=0, size=2, color="red")+
	scale_y_continuous(breaks=c(0,40,80,120))+
	facet_grid(rows=vars(model), cols=vars(annotation), axes="all_x")+
	theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), text=element_text(size=6),
		panel.grid.minor=element_blank(), axis.title.x=element_blank())

ggsave(file="plots/publication/supp_covariate-selection/DE-model_simple-overlaps.pdf", p1, width=4, height=3.5)


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
