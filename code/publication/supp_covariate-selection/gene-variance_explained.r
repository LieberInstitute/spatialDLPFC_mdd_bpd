setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(edgeR)
        library(scater)
        library(dplyr)
        library(ggplot2)
})

set.seed(123)

#need to remove age from colData because of significant digits change
demo = read.csv("processed-data/publication/supp_tables/demographics.csv")



load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
spe_sm <- spe_pseudo

new.cdata = merge(colData(spe_sm)[,c("sample_id","brnum","condition","sex","smoothed_k9_1663","nspots","chrM_ratio","pc3")],
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_sm$pc3, new.cdata$pc3))
colData(spe_sm) <- new.cdata

l1 <- lapply(levels(spe_sm$smoothed_k9_1663), function(x) {
	tmp <- spe_sm[,spe_sm$smoothed_k9_1663==x]
	var.m = as.data.frame(getVarianceExplained(tmp, variables=c("nspots","chrM_ratio","pc3","age","RIN","BMI","Smoking")))
	var.m$cluster = x
	return(var.m)
})

df1 = do.call(rbind, l1) %>% mutate(cluster=factor(cluster, levels=levels(spe_sm$smoothed_k9_1663)))


load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
spe_se <- spe_pseudo

new.cdata = merge(colData(spe_se)[,c("sample_id","brnum","condition","sex","seurat_label","nspots","chrM_ratio","pc3")],
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_se$pc3, new.cdata$pc3))
colData(spe_se) <- new.cdata

l2 <- lapply(levels(spe_se$seurat_label), function(x) {
	tmp <- spe_se[,spe_se$seurat_label==x]
	var.m = as.data.frame(getVarianceExplained(tmp, variables=c("nspots","chrM_ratio","pc3","age","RIN","BMI","Smoking")))
	var.m$cluster = x
	return(var.m)
})

df2 = do.call(rbind, l2) %>% mutate(cluster=factor(cluster, levels=levels(spe_se$seurat_label)))

df_both = bind_rows(mutate(df1, annotation="domain-SP"), mutate(df2, annotation="domain-CT")) %>% 
	mutate(cluster=factor(cluster, levels=c("Micro.Vasc","Astro","L1","L2","L2.3","L3.4","L4","Inhb","L5","L6","WM","Oligo"),
		labels=c("M/V","Ast","L1","L2","L2/3","L3/4","L4","Inb","L5","L6","WM","Olg")),
	annotation=factor(annotation, levels=c("domain-SP","domain-CT")))



col.pal = RColorBrewer::brewer.pal(n=6, "YlOrRd")

df_tally <- do.call(rbind, lapply(c("nspots","chrM_ratio","pc3","age","RIN","BMI","Smoking"), function(x) {
	tmp <- df_both
	tmp$plot.me = cut(tmp[[x]], breaks=c(0,1,5,10,25,50,100), include.lowest=T)

	return(group_by(tmp, plot.me, annotation, cluster) %>% tally() %>% mutate(covar=x))
})) %>% mutate(covar=factor(covar, levels=c("pc3","nspots","chrM_ratio","age","BMI","RIN","Smoking")))

p1 <- ggplot(df_tally, aes(x=cluster, y=n, fill=plot.me))+
	geom_bar(stat="identity", position="stack")+
	scale_fill_brewer("% gene\nvar. (bin)", palette="Greys")+
	facet_grid(cols=vars(annotation), rows=vars(covar), scales="free_x", space="free_x")+
	theme_bw()

ggsave(file="plots/publication/supp_covariate-selection/gene-level_variance_explained.pdf", 
	p1, height=6, width=4)




cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
