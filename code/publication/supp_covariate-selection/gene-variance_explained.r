setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(edgeR)
	library(pheatmap)
        library(scater)
        library(dplyr)
        library(ggplot2)
	library(gridExtra)
})

set.seed(123)


demo = read.csv("processed-data/publication/supp_tables/demographics.csv")


# domain-SP
cat("\ndomain-SP...\n")
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
spe_sm <- spe_pseudo

new.cdata = merge(colData(spe_sm)[,c("sample_id","brnum","condition","sex","smoothed_k9_1663","nspots","chrM_ratio","pc3","slide", "seq")],
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_sm$pc3, new.cdata$pc3))
colData(spe_sm) <- new.cdata

# variance explained correlation
covars = c("sample_id", "nspots", "chrM_ratio", "pc3","age", "BMI", "RIN", "Smoking", "slide", "seq")
var.m = getVarianceExplained(spe_sm, variables=covars, exprs_values="logcounts")
summary(var.m)

cor.var.m = cor(var.m, method="pearson")
cor.var.m[cor.var.m==1] = NA

ordered = c("sample_id","slide","pc3","seq","nspots","age","Smoking","BMI","chrM_ratio","RIN")

phm1 = pheatmap(cor.var.m[ordered, ordered], angle_col=90, main="domain-SP gene var corr",
                cluster_cols=F, cluster_rows=F,
                color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100), na_col="grey50",
                breaks=seq(-1, 1, length.out=101))


# binned explanation
l1 <- lapply(levels(spe_sm$smoothed_k9_1663), function(x) {
	tmp <- spe_sm[,spe_sm$smoothed_k9_1663==x]
	var.m = as.data.frame(getVarianceExplained(tmp, variables=c("nspots","chrM_ratio","pc3","age","RIN","BMI","Smoking")))
	var.m$cluster = x
	return(var.m)
})

df1 = do.call(rbind, l1) %>% mutate(cluster=factor(cluster, levels=levels(spe_sm$smoothed_k9_1663)))


# domain-CT annotation
cat("\n\ndomain-CT...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
colnames(colData(spe_pseudo))[grep("subsets_mito_percent", colnames(colData(spe_pseudo)))] = "chrM_ratio"
spe_pseudo$pc3 = reducedDim(spe_pseudo, "PCA_1663")[,3]
rowData(spe_pseudo)$high_expr_group_sample_id2 <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster2 <- filterByExpr(spe_pseudo, group = spe_pseudo$seurat_label)
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id2==T & rowData(spe_pseudo)$high_expr_group_cluster2==T,]
spe_se <- spe_pseudo

new.cdata = merge(colData(spe_se)[,c("sample_id","brnum","condition","sex","seurat_label","nspots","chrM_ratio","pc3","slide", "seq")],
        demo[,c("sample_id","brnum","sex","age","RIN","BMI","Smoking")], sort=F)
stopifnot(identical(spe_se$pc3, new.cdata$pc3))
colData(spe_se) <- new.cdata

# correlation
var.m2 = getVarianceExplained(spe_se, variables=covars, exprs_values="logcounts")
summary(var.m2)

cor.var.m2 = cor(var.m2, method="pearson")
cor.var.m2[cor.var.m2==1] = NA

phm2 = pheatmap(cor.var.m2[ordered, ordered], angle_col=90, main="domain-CT gene var corr",
                cluster_cols=F, cluster_rows=F,
                color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100), na_col="grey50",
                breaks=seq(-1, 1, length.out=101))


ggsave(file="plots/publication/supp_covariate-selection/gene-level-variance_correlation.pdf",
	arrangeGrob(grobs=list(phm1[[4]], phm2[[4]]), ncol=1, top=NULL),
	height=5, width=3)

# binned explanation
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
	geom_bar(stat="identity", position="stack", color="black", linewidth=.1)+
	scale_fill_brewer("% gene\nvar. (bin)", palette="Greys")+
	facet_grid(cols=vars(annotation), rows=vars(covar), scales="free_x", space="free_x")+
	theme_bw()+theme(text=element_text(size=6))

ggsave(file="plots/publication/supp_covariate-selection/gene-level-variance_by-domain_binned.pdf", 
	p1, height=5, width=3)




cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
