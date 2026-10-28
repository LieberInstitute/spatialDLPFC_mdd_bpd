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

load("processed-data/08_compare_DEGs/SZBDMulti-seq/sce_SZBD_pseudo_indivID-channel-low-res_norm-filt.Rdata")
colnames(colData(sce_pseudo))[grep("_percent", colnames(colData(sce_pseudo)))] = "ribo_ratio"

colData(sce_pseudo)$pc3_svg = reducedDim(sce_pseudo, "PCA_1663")[,"PC3"]
colData(sce_pseudo)$pc7_svg = reducedDim(sce_pseudo, "PCA_1663")[,"PC7"]
#colData(sce_pseudo)$pc9_svg = reducedDim(sce_pseudo, "PCA_1663")[,"PC9"]
colData(sce_pseudo)$pc8_hvg = reducedDim(sce_pseudo, "PCA_HVG")[,"PC8"]
#colData(sce_pseudo)$pc10_hvg = reducedDim(sce_pseudo, "PCA_HVG")[,"PC10"]
#colData(sce_pseudo)$pc11_hvg = reducedDim(sce_pseudo, "PCA_HVG")[,"PC11"]

# demographic data anova
demo = distinct(as.data.frame(colData(sce_pseudo)), individualID, sex, condition, Age_death, PMI)

#lm
cat("\n\nLinear regression of dx*sex on donor covars...")
cat("\n>>> Age...\n")
m1 = lm(Age_death ~ condition*sex, data=demo)
anova(m1)

cat("\n>>> PMI...\n")
m2 = lm(PMI ~ condition*sex, data=demo)
anova(m2)

# now plotting
p1 <- ggplot(demo, aes(x=condition, y=Age_death))+geom_boxplot(outliers=F)+
	ggbeeswarm::geom_beeswarm()+
	facet_grid(cols=vars(sex))+
	ylim(0,100)+ggtitle("Age ~ condition = sig. with ANOVA")

p2 <- ggplot(demo, aes(x=condition, y=PMI))+geom_boxplot(outliers=F)+
        ggbeeswarm::geom_beeswarm()+
        facet_grid(cols=vars(sex))+
	ylim(0,42)+ggtitle("PMI is NS by ANOVA")


covars = c("individualID","Channel","batch","Age_death","PMI","condition","sex","ncells","ribo_ratio","detected","pc8_hvg","pc3_svg","pc7_svg")

var.m = getVarianceExplained(sce_pseudo, variables=covars, exprs_values="logcounts")

cat("\nTop 5 genes with highest % variance explained by each covar:\n")
sapply(colnames(var.m), function(x) rowData(sce_pseudo)[names(sort(var.m[,x], decreasing=T)[1:5]),"gene_name"])

cor.var.m = cor(var.m, method="pearson")
cor.var.m[cor.var.m==1] = NA

phm = pheatmap(cor.var.m, angle_col=90, fontsize=8, 
	color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100),
	na_col="grey50", treeheight_row=12, treeheight_col=0,
	breaks=seq(-1, 1, length.out=101))

l1 <- lapply(levels(sce_pseudo$seurat_low.res), function(x) {
        tmp <- sce_pseudo[,sce_pseudo$seurat_low.res==x]
        var.m = as.data.frame(getVarianceExplained(tmp, variables=c("ncells","Age_death","batch","pc8_hvg","pc3_svg","pc7_svg")))
        var.m$cluster = x
        return(var.m)
})

df1 = do.call(rbind, l1) %>% mutate(cluster=factor(cluster, levels=levels(sce_pseudo$seurat_low.res)))

df_tally <- do.call(rbind, lapply(c("ncells","Age_death","batch","pc8_hvg","pc3_svg","pc7_svg"), function(x) {
        tmp <- df1
        tmp$plot.me = cut(tmp[[x]], breaks=c(0,1,5,10,25,50,100), include.lowest=T)

        return(group_by(tmp, plot.me, cluster) %>% tally() %>% mutate(covar=x))
})) %>% mutate(covar=factor(covar, levels=c("ncells","Age_death","batch","pc8_hvg","pc3_svg","pc7_svg")))

p3 <- ggplot(df_tally, aes(x=cluster, y=n, fill=plot.me))+
        geom_bar(stat="identity", position="stack", color="black", linewidth=.1)+
        scale_fill_brewer("% gene\nvar. (bin)", palette="Greys")+
        facet_grid(rows=vars(covar), scales="free_x", space="free_x")+
        theme_bw()+theme(text=element_text(size=6))

ggsave(file="plots/08_compare_DEGs/SZBDMulti-seq/covariate-selection_channel_demo_gene-level-variance.png",
	arrangeGrob(p1, p2, phm[[4]], p3, layout_matrix=cbind(c(1,1,2,2),c(3,3,4,4))),
	height=8, width=8, bg="white")
cat("\nCovariate selection plots saved to: plots/08_compare_DEGs/SZBDMulti-seq/covariate-selection_channel_demo_gene-level-variance.png\n")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
