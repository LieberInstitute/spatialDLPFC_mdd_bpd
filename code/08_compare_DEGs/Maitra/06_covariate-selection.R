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

load("processed-data/08_compare_DEGs/Maitra/sce_Maitra_pseudo_sample-low-res_norm-filt.Rdata")
colnames(colData(sce_pseudo))[15] = "ribo_ratio"

colData(sce_pseudo)$pc3_svg = reducedDim(sce_pseudo, "PCA_1663")[,"PC3"]
colData(sce_pseudo)$pc3_hvg = reducedDim(sce_pseudo, "PCA_HVG")[,"PC3"]

covars = c("Sample","Batch","Chemistry","Sequencing","condition","sex","ncells","ribo_ratio","detected","pc3_svg","pc3_hvg")

var.m = getVarianceExplained(sce_pseudo, variables=covars, exprs_values="logcounts")
cor.var.m = cor(var.m, method="pearson")
cor.var.m[cor.var.m==1] = NA

phm = pheatmap(cor.var.m, angle_col=90, fontsize=8, 
	color=colorRampPalette(RColorBrewer::brewer.pal(n=7, "RdBu")[7:1])(100),
	na_col="grey50", treeheight_row=12, treeheight_col=12,
	breaks=seq(-1, 1, length.out=101))

l1 <- lapply(levels(sce_pseudo$seurat_low.res), function(x) {
        tmp <- sce_pseudo[,sce_pseudo$seurat_low.res==x]
        var.m = as.data.frame(getVarianceExplained(tmp, variables=c("ncells","ribo_ratio","pc3_svg","pc3_hvg")))
        var.m$cluster = x
        return(var.m)
})

df1 = do.call(rbind, l1) %>% mutate(cluster=factor(cluster, levels=levels(sce_pseudo$seurat_low.res)))

df_tally <- do.call(rbind, lapply(c("ncells","ribo_ratio","pc3_svg","pc3_hvg"), function(x) {
        tmp <- df1
        tmp$plot.me = cut(tmp[[x]], breaks=c(0,1,5,10,25,50,100), include.lowest=T)

        return(group_by(tmp, plot.me, cluster) %>% tally() %>% mutate(covar=x))
})) %>% mutate(covar=factor(covar, levels=c("ncells","pc3_hvg","pc3_svg","ribo_ratio")))

p1 <- ggplot(df_tally, aes(x=cluster, y=n, fill=plot.me))+
        geom_bar(stat="identity", position="stack", color="black", linewidth=.1)+
        scale_fill_brewer("% gene\nvar. (bin)", palette="Greys")+
        facet_grid(rows=vars(covar), scales="free_x", space="free_x")+
        theme_bw()+theme(text=element_text(size=6))

ggsave(file="plots/08_compare_DEGs/Maitra/covariate-selection_gene-level-variance.png",
	arrangeGrob(phm[[4]], p1, layout_matrix=matrix(c(1,2,2), ncol=1)),
	height=8, width=4, bg="white")
cat("\nCovariate selection plots saved to: plots/08_compare_DEGs/Maitra/covariate-selection_gene-level-variance.png\n")


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
