setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(sva)
	library(edgeR)
	library(scater)
	library(ggplot2)
})
set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata")

#compute adjusted counts (batch correction)
#want to preserve biological variation belonging to diagnosis, sex, and cluster, and want to remove variation belonging to sample_id
#according to page 10: https://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
#need to make contrast matrix of diagnosis, sex, and cluster, as these will be treated as covariates (# covar levels = n-1; therefore set condition first for NTC)
cat("\ncovars included (to preserve): condition, sex, precast_k9_1663\n")
m1 = model.matrix(~condition + sex + precast_k9_1663, data=colData(spe_pseudo))
#covar_matrix = as.matrix(m1[,-1])
cat("\nbatch variable to remove: slide\n")
adjusted_counts = ComBat_seq(counts(spe_pseudo), batch=spe_pseudo$slide, group=NULL, covar_mod=m1)
dimnames(adjusted_counts) <- dimnames(spe_pseudo)
assay(spe_pseudo, "counts") <- NULL
assay(spe_pseudo, "logcounts") <- NULL
assay(spe_pseudo, "adjusted_counts") <- adjusted_counts

#quick save in case of code error downstream
save(spe_pseudo, file="processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata")
#load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata")
adjusted_counts = assay(spe_pseudo, "adjusted_counts")

#normalize adjusted counts
dge = DGEList(counts=adjusted_counts)
x <- cpm(calcNormFactors(dge), log = TRUE, prior.count = 2)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x
rm(x)

#pca
reducedDim(spe_pseudo, "PCA_1663") <- NULL
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
sapply(geneList, length)
n1663.ids = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% geneList$qual_genes]
length(n1663.ids)

exp.vars = c("precast_k9_1663","sample_id",
             "slide","round",
             "condition","sex",
             "age","PMI","RIN",
             "sum","detected","nspots")
#colors based on RColorBrewer::brewer.pal(n=length(exp.vars), "Set3")
exp.vars.colors = c("#FB8072", "#80B1D3",
                    "#BC80BD", "#B3DE69",
                    "#8DD3C7", "#FDB462",
                    "#FCCDE5", "#BEBADA", "#FFFFB3",
                    "#CCEBC5", "#D9D9D9","black")
names(exp.vars.colors) = exp.vars
cpList <- readRDS("plots/colorPalettes.rds")

set.seed(123)
spe_pseudo <- runPCA(spe_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1663")

p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=20, variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "precast_k9_1663", point_alpha=1)+
	scale_color_manual("", values=cpList$earthy.pal2)
p4 <- plotPCA(spe_pseudo, dimred="PCA_1663", ncomponents=4, colour_by = "precast_k9_1663")+
	scale_color_manual("", values=cpList$earthy.pal2)

#retroactively modify linewidth and point size
change_linewidth = function(.plot) {
        m = ggplot_build(.plot)
        m$data[[1]]$size = 2
        m$data[[2]]$linewidth = 1
        m1 = ggplot_gtable(m)
        return(m1)
}
q1 <- change_linewidth(p1)
q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

ggsave("plots/06_pseudobulk/sample-n1663-k9_combat-seq_PCA-1663-eval.png", gridExtra::grid.arrange(q1, q2, p3, q4, layout_matrix=cbind(c(1,3,3),c(2,4,4))),
	bg="white", height=8, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/sample-n1663-k9_combat-seq_PCA-1663-eval.png\n")

#gene-level variance explained correlation
var.m = getVarianceExplained(spe_pseudo, variables=exp.vars)
#par(mfrow=c(4,3))
#for(i in colnames(var.m)) {
#  qqnorm(var.m[,i], pch = 1, frame = FALSE, main=i)
#  qqline(var.m[,i], col = "steelblue", lwd = 2) 
#}
#decidedly not normal distribution
cor.var.m = cor(var.m, method="spearman")
col_annot = data.frame(colMeans(var.m))
colnames(col_annot) = "percVar"
ann_colors = list(
  percVar = colorRampPalette(c("white", "purple3", "black"), bias=1)(10)
)

hmp = pheatmap::pheatmap(cor.var.m, annotation_col = col_annot, annotation_colors = ann_colors,
                         annotation_names_col=FALSE, annotation_legend=T)
ggsave("plots/06_pseudobulk/sample-n1663-k9_variance-explained_experimental-design_combat-seq_heatmap.png",
       hmp[[4]], bg="white", height=7, width=9, units="in")

save(spe_pseudo, file="processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata")
cat("\nSpe with combat-seq adjusted counts and reduced dims saved to: processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata\n")

#update spe tracker
write(c(paste("************* Combat-Seq correction of pseudobulked spe on",format(Sys.time())),
        "************* Old file location: processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt.Rdata",
        "************* New file location: processed-data/06_pseudobulk/spe_n119_pseudo_sample-n1663-k9_norm-filt-combat-seq.Rdata",
        "************* Source code: code/06_pseudobulk/03_combat-seq.r",
        "*************","*************","*************"), "spe_tracker_current.txt", append=TRUE)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
