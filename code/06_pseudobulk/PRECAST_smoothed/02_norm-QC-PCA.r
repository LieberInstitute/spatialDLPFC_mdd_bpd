setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(edgeR)
	library(scuttle)
	library(scater)
	#library(batchelor)
	#library(BiocParallel)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7.Rdata")
dim(spe_pseudo) #

#filter by expression before recalculating norm counts
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k7_1626)

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))

spe_keep.genes = rownames(spe_pseudo)[rowData(spe_pseudo)$high_expr_group_sample_id & rowData(spe_pseudo)$high_expr_group_cluster]
cat("\nGenes to keep from edgeR filtered PRECAST (smoothed) gene set:", length(spe_keep.genes), "\n")

#include any genes from MBv seurat label transfer
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res.Rdata")
rowData(sce_pseudo)$high_expr_group_individualID <- filterByExpr(sce_pseudo, group = sce_pseudo$individualID)
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$seurat_low.res)
sce_keep.genes = rownames(sce_pseudo)[rowData(sce_pseudo)$high_expr_group_individualID & rowData(sce_pseudo)$high_expr_group_cluster]
cat("Genes to keep from edgeR filtered SZBDMulti-seq gene set:", length(sce_keep.genes), "\n")


cat("SZBDMulti-seq genes completely missing from PRECAST (smoothed) gene set:", length(setdiff(sce_keep.genes, rownames(spe_pseudo))),"\n")
sce_keep.genes = intersect(sce_keep.genes, rownames(spe_pseudo))
cat("Revised genes to keep from edgeR filtered SZBDMulti-seq gene set:", length(sce_keep.genes), "\n")

rowData(spe_pseudo)$high_expr_SZBDMultiseq = rownames(spe_pseudo) %in% sce_keep.genes

keep.genes = union(spe_keep.genes, sce_keep.genes)
stopifnot(length(setdiff(keep.genes, rownames(spe_pseudo)))==0)

cat("\nNumber of genes to keep (union of PRECAST (smoothed) and eligible SZBDMuli-seq):", length(keep.genes), "\n")
spe_pseudo <- spe_pseudo[keep.genes, ]

#calculate QC metrics, then remove MT genes prior to normalizing

#pseudobulk qc
mt.genes = rownames(spe_pseudo)[grep("MT-",rowData(spe_pseudo)$gene_name)]
cat("\nNumber of mito genes:\n")
length(mt.genes) #13

ribo.genes = rownames(spe_pseudo)[grep("RPS|RPL",rowData(spe_pseudo)$gene_name)]
cat("\nNumber of ribo genes:\n")
length(ribo.genes) #180

spe_pseudo <- addPerCellQC(spe_pseudo, subsets = list(mito = mt.genes, ribo=ribo.genes))

#first go around i removed spots with very low ncells (<=7) but even so a few spots with v low detected genes generate a top 4 PC of low detected genes
#this PC disappeared after removing spots with <8k detected genes
#so second go we are starting with that step  
cdata= as.data.frame(colData(spe_pseudo))
#ggplot(cdata, aes(x=combined_cluster, y=detected))+
#  geom_boxplot()+scale_y_log10()+
#  geom_point(data=filter(cdata, detected<8000), color="red3")+
#  theme_bw()
p1 <- ggplot(cdata, aes(x=condition, y=nspots))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_wrap(vars(smoothed_k7_1626), ncol=6, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Number of spots per pseudobulked sample", y="nspots")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p2 <- ggplot(cdata, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k7_1626))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k7_1626))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k7_1626))+
  scale_shape_manual(values=c(19,1))+ylim(0,45)+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_unfiltered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
        bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_unfiltered_QC-metrics.png\n")

#remove MT- genes prior to calculating norm factors (this improves histogram of norm factors)
spe_pseudo = spe_pseudo[-grep("MT-", rowData(spe_pseudo)$gene_name),]
tmp = calcNormFactors(spe_pseudo)
x = cpm(tmp, log=T, prior.count=2)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x


#PCA before filtering reveals one component dominated by low detected genes samples
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_conservative_geneList.rds")
n1663.ids = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% geneList$qual_genes]
length(n1663.ids)

exp.vars = c("smoothed_k7_1626","sample_id",
             "slide","seq",
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

spe_pseudo <- runPCA(spe_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1626")
p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1626", variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1626", npcs_to_plot=20, variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(spe_pseudo, dimred="PCA_1626", ncomponents=2, colour_by = "smoothed_k7_1626", point_alpha=1)+
        scale_color_manual("", values=cpList$smoothed.bright)
p4 <- plotPCA(spe_pseudo, dimred="PCA_1626", ncomponents=4, colour_by = "smoothed_k7_1626")+
        scale_color_manual("", values=cpList$smoothed.bright)

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

#barplot of samples
cdata = as.data.frame(colData(spe_pseudo))
cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))

bp1 <- ggplot(group_by(cdata, cond_sex, smoothed_k7_1626) %>% summarise(n_total=sum(nspots)), 
       aes(x=cond_sex, y=n_total, fill=smoothed_k7_1626))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("# spots", labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("smoothed\nPRECAST\ncluster", values=cpList$smoothed.bright)+
  labs(title="Cluster abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(),
	legend.position="bottom")

#PCA colored by detected
p5 <- plotPCA(spe_pseudo, dimred="PCA_1626", ncomponents=4, colour_by = "detected")+
        scale_color_viridis_c("detected\ngenes")
q5 = ggplot_build(p5)
q5$data[[2]]$size = 1
q5 = ggplot_gtable(q5)

ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_unfiltered_PCA-1626-eval.png",
	gridExtra::grid.arrange(q1, q2, p3, q4, bp1, q5, layout_matrix=cbind(c(1,3,3,5,5),c(2,4,4,6,6))),
        bg="white", height=12, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_unfiltered_PCA-1626-eval.png\n")

save(spe_pseudo, file="processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm.Rdata")
cat("\nFiltered, normalized pseudobulk spe saved to: processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm.Rdata\n")

##### stop early to evaluate filters then run again
#stop("Evaluate plots then continue")
#####


#apply sample filter
cat("\nSample number per group after removing samples with <5k detected genes:\n")
cdata2 = cdata[cdata$detected>5000,]
table(cdata2[,c("smoothed_k7_1626","condition","sex")])

cat("\n\nSample number per group after removing L1 with <8k detected genes:\n")
remove.sample = cdata2$detected<8000 & cdata2$smoothed_k7_1626=="L1"
cdata2 = cdata2[!remove.sample,]
table(cdata2[,c("smoothed_k7_1626","condition","sex")])

cat("\nFilter out spots with low # detected genes...\n")
spe_pseudo = spe_pseudo[,spe_pseudo$detected>5000]
remove.sample = spe_pseudo$detected<8000 & spe_pseudo$smoothed_k7_1626=="L1"
spe_pseudo = spe_pseudo[,!remove.sample]
dim(spe_pseudo)

#renorm after removing samples
tmp = calcNormFactors(spe_pseudo)
x = cpm(tmp, log=T, prior.count=2)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x

#boxplots of QC metrics by condition and sex
p1 <- ggplot(cdata2, aes(x=condition, y=nspots))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+                     
  facet_wrap(vars(smoothed_k7_1626), ncol=6, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Number of spots per pseudobulked sample", y="nspots")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p2 <- ggplot(cdata2, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k7_1626))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata2, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k7_1626))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata2, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k7_1626))+
  scale_shape_manual(values=c(19,1))+ylim(0,45)+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_filtered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
	bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_filtered_QC-metrics.png\n")


#pca after filtering
spe_pseudo <- runPCA(spe_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1626")
p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1626", variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1626", npcs_to_plot=20, variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(spe_pseudo, dimred="PCA_1626", ncomponents=2, colour_by = "smoothed_k7_1626", point_alpha=1)+
	scale_color_manual("", values=cpList$smoothed.bright)
p4 <- plotPCA(spe_pseudo, dimred="PCA_1626", ncomponents=4, colour_by = "smoothed_k7_1626")+
	scale_color_manual("", values=cpList$smoothed.bright)

#retroactively modify linewidth and point size
q1 <- change_linewidth(p1)
q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

#barplot
cdata =	as.data.frame(colData(spe_pseudo))
cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))

bp1 <- ggplot(group_by(cdata, cond_sex, smoothed_k7_1626) %>% summarise(n_total=sum(nspots)),
       aes(x=cond_sex, y=n_total, fill=smoothed_k7_1626))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("#	spots",	labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("smoothed\nPRECAST\ncluster", values=cpList$smoothed.bright)+
  labs(title="Cluster abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(), 
        legend.position="bottom")

#PCA colored by detected
p5 <- plotPCA(spe_pseudo, dimred="PCA_1626", ncomponents=4, colour_by = "detected")+
        scale_color_viridis_c("detected\ngenes")
q5 = ggplot_build(p5)
q5$data[[2]]$size = 1
q5 = ggplot_gtable(q5)

ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_filtered_PCA-1626-eval.png", 
	gridExtra::grid.arrange(q1, q2, p3, q4, bp1, q5, layout_matrix=cbind(c(1,3,3,5,5),c(2,4,4,6,6))),
	bg="white", height=12, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_filtered_PCA-1626-eval.png\n")

#look for correlation between important experimental design variables
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
hmp = pheatmap::pheatmap(cor.var.m,
	annotation_col = col_annot, annotation_colors = ann_colors,
	annotation_names_col=FALSE, annotation_legend=T)
ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1626-k7_conservative_variance-explained_experimental-design_heatmap.png",
       hmp[[4]], bg="white", height=7, width=9, units="in")


#save csv with avg expr and quantile just like i did with nnSVG filt genes
avg1 = rowMeans(logcounts(spe_pseudo))
q.decile = quantile(avg1, prob=seq(0,1,.1))
avg.expr = cbind.data.frame("gene_name"=rowData(spe_pseudo)[names(avg1),"gene_name"],
                 "avg_expr"=avg1,
                 "decile" = cut(avg1, breaks=c(0,q.decile[2:11]), labels=F))
write.csv(avg.expr, "processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1626-k7_conservative_filtered-genes_avg-logcounts.csv", row.names=T)
cat("\nAverage expression of",paste0(length(keep.genes)),"genes after filtering saved to: processed-data/06_pseudobulk/PRECAST_smoothed/pseudobulk-sample-smoothed-n1626-k7_conservative_filtered-genes_avg-logcounts.csv\n")

save(spe_pseudo, file="processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm-filt.Rdata")
cat("\nFiltered, normalized pseudobulk spe saved to: processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm-filt.Rdata\n")

#update spe tracker
write(c(paste("********** QC filtered and normalized pseudobulked spe on",format(Sys.time())),
        "********** Old file location: processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7.Rdata",
        "********** New file location: processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_conservative_pseudo_sample-smoothed-n1626-k7_norm-filt.Rdata",
        "********** Source code: code/06_pseudobulk/PRECAST_smoothed/02_norm-QC-PCA.r",
        "**********","**********","**********"), "spe_tracker_current.txt", append=TRUE)

cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
