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
cpList$smoothed.bright = c(cpList$smoothed.bright, "low UMI"="grey50")

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-with-lowUMI_sample-smoothed-n1663-k9.Rdata")
dim(spe_pseudo) 

#use the gene set from the without low UMI object b/c with low UMI cluster there are fewer genes that pass 
spe_save = spe_pseudo

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9.Rdata")
dim(spe_pseudo) 

#filter by expression before recalculating norm counts
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$smoothed_k9_1663)

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))

keep.genes = rowData(spe_pseudo)$high_expr_group_sample_id & rowData(spe_pseudo)$high_expr_group_cluster
table(keep.genes)

spe_pseudo <- spe_save[keep.genes, ]

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

p1 <- ggplot(cdata, aes(x=condition, y=nspots))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_wrap(vars(smoothed_k9_1663), ncol=7, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Number of spots per pseudobulked sample", y="nspots")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank(),
	axis.text.y=element_text(size=6))

p2 <- ggplot(cdata, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k9_1663))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k9_1663))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(smoothed_k9_1663))+
  scale_shape_manual(values=c(19,1))+ylim(0,45)+
  scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1663-k9_with-lowUMI_unfiltered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
        bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1663-k9_with-lowUMI_unfiltered_QC-metrics.png\n")

#remove MT- genes prior to calculating norm factors (this improves histogram of norm factors)
spe_pseudo = spe_pseudo[-grep("MT-", rowData(spe_pseudo)$gene_name),]
tmp = calcNormFactors(spe_pseudo)
x = cpm(tmp, log=T, prior.count=2)
stopifnot(min(x)>0)
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x


#PCA before filtering reveals one component dominated by low detected genes samples
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
n1663.ids = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% geneList$qual_genes]
length(n1663.ids)

exp.vars = c("smoothed_k9_1663","sample_id",
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

spe_pseudo <- runPCA(spe_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1663")
p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=20, variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "smoothed_k9_1663", point_alpha=1)+
        scale_color_manual("", values=cpList$smoothed.bright)
p4 <- plotPCA(spe_pseudo, dimred="PCA_1663", ncomponents=4, colour_by = "smoothed_k9_1663")+
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

bp1 <- ggplot(group_by(cdata, cond_sex, smoothed_k9_1663) %>% summarise(n_total=sum(nspots)), 
       aes(x=cond_sex, y=n_total, fill=smoothed_k9_1663))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("# spots", labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("PRECAST\n(smoothed)", values=cpList$smoothed.bright)+
  labs(title="Cluster abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(),
	legend.position="bottom")

#PCA colored by detected
p5 <- plotPCA(spe_pseudo, dimred="PCA_1663", ncomponents=4, colour_by = "detected")+
        scale_color_viridis_c("detected\ngenes")
q5 = ggplot_build(p5)
q5$data[[2]]$size = 1
q5 = ggplot_gtable(q5)

ggsave("plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1663-k9_with-lowUMI_unfiltered_PCA-1663-eval.png",
	gridExtra::grid.arrange(q1, q2, p3, q4, bp1, q5, layout_matrix=cbind(c(1,3,3,5,5),c(2,4,4,6,6))),
        bg="white", height=12, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/PRECAST_smoothed/sample-smoothed-n1663-k9_with-lowUMI_unfiltered_PCA-1663-eval.png\n")


save(spe_pseudo, file="processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-with-lowUMI_sample-smoothed-n1663-k9_norm.Rdata")
##### stop early to evaluate filters then run again
stop("Evaluate plots then continue")
#####
