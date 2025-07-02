setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(edgeR)
	library(scuttle)
	library(scater)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")
low.res.pal = c("Astro"="#cfa45c","Micro/Vasc"="#911223",
                "Inhb"="#9377AC",
                "L2"="#5D9940","L3"="#5095CD",
                "L4"="#85A0A0",
                "L5"="#ddc94e","L6"="#E45C5F",
                "Oligo"="#D1C4B0")
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res.Rdata")
dim(sce_pseudo) #

sce_pseudo$Age_death = as.numeric(as.character(sce_pseudo$Age_death))

#filter by expression before recalculating norm counts
rowData(sce_pseudo)$high_expr_group_individualID <- filterByExpr(sce_pseudo, group = sce_pseudo$individualID)
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$seurat_low.res)

with(rowData(sce_pseudo), table(high_expr_group_individualID, high_expr_group_cluster))

keep.genes = rowData(sce_pseudo)$high_expr_group_individualID & rowData(sce_pseudo)$high_expr_group_cluster
table(keep.genes)
sce_pseudo <- sce_pseudo[keep.genes, ]

#calculate QC metrics, then remove MT genes prior to normalizing

#pseudobulk qc
mt.genes = rownames(sce_pseudo)[grep("MT-",rowData(sce_pseudo)$gene_name)]
cat("\nNumber of mito genes:\n")
length(mt.genes) #13

ribo.genes = rownames(sce_pseudo)[grep("RPS|RPL",rowData(sce_pseudo)$gene_name)]
cat("\nNumber of ribo genes:\n")
length(ribo.genes) #180

sce_pseudo <- addPerCellQC(sce_pseudo, subsets = list(mito = mt.genes, ribo=ribo.genes))

#first go around i removed spots with very low ncells (<=7) but even so a few spots with v low detected genes generate a top 4 PC of low detected genes
#this PC disappeared after removing spots with <8k detected genes
#so second go we are starting with that step  
cdata= as.data.frame(colData(sce_pseudo)) %>% 
	mutate(sex=factor(Biological_Sex, levels=c("female","male"), labels=c("F","M")),
		condition=factor(Disorder, levels=c("control"), labels=c("con"))
	)
#ggplot(cdata, aes(x=combined_cluster, y=detected))+
#  geom_boxplot()+scale_y_log10()+
#  geom_point(data=filter(cdata, detected<8000), color="red3")+
#  theme_bw()
p1 <- ggplot(cdata, aes(x=condition, y=ncells))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_wrap(vars(seurat_low.res), ncol=8, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  theme_bw()+labs(title="Number of nuclei per pseudobulked sample", y="ncells")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p2 <- ggplot(cdata, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+#ylim(0,45)+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/06_pseudobulk/SZBDMulti-seq/control-low-res_unfiltered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
        bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/SZBDMulti-seq/control-low-res_unfiltered_QC-metrics.png\n")

#remove MT- genes prior to calculating norm factors (this improves histogram of norm factors)
### MT genes already removed from source dataset
#sce_pseudo = sce_pseudo[-grep("MT-", rowData(sce_pseudo)$gene_name),]
tmp = calcNormFactors(sce_pseudo)
x = cpm(tmp, log=T, prior.count=2)
warning("After norm. with prior.count=2, min counts still <0.\nTried increasing up to prior.count=6 and still <0.")
#stopifnot(min(x)>0)
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x


#PCA before filtering
varfeat.id = rownames(sce_pseudo)[rowData(sce_pseudo)$VarFeat_MBv.filtered]
cat("\nUsed VariableFeatures for PCA...\n")
length(varfeat.id)

exp.vars = c("seurat_low.res","individualID",
             #"slide","seq",
             "Biological_Sex",
             "Age_death",#"PMI","RIN",
             "sum","detected","ncells")
#colors based on RColorBrewer::brewer.pal(n=length(exp.vars), "Set3")
exp.vars.colors = c("#FB8072", "#80B1D3",
                    #"#BC80BD", "#B3DE69",
                    "#FDB462",
                    "#FCCDE5", #"#BEBADA", "#FFFFB3",
                    "#CCEBC5", "#D9D9D9","black")
names(exp.vars.colors) = exp.vars


sce_pseudo <- runPCA(sce_pseudo, subset_row=varfeat.id, exprs_values="logcounts", name="PCA_VarFeat")
p1 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_VarFeat", variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_VarFeat", npcs_to_plot=20, variables=exp.vars)+
        scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(sce_pseudo, dimred="PCA_VarFeat", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
        scale_color_manual("", values=low.res.pal)
p4 <- plotPCA(sce_pseudo, dimred="PCA_VarFeat", ncomponents=4, colour_by = "seurat_low.res")+
        scale_color_manual("", values=low.res.pal)

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
#cdata = as.data.frame(colData(sce_pseudo))
#cdata$cond_sex = factor(paste(cdata$condition, cdata$sex), levels=c("NTC M","NTC F","MDD M","MDD F","BPD M","BPD F"))

bp1 <- ggplot(group_by(cdata, sex, seurat_low.res) %>% summarise(n_total=sum(ncells)), 
       aes(x=sex, y=n_total, fill=seurat_low.res))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("# nuclei", labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("Seurat\nlabels", values=low.res.pal)+
  labs(title="Cell type abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(),
	legend.position="bottom")

#PCA colored by detected
p5 <- plotPCA(sce_pseudo, dimred="PCA_VarFeat", ncomponents=4, colour_by = "detected")+
        scale_color_viridis_c()
q5 = ggplot_build(p5)
q5$data[[2]]$size = 1
q5 = ggplot_gtable(q5)

ggsave("plots/06_pseudobulk/SZBDMulti-seq/control-low-res_unfiltered_PCA-VarFeat-eval.png",
	gridExtra::grid.arrange(q1, q2, p3, q4, bp1, q5, layout_matrix=cbind(c(1,3,3,5,5),c(2,4,4,6,6))),
        bg="white", height=12, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/SZBDMulti-seq/control-low-res_unfiltered_PCA-VarFeat-eval.png\n")


##### stop early to evaluate filters then run again
#stop("Evaluate plots then continue")
#####

#apply sample filter
cat("\nSample number per group after removing samples with <5k detected genes:\n")
cdata2 = cdata[cdata$detected>5000,]
table(cdata2[,c("seurat_low.res","sex")])

cat("\nFilter out spots with low # detected genes...\n")
sce_pseudo = sce_pseudo[,sce_pseudo$detected>5000]
dim(sce_pseudo)

#renorm after removing samples
tmp = calcNormFactors(sce_pseudo)
x = cpm(tmp, log=T, prior.count=2)
#stopifnot(min(x)>0)
dimnames(x) <- dimnames(sce_pseudo)
logcounts(sce_pseudo) <- x

#boxplots of QC metrics by condition and sex
p1 <- ggplot(cdata2, aes(x=condition, y=ncells))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+                     
  facet_wrap(vars(seurat_low.res), ncol=8, scales="free_y")+
  scale_shape_manual(values=c(19,1))+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Number of nuclei per pseudobulked sample", y="ncells")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p2 <- ggplot(cdata2, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p3 <- ggplot(cdata2, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

p4 <- ggplot(cdata2, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(shape=sex))+
  geom_boxplot(alpha=.4, linewidth=.5, outliers=F)+
  facet_grid(cols=vars(seurat_low.res))+
  scale_shape_manual(values=c(19,1))+#ylim(0,45)+
  #scale_color_manual("diagnosis", values=cpList$dx.pal)+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA), panel.grid.minor=element_blank())

ggsave("plots/06_pseudobulk/SZBDMulti-seq/control-low-res_filtered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
	bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/SZBDMulti-seq/control-low-res_filtered_QC-metrics.png\n")


#pca after filtering
sce_pseudo <- runPCA(sce_pseudo, subset_row=varfeat.id, exprs_values="logcounts", name="PCA_VarFeat")
p1 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_VarFeat", variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(sce_pseudo, dimred="PCA_VarFeat", npcs_to_plot=20, variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(sce_pseudo, dimred="PCA_VarFeat", ncomponents=2, colour_by = "seurat_low.res", point_alpha=1)+
	scale_color_manual("", values=low.res.pal)
p4 <- plotPCA(sce_pseudo, dimred="PCA_VarFeat", ncomponents=4, colour_by = "seurat_low.res")+
	scale_color_manual("", values=low.res.pal)

#retroactively modify linewidth and point size
q1 <- change_linewidth(p1)
q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

#barplot
cdata =	as.data.frame(colData(sce_pseudo))
cdata$sex = factor(cdata$Biological_Sex, levels=c("female","male"), labels=c("F","M"))

bp1 <- ggplot(group_by(cdata, sex, seurat_low.res) %>% summarise(n_total=sum(ncells)),
       aes(x=sex, y=n_total, fill=seurat_low.res))+
  geom_bar(stat="identity", position="stack", width=.7)+
  scale_y_continuous("#	nuclei", labels=function(x) paste0(x/1000,"k"))+
  scale_fill_manual("Seurat\nlabels", values=low.res.pal)+
  labs(title="Cell type abundance")+
  theme_bw()+theme(text=element_text(size=10), axis.title.x=element_blank(), 
        legend.position="bottom")

#PCA colored by detected
p5 <- plotPCA(sce_pseudo, dimred="PCA_VarFeat", ncomponents=4, colour_by = "detected")+
        scale_color_viridis_c()
q5 = ggplot_build(p5)
q5$data[[2]]$size = 1
q5 = ggplot_gtable(q5)

ggsave("plots/06_pseudobulk/SZBDMulti-seq/control-low-res_filtered_PCA-VarFeat-eval.png", 
	gridExtra::grid.arrange(q1, q2, p3, q4, bp1, q5, layout_matrix=cbind(c(1,3,3,5,5),c(2,4,4,6,6))),
	bg="white", height=12, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/SZBDMulti-seq/control-low-res_filtered_PCA-VarFeat-eval.png\n")

#look for correlation between important experimental design variables
var.m = getVarianceExplained(sce_pseudo, variables=exp.vars)
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
ggsave("plots/06_pseudobulk/SZBDMulti-seq/control-low-res_variance-explained_experimental-design_heatmap.png",
       hmp[[4]], bg="white", height=7, width=9, units="in")


#save csv with avg expr and quantile just like i did with nnSVG filt genes
avg1 = rowMeans(logcounts(sce_pseudo))
q.decile = quantile(avg1, prob=seq(0,1,.1))
avg.expr = cbind.data.frame("gene_name"=rowData(sce_pseudo)[names(avg1),"gene_name"],
                 "avg_expr"=avg1,
                 "decile" = cut(avg1, breaks=c(0,q.decile[2:11]), labels=F))
write.csv(avg.expr, "processed-data/06_pseudobulk/SZBDMulti-seq/pseudobulk-control-low-res_filtered-genes_avg-logcounts.csv", row.names=T)
cat("\nAverage expression of",paste0(length(keep.genes)),"genes after filtering saved to: processed-data/06_pseudobulk/SZBDMulti-seq/pseudobulk-control-low-res_filtered-genes_avg-logcounts.csv\n")

save(sce_pseudo, file="processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res_norm-filt.Rdata")
cat("\nFiltered, normalized pseudobulk sce saved to: processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res_norm-filt.Rdata\n")


cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
