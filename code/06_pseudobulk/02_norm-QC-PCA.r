setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DelayedArray)
	library(edgeR)
	library(scuttle)
	library(scater)
	#library(batchelor)
	#library(BiocParallel)
	#library(dplyr)
	library(ggplot2)
})
set.seed(123)

load("processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus.Rdata")
dim(spe_pseudo) #28965   952

#filter by expression before recalculating norm counts
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$combined_cluster)

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#..........................high_expr_group_cluster
#high_expr_group_sample_id FALSE  TRUE
#....................FALSE 13445     0
#....................TRUE   3123 12397

keep.genes = rowData(spe_pseudo)$high_expr_group_sample_id & rowData(spe_pseudo)$high_expr_group_cluster
table(keep.genes)
spe_pseudo <- spe_pseudo[keep.genes, ]

x <- cpm(calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
dimnames(x) <- dimnames(spe_pseudo)
logcounts(spe_pseudo) <- x
rm(x)

#save csv with avg expr and quantile just like i did with nnSVG filt genes
avg1 = rowMeans(logcounts(spe_pseudo))
q.decile = quantile(avg1, prob=seq(0,1,.1))
avg.expr = cbind.data.frame("gene_name"=rowData(spe_pseudo)[names(avg1),"gene_name"],
                 "avg_expr"=avg1,
                 "decile" = cut(avg1, breaks=c(0,q.decile[2:11]), labels=F))
write.csv(avg.expr, "processed-data/06_pseudobulk/filtered-genes_avg-logcounts.csv", row.names=T)
cat("\nAverage expression of",paste0(length(keep.genes)),"genes after filtering saved to: processed-data/06_pseudobulk/filtered-genes_avg-logcounts.csv\n")

################################################
################ SANITY CHECK OF NORM. VALUES
#
#spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
#clusters = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
#identical(rownames(colData(spe)), rownames(clusters))
#spe$combined_cluster= factor(clusters$combined_cluster, levels=c("Vasc","L1","L2","L3","GABA","L5","L6","WM"))
#
#spe_sub = spe[,spe$sample_id=="V13B23-309_A1" & !is.na(spe$combined_cluster)]
#spe_sub$mobp_log = logcounts(spe_sub)[rowData(spe_sub)$gene_name=="MOBP",]
#p = make_escheR(spe_sub) |> add_ground(var="combined_cluster") |> add_fill(var="mobp_log")
#p+scale_color_manual(values=precast.colorList[["n1663_k9"]][["colors"]][1:8])+
#  scale_fill_gradient(low="white",high="black")
#
#spe_sub$mobp_raw = counts(spe_sub)[rowData(spe_sub)$gene_name=="MOBP",]
#group_by(as.data.frame(colData(spe_sub)), combined_cluster) %>% 
#  summarise(nspots=n(), sum_umi=sum(sum_umi), sum_mobp_raw=sum(mobp_raw),
#            mobp_raw.over.sum= sum_mobp_raw/sum_umi)
#
#spe_pseudo_sub = spe_pseudo[,spe_pseudo$sample_id=="V13B23-309_A1"]
#spe_pseudo_sub$mobp_raw = counts(spe_pseudo_sub)[rowData(spe_pseudo_sub)$gene_name=="MOBP",]
#spe_pseudo_sub$mobp_log = logcounts(spe_pseudo_sub)[rowData(spe_pseudo_sub)$gene_name=="MOBP",]
#colData(spe_pseudo_sub)[spe_pseudo_sub$sample_id=="V13B23-309_A1",c("combined_cluster", "ncells", 
#                                                                    "sum", "mobp_raw", "mobp_log")]
#
#so the raw mobp are identical (obviously) and the sums aren't really that different
#x$samples[x$samples$sample_id=="V13B23-309_A1",c("combined_cluster", "ncells", "sum")]
#what messing with me is the log2 function on the cpm
#x2 = cpm(counts(spe_pseudo_sub))
#spe_pseudo_sub$mobp_cpm = x2[rowData(spe_pseudo_sub)$gene_name=="MOBP",]
#spe_pseudo_sub$mobp_lg2p1cpm = log2(x2[rowData(spe_pseudo_sub)$gene_name=="MOBP",]+1)
#
#p1 <- ggplot(as.data.frame(colData(spe_pseudo_sub)), aes(x=combined_cluster, y=mobp_cpm))+
#  geom_bar(stat="identity")+labs(title="CPM", y="MOBP")+theme_bw()
#
#p2 <- ggplot(as.data.frame(colData(spe_pseudo_sub)), aes(x=combined_cluster, y=mobp_cpm))+
#  geom_bar(stat="identity")+labs(title="CPM, y axis = log2 scale", y="MOBP")+theme_bw()+
#  scale_y_continuous(trans="log2")
#
#p3 <-  ggplot(as.data.frame(colData(spe_pseudo_sub)), aes(x=combined_cluster, y=mobp_lg2p1cpm))+
#  geom_bar(stat="identity")+labs(title="log2(CPM+1)", y="MOBP")+theme_bw()
#
#p4 <-  ggplot(as.data.frame(colData(spe_pseudo_sub)), aes(x=combined_cluster, y=mobp_log))+
#  geom_bar(stat="identity")+labs(title="calcNorm factors + cpm(log=T, prior.count=1)", y="MOBP")+theme_bw()
#
#gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1)
################################################


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
p1 <- ggplot(cdata, aes(x=condition, y=ncells))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Number of spots per pseudobulked sample", y="ncells")+
  theme(strip.background=element_rect(fill=NA, color=NA))

p2 <- ggplot(cdata, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA))

p3 <- ggplot(cdata, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA))

p4 <- ggplot(cdata, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_shape_manual(values=c(19,1))+ylim(0,45)+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA))

ggsave("plots/06_pseudobulk/unfiltered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
        bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/unfiltered_QC-metrics.png\n")

#PCA before filtering reveals one component dominated by low detected genes samples
source("code/05_clustering/PRECAST/PRECAST_colorLists.r")
geneList <- readRDS("processed-data/04_feature_selection/nnSVG-eval_geneList.rds")
n1663.ids = rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name %in% geneList$qual_genes]
length(n1663.ids)

exp.vars = c("combined_cluster","sample_id",
             "slide","round",
             "condition","sex",
             "age","PMI","RIN",
             "sum","detected","ncells")
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
p3 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "combined_cluster", point_alpha=1)+
        scale_color_manual("", values=precast.colorList[["n1663_k9"]][["colors"]][1:8])
p4 <- plotPCA(spe_pseudo, dimred="PCA_1663", ncomponents=4, colour_by = "combined_cluster")+
        scale_color_manual("", values=precast.colorList[["n1663_k9"]][["colors"]][1:8])

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

ggsave("plots/06_pseudobulk/PCA-1663_unfiltered_eval.png", gridExtra::grid.arrange(q1, q2, p3, q4, layout_matrix=cbind(c(1,3,3),c(2,4,4))),
        bg="white", height=8, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/PCA-1663_unfiltered_eval.png\n")


cdata2 = cdata[cdata$detected>8000,]
table(cdata2[,c("combined_cluster","condition","sex")])
#FEMALE
#.................condition
#combined_cluster NTC MDD BPD
#............Vasc  20  20  20
#............L1    20  20  20
#............L2    20  20  20
#............L3    20  20  20
#............GABA  19  19  20
#............L5    20  20  20
#............L6    20  19  20
#............WM    20  20  20

#MALE
#.................condition
#combined_cluster NTC MDD BPD
#............Vasc  20  19  20
#............L1    19  19  20
#............L2    20  19  20
#............L3    20  19  20
#............GABA  18  17  19
#............L5    20  19  20
#............L6    20  19  20
#............WM    20  18  20

cat("\nFilter out spots with low # detected genes...\n")
spe_pseudo = spe_pseudo[,spe_pseudo$detected>8000]
dim(spe_pseudo)

#boxplots of QC metrics by condition and sex
p1 <- ggplot(cdata2, aes(x=condition, y=ncells))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+                     
  facet_grid(cols=vars(combined_cluster))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Number of spots per pseudobulked sample", y="ncells")+
  theme(strip.background=element_rect(fill=NA, color=NA))

p2 <- ggplot(cdata2, aes(x=condition, y=sum))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_y_log10()+scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Library size per pseudobulked sample", y="sum UMI (log10 scale)")+
  theme(strip.background=element_rect(fill=NA, color=NA))

p3 <- ggplot(cdata2, aes(x=condition, y=detected))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Unique detected genes per pseudobulked sample", y="detected")+
  theme(strip.background=element_rect(fill=NA, color=NA))

p4 <- ggplot(cdata2, aes(x=condition, y=subsets_mito_percent))+
  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
  facet_grid(cols=vars(combined_cluster))+
  scale_shape_manual(values=c(19,1))+ylim(0,45)+
  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
  theme_bw()+labs(title="Fraction of chrM reads per pseudobulked sample", y="subsets_mito_percent")+
  theme(strip.background=element_rect(fill=NA, color=NA))

#p4 <- ggplot(cdata, aes(x=condition, y=subsets_ribo_percent))+
#  ggbeeswarm::geom_quasirandom(aes(color=condition, shape=sex))+
#  geom_boxplot(color="grey50", alpha=.5, linewidth=1, outliers=F)+
#  facet_grid(cols=vars(combined_cluster))+
#  scale_shape_manual(values=c(19,1))+ylim(0,13)+
#  scale_color_manual(values=c("black","#9e771b","#1b9e77"))+
#  theme_bw()+labs(title="Fraction of RPS|RPL reads per pseudobulked sample", y="subsets_ribo_percent")+
#  theme(strip.background=element_rect(fill=NA, color=NA))

ggsave("plots/06_pseudobulk/filtered_QC-metrics.png", gridExtra::grid.arrange(p1, p2, p3, p4, ncol=1),
	bg="white", width=12, height=12, units="in")
cat("\nQC plots saved to: plots/06_pseudobulk/filtered_QC-metrics.png\n")


#pca after filtering
spe_pseudo <- runPCA(spe_pseudo, subset_row=n1663.ids, exprs_values="logcounts", name="PCA_1663")
p1 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p2 <- plotExplanatoryPCs(spe_pseudo, dimred="PCA_1663", npcs_to_plot=20, variables=exp.vars)+
	scale_y_continuous()+scale_color_manual("", values=exp.vars.colors)
p3 <- plotReducedDim(spe_pseudo, dimred="PCA_1663", ncomponents=2, colour_by = "combined_cluster", point_alpha=1)+
	scale_color_manual("", values=precast.colorList[["n1663_k9"]][["colors"]][1:8])
p4 <- plotPCA(spe_pseudo, dimred="PCA_1663", ncomponents=4, colour_by = "combined_cluster")+
	scale_color_manual("", values=precast.colorList[["n1663_k9"]][["colors"]][1:8])

#retroactively modify linewidth and point size
q1 <- change_linewidth(p1)
q2 <- change_linewidth(p2)
q4 = ggplot_build(p4)
q4$data[[2]]$size = 1
q4 = ggplot_gtable(q4)

ggsave("plots/06_pseudobulk/PCA-1663_filtered_eval.png", gridExtra::grid.arrange(q1, q2, p3, q4, layout_matrix=cbind(c(1,3,3),c(2,4,4))),
	bg="white", height=8, width=12, units="in")
cat("\nPCA eval plots saved to: plots/06_pseudobulk/PCA-1663_filtered_eval.png\n")

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
ggsave("plots/06_pseudobulk/variance-explained_experimental-design_heatmap.png",
       hmp[[4]], bg="white", height=7, width=9, units="in")

#batch correction
#samp.data = distinct(as.data.frame(colData(spe_pseudo)[,c("sample_id","brnum","condition","sex","round","age","PMI","RIN")]))
#ntcM = filter(samp.data, condition=="NTC", sex=="M") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
#ntcF = filter(samp.data, condition=="NTC", sex=="F") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
#mddM = filter(samp.data, condition=="MDD", sex=="M") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
#mddF = filter(samp.data, condition=="MDD", sex=="F") %>% mutate(round=factor(round, levels=c("r2","r2_1","r1","r3"))) %>% arrange(round) %>% pull(brnum)
#bpdM = filter(samp.data, condition=="BPD", sex=="M") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
#bpdF = filter(samp.data, condition=="BPD", sex=="F") %>% mutate(round=factor(round, levels=c("r2","r1","r3"))) %>% arrange(round) %>% pull(brnum)
#m.order = list(ntcM, ntcF, mddM, mddF, bpdM, bpdF)
#sapply(m.order, length)
#set.seed(123)
#mnn <- reducedMNN(reducedDim(spe_pseudo,'PCA_1663'), batch=spe_pseudo$brnum, k=5, merge.order=m.order,
#                  BPPARAM=SerialParam())
#
##check merge order
#lost.var = mnn@metadata$merge.info$lost.var
#
#m.order2 = c(100:119,80:99,61:79,41:60,21:40,1:20)
#names(m.order2) = c(ntcM, ntcF, mddM, mddF, bpdM, bpdF)
#lost.var2 = lost.var[,names(m.order2)[m.order2]]
#rownames(lost.var2) = 1:118
#
#samp.data = distinct(as.data.frame(colData(spe_pseudo)[,c("sample_id","brnum","condition","sex","round","age","PMI","RIN")]))
#col_annot = samp.data[,c("condition","sex")]
#rownames(col_annot) = samp.data$brnum
#annot_colors = list("condition"=c(NTC="black",MDD="#9e771b",BPD="#1b9e77"), "sex"=c(M="white",F="black"))
#
#row_annot = cbind.data.frame("batch.size"=mnn@metadata$merge.info[,"batch.size"])
#hmp = pheatmap::pheatmap(lost.var2, cluster_rows=F, cluster_cols=F, annotation_row = row_annot, show_rownames = F,
#                   annotation_col=col_annot, annotation_colors=annot_colors) 
#ggsave("plots/06_pseudobulk/MNN-1663_merge-order.png", hmp[[4]], bg="white", width=9, height=9, units="in")
#cat("\nMNN merge order heatmap saved to: plots/06_pseudobulk/MNN-1663_merge-order.png\n")
#
##check impact of sample_id correction
#reducedDim(spe_pseudo, "MNN_1663") <- mnn$corrected
#p1 <- plotExplanatoryPCs(spe_pseudo, dimred="MNN_1663", variables=exp.vars)+
#	scale_y_continuous()+scale_color_manual(values=exp.vars.colors)
#p2 <- plotExplanatoryPCs(spe_pseudo, dimred="MNN_1663", npcs_to_plot=20, variables=exp.vars)+
#	scale_y_continuous()+scale_color_manual(values=exp.vars.colors)
#p3 <- plotReducedDim(spe_pseudo, dimred="MNN_1663", ncomponents=2, colour_by = "combined_cluster", point_alpha=1)+
#	scale_color_manual(values=precast.colorList[["n1663_k9"]][["colors"]][1:8])
#p4 <- plotPCA(spe_pseudo, dimred="MNN_1663", ncomponents=4, colour_by = "combined_cluster")+
#	scale_color_manual(values=precast.colorList[["n1663_k9"]][["colors"]][1:8])
#
#ggsave("plots/06_pseudobulk/MNN-1663_eval.png", gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix=cbind(c(1,3,3),c(2,4,4))),
#        bg="white", height=8, width=12, units="in")
#cat("\nMNN eval plots saved to: plots/06_pseudobulk/MNN-1663_eval.png\n")


save(spe_pseudo, file="processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata")
cat("\nFiltered, normalized pseudobulk spe saved to: processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata\n")

#update spe tracker
write(c(paste("********** QC filtered and normalized pseudobulked spe on",format(Sys.time()),"EST"),
        "********** Old file location: processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus.Rdata",
        "********** New file location: processed-data/06_pseudobulk/spe_n119_pseudo_sample-comb-clus_norm-filt.Rdata",
        "********** Source code: code/06_pseudobulk/02_norm-QC-PCA.r",
        "**********","**********","**********"), "spe_tracker_current.txt", append=TRUE)

cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
