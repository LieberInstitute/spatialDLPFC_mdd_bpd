setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SingleCellExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

#load in sce objects and fix colnames
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_seurat-low-res.Rdata")
colnames(sce_summ) <- sce_summ$seurat_low.res

load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
precast_levels= c("L1","L2","L3.4","L5","L6","WM")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$smoothed_k9_1663),
                            levels=as.character(outer(cond_sex, precast_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id
spe_summ_sm <- spe_summ

load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb")
spe_summ$seurat_label = factor(spe_summ$seurat_label, levels=seurat_levels)
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                                 levels=as.character(outer(cond_sex, seurat_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id
spe_summ_pc30 <- spe_summ

#load in functions, gene list, and color palettes
source("code/06_pseudobulk/custom_functions.r")
source("code/06_pseudobulk/dlpfc_genes.r")
cpList <- readRDS("plots/colorPalettes.rds")

#load in enrichment results
enrich.df_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_all-results.csv")
enrich.df_sn$seurat_label = factor(enrich.df_sn$seurat_label, levels=names(cpList$low.res.bright))
t_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_t-stat.csv", row.names=1)
lf_sn <- read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/layer-enrichment_control-low-res_logFC.csv", row.names=1)


enrich.df_smooth <- read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_all-results.csv")
enrich.df_smooth$smoothed = factor(enrich.df_smooth$smoothed, levels=names(cpList$smoothed.bright))
t_sm <-	read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_t-stat.csv", row.names=1)
lf_sm <- read.csv("processed-data/06_pseudobulk/PRECAST_smoothed/layer-enrichment_smoothed-k9-1663_logFC.csv", row.names=1)


enrich.df_pc30 <- read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_all-results.csv")
enrich.df_pc30$seurat_label = factor(enrich.df_pc30$seurat_label, levels=names(cpList$transfer.bright))
t_pc30 <- read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_t-stat.csv", row.names=1)
lf_pc30 <- read.csv("processed-data/06_pseudobulk/Seurat/layer-enrichment_seurat-pc30_logFC.csv", row.names=1)



#aligned volcano plots
snList = VolcanoPlot(enrich.df_sn, "seurat_label")
sn.v = arrangeGrob(grobs=snList, ncol=1, top="SZBDMulti-seq")

smList = VolcanoPlot(enrich.df_smooth, "smoothed")
#gridExtra won't populate the last entry with an NA so make an empty plot 
smList[["blank"]]  = ggplot(enrich.df_smooth, aes(x=logFC, y=-log10(adj.P.Val)))+geom_blank()+theme_void()
sm.v= arrangeGrob(grobs=smList,
	layout_matrix=as.matrix(c(NA,1,2,3,NA,4,5,6,7)),
        top="PRECAST (smoothed)")

pc30List = VolcanoPlot(enrich.df_pc30, "seurat_label")
pc30.v = arrangeGrob(grobs=pc30List, 
	layout_matrix=as.matrix(c(1,2,3,NA,4,5,6,7,8)),
	top="MBv label transfer")

ggsave(file="plots/06_pseudobulk/compare-enrichment_volcano.png",
	grid.arrange(sn.v, sm.v, pc30.v, ncol=3,
		top="Black: top 50 t-stat and top 50 logFC"),
	bg="white", width=8.5, height=11)



#t-stat correlation plots

#load in spatialDLPFC results
layer_modeling_results <- spatialLIBD::fetch_data(type = "modeling_results")
### i guess these are from the Maynard paper.... (aka manual annotations)
#### https://www.bioconductor.org/packages/release/data/experiment/vignettes/spatialLIBD/inst/doc/spatialLIBD.html#spatiallibd-functions
#### We already covered fetch_data() which allows you to download the Human DLPFC Visium data from LIBD researchers and colleagues (Maynard, Collado-Torres, Weber et al., 2021).
t1 = layer_modeling_results$enrichment[,c(grep("t_stat",colnames(layer_modeling_results$enrichment), value=T),"ensembl")]
colnames(t1) = gsub("t_stat_","", colnames(t1))
colnames(t1) = gsub("Layer","L", colnames(t1))
m1 = as.matrix(t1[,1:7])
rownames(m1) = t1$ensembl
cat("\n\nDimensions of spatialDLPFC t stat matrix:\n")
dim(m1) #22331 7


#top 50 t stat
top50.list <- list("SZBDMulti-seq"=getTopGenes(t_sn, top_n=50),
                    "PRECAST (smoothed)"=getTopGenes(t_sm, top_n=50),
                    "MBv label transfer"=getTopGenes(t_pc30, top_n=50)
)
cat("\n\nNumber of marker genes identified by top 50 t statistic (per dataset):\n\n")
sapply(top50.list, length)


top50.lf.list <- list("SZBDMulti-seq"=getTopGenes(lf_sn, top_n=50),
                      "PRECAST (smoothed)"=getTopGenes(lf_sm, top_n=50),
                      "MBv label transfer"=getTopGenes(lf_pc30, top_n=50)
) 
cat("\n\nNumber	of marker genes	identified by top 50 logFC (per dataset):\n\n")
sapply(top50.lf.list, length)


top50.both.list = lapply(names(top50.lf.list), function(x) union(top50.lf.list[[x]], top50.list[[x]]))
names(top50.both.list) <- names(top50.lf.list)
cat("\n\nNumber of marker genes identified by union of top 50 t stat and top 50 logFC (per dataset):\n\n")
sapply(top50.both.list, length)

#list for factor re-ordering
orderList <- list("spatialDLPFC"=c("L1","L2","L3","L4","L5","L6","WM"),
                  "SZBDMulti-seq"=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
                  "PRECAST (smoothed)"=c("L1","L2","L3.4","L5","L6","WM"),
                  "MBv label transfer"=c("Micro.Vasc","Astro","L2.3","L4","L5","L6","Oligo","Inhb"))

#generate cor plot heatmaps via custom function
p1 <- plotCorHeatmap(query_list_level= "SZBDMulti-seq", query_stats= t_sn,
               reference_list_level= "spatialDLPFC", reference_stats= m1)

p2 <- plotCorHeatmap(query_list_level= "PRECAST (smoothed)", query_stats= t_sm,
                     reference_list_level= "spatialDLPFC", reference_stats= m1)

p3 <- plotCorHeatmap(query_list_level= "MBv label transfer", query_stats= t_pc30,
                     reference_list_level= "spatialDLPFC", reference_stats= m1)

ggsave(file="plots/06_pseudobulk/t-stat-cor-heatmap_all-vs-spatialDLPFC.png",
       grid.arrange(p1, p2, p3, ncol=3),
       bg="white", width=9, height=3)
cat("\nSaved t stat correlation plots to: plots/06_pseudobulk/t-stat-cor-heatmap_all-vs-spatialDLPFC.png\n")

p4 <- plotCorHeatmap(query_list_level= "SZBDMulti-seq", query_stats= t_sn,
                     reference_list_level= "MBv label transfer", reference_stats= t_pc30,
                     .coord_flip=T)

p5 <- plotCorHeatmap(query_list_level= "MBv label transfer", query_stats= t_pc30,
                     reference_list_level= "SZBDMulti-seq", reference_stats= t_sn)

ggsave(file="plots/06_pseudobulk/t-stat-cor-heatmap_MBv-label-transfer-vs-SZBDMulti-seq.png",
       grid.arrange(p4, p5, ncol=2),
       bg="white", width=6, height=3) 
cat("\nSaved t stat correlation plots to: plots/06_pseudobulk/t-stat-cor-heatmap_MBv-label-transfer-vs-SZBDMulti-seq.png\n")


p6 <- plotCorHeatmap(query_list_level= "PRECAST (smoothed)", query_stats= t_sm,
                     reference_list_level= "MBv label transfer", reference_stats= t_pc30,
                     .coord_flip=T)

p7 <- plotCorHeatmap(query_list_level= "MBv label transfer", query_stats= t_pc30,
                     reference_list_level= "PRECAST (smoothed)", reference_stats= t_sm)

ggsave(file="plots/06_pseudobulk/t-stat-cor-heatmap_MBv-label-transfer-vs-PRECAST-smoothed.png",
       grid.arrange(p6, p7, ncol=2),
       bg="white", width=6, height=3) 
cat("\nSaved t stat correlation plots to: plots/06_pseudobulk/t-stat-cor-heatmap_MBv-label-transfer-vs-PRECAST-smoothed.png\n")



#top 10 marker dotplots

#need to load rdata
# load in unfiltered spe object which should have whole gene universe
### the information here (https://www.synapse.org/Synapse:syn22963646) indicates that they used
### the hg38 genome for their reference which is what we used too (?)
spe <- HDF5Array::loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
rdata = as.data.frame(rowData(spe))
dim(rdata) #36601 7
rm(spe)

#make dotplot dataframes
top.genes.list = list(filter(enrich.df_sn, logFC>0) %>% group_by(seurat_label) %>% slice_min(n=10, adj.P.Val) %>% mutate(fill_color=seurat_label),
                      filter(enrich.df_smooth, logFC>0) %>% group_by(smoothed) %>% slice_min(n=10, adj.P.Val) %>% mutate(fill_color=smoothed),
                      filter(enrich.df_pc30, logFC>0) %>% group_by(seurat_label) %>% slice_min(n=10, adj.P.Val) %>% mutate(fill_color=seurat_label))
names(top.genes.list) <- c("SZBDMulti-seq","PRECAST (smoothed)","MBv label transfer")

rm.genes = c(setdiff(top.genes.list[[1]]$gene_id, rowData(spe_summ_sm)$gene_id), 
	setdiff(top.genes.list[[1]]$gene_id, rowData(spe_summ_pc30)$gene_id),
             setdiff(top.genes.list[[2]]$gene_id, rowData(sce_summ)$gene_id), 
	setdiff(top.genes.list[[2]]$gene_id, rowData(spe_summ_pc30)$gene_id),
             setdiff(top.genes.list[[3]]$gene_id, rowData(sce_summ)$gene_id), 
	setdiff(top.genes.list[[3]]$gene_id, rowData(spe_summ_sm)$gene_id)
)
gene_f = setdiff(do.call(c, lapply(top.genes.list, function(x) x$gene_id)), rm.genes)

top.df.list <- lapply(top.genes.list, function(i) {
  gene_f = filter(i, !gene_id %in% rm.genes)$gene_name
  dplot.sn <- dotplotDF(sce_summ, setdiff(i$gene_id, rm.genes)) %>%
    mutate(clusters= factor(clusters, levels=levels(sce_summ$seurat_low.res)),
           gene_name_f=factor(gene_name, levels=rev(gene_f))) %>%
    left_join(i[,c("gene_name","fill_color")], by=c("gene_name"))
  dplot.sm <- dotplotDF(spe_summ_sm, setdiff(i$gene_id, rm.genes), 
                        summarize_groups=T, cluster_labels="smoothed_k9_1663") %>%
    mutate(gene_name_f= factor(gene_name, levels=rev(gene_f))) %>%
    left_join(i[,c("gene_name","fill_color")], by=c("gene_name"))
  dplot.pc30 <- dotplotDF(spe_summ_pc30, setdiff(i$gene_id, rm.genes), 
                          summarize_groups=T, cluster_labels="seurat_label") %>%
    mutate(gene_name_f= factor(gene_name, levels=rev(gene_f))) %>%
    left_join(i[,c("gene_name","fill_color")], by=c("gene_name"))
  return(list("SZBDMulti-seq"=dplot.sn,
              "PRECAST (smoothed)"=dplot.sm,
              "MBv label transfer"=dplot.pc30))
})

#for all 3 annotation strategies, plot top 10 most sig SZBDMulti-seq genes
plots1 = dotplotList("SZBDMulti-seq", grob_prop=.5)
pdf(file="plots/06_pseudobulk/top-10-sig-each-annot_SZBDMulti-seq-genes_dotplot.pdf", 
    height=11, width=6)
marrangeGrob(plots1, ncol=1, nrow=1, top=NULL)
dev.off()

#for all 3 annotation strategies, plot top 10 most sig PRECAST (smoothed) genes
plots2 = dotplotList("PRECAST (smoothed)", grob_prop=.7)
pdf(file="plots/06_pseudobulk/top-10-sig-each-annot_smoothed-k9-1663-genes_dotplot.pdf", height=11, width=6)
marrangeGrob(plots2, ncol=1, nrow=1, top=NULL)
dev.off()

#for all 3 annotation strategies, plot top 10 most sig MBv label transfer genes
plots3 = dotplotList("MBv label transfer", grob_prop=.9)
pdf(file="plots/06_pseudobulk/top-10-sig-each-annot_MBv-seurat-pc30-genes_dotplot.pdf", height=11, width=6)
marrangeGrob(plots3, ncol=1, nrow=1, top=NULL)
dev.off()

cat("\n\nSaved dotplots of all 3 annotation strategies for...",
	"\n>> Top 10 most sig. SZBDMulti-seq genes: plots/06_pseudobulk/top-10-sig-each-annot_SZBDMulti-seq-genes_dotplot.pdf",
	"\n>> Top 10 most sig. PRECAST (smoothed) genes: plots/06_pseudobulk/top-10-sig-each-annot_smoothed-k9-1663-genes_dotplot.pdf",
	"\n>> Top 10 most sig. MBv label transfer genes: plots/06_pseudobulk/top-10-sig-each-annot_MBv-seurat-pc30-genes_dotplot.pdf\n\n") 



#spot level annotation assignment for MBv label transfer and PRECAST (smoothed)
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata$smoothed_k9_1663_f = factor(cdata$smoothed_k9_1663, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"), 
                                labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))
                                #labels=c("L1","L2","L3.4","L5","L6","WM","low UMI","dropped","dropped"))

res.pc30 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(res.pc30)))

cdata$seurat_label = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"),
                           labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","L5","L6","Oligo","Inhb"))

cdata$prediction.score.max = res.pc30$prediction.score.max

##filter out low UMI spots
#cdata = filter(cdata, smoothed_k9_1663_f!="low UMI")

spots.df = group_by(cdata, smoothed_k9_1663_f, seurat_label) %>% 
  summarise(n=n(), avg.score = mean(prediction.score.max)) %>%
  mutate(smoothed_k9_1663_f= factor(smoothed_k9_1663_f, levels=rev(levels(cdata$smoothed_k9_1663_f))))

p <- ggplot(spots.df,  aes(x=seurat_label, y=smoothed_k9_1663_f, size=n, color=avg.score))+
  geom_count()+
  scale_color_gradientn("Avg. Seurat\nprediction\nscore",
                        colors=RColorBrewer::brewer.pal(n=5, "Purples"),
                        limits=c(0,1))+
  scale_size("# spots", range=c(0,6), breaks=c(10,30,50)*1000,
             labels=function(x) paste0(x/1000,"k"))+
  scale_x_discrete("MBv label transfer", labels=c("M/V","Astro","L2.3","L4","L5","L6","Oligo","Inhb"))+
  scale_y_discrete("PRECAST (smoothed)")+
  labs(title="Spot-level annotation")+#, 
	##added this
	#subtitle="no low UMI")+
  theme_minimal()+theme(panel.grid= element_blank(), aspect.ratio=1,
                        legend.key.size= unit(10, "pt"), legend.title = element_text(size=8),
                        legend.text = element_text(size=7),
                        plot.title=element_text(size=10), plot.subtitle = element_text(size=8),
                        axis.text.x= element_text(angle=45, hjust=1))

ggsave(file="plots/06_pseudobulk/spot-level-annotation-dotplot_MBv-label-transfer-vs-PRECAST-smoothed.png",
       p,
       bg="white", width=3.4, #height=3.4)
	height=3) 

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
