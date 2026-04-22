
setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  #library(HDF5Array)
  #library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(scater)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
#load in dlpfc marker genes (manually curated)
source("code/06_pseudobulk/dlpfc_genes.r")
dlpfc.genes = dlpfc.genes[c(1:5,9,6:8)]

##### smoothed version first #####

#load in sce for heatmap and dotplots
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
precast_levels= c("L1","L2","L3.4","L5","L6","WM")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$smoothed_k9_1663),
                            levels=as.character(outer(cond_sex, precast_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

#dlpfc.genes.ids = rownames(spe_summ)[rowData(spe_summ)$gene_name %in% unlist(dlpfc.genes)]
#cat("\nNumber of dlPFC genes to highlight:", length(dlpfc.genes.ids),"\n\n")

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")
sm.df = dotplotDF(spe_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
	clusters=factor(clusters, levels=levels(spe_summ$sample_id), labels=gsub(" ","\n", levels(spe_summ$sample_id))))
sm.df$condition = factor(substr(sm.df$clusters, start=0, stop=3), levels=c("NTC","MDD","BPD")) 
sm.df$sex = factor(substr(sm.df$clusters, start=5, stop=5), levels=c("F","M"))
#dlpfc marker dotplot
c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
sm.df2 = mutate(sm.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))

#make labels
xlbs = levels(sm.df$clusters)
## keep only the odd numbered dx indicators
odd1 = seq(1, length(xlbs), by=2)
xlbs[-odd1] = substr(xlbs[-odd1], start=4, stop=11)

p1 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f))+
  geom_tile(data=sm.df2, aes(fill=fill_color), alpha=.5)+
  scale_fill_manual(values=cpList$low.res.light, guide="none")+
  geom_count(aes(shape=sex, color=mean_expr_scaled, size=prop_spots))+
  scale_shape_manual(values=c(20,18))+
  scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=xlbs)+
  scale_size(range=c(1,5), limits=c(0,1), breaks=c(0,.5,1))+
  guides(size = guide_legend(override.aes = list(shape = 20)),
	 shape = guide_legend(overrisde.aes = list(size=5)))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="common marker genes")+
  theme_minimal()+theme(axis.title.x=element_blank(), axis.text.x=element_text(size=7, hjust=0),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"))
#                        axis.title.y=element_text(margin=margin(0,20,0,40,"pt")))

ggsave(file="plots/publication/Figure1/dlpfc-genes_smoothed_dx-sex-dotplot.pdf", p1, height=7, width=7)

#version not by dx*sex
sm.df = dotplotDF(spe_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=T, cluster_labels="smoothed_k9_1663", row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
         clusters=factor(clusters, levels=precast_levels))

c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
sm.df2 = mutate(sm.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))

p2 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f))+
  geom_tile(data=sm.df2, aes(fill=fill_color), alpha=.5)+
  scale_fill_manual(values=cpList$low.res.light, guide="none")+
  geom_count(aes(color=mean_expr_scaled, size=prop_spots))+
  #scale_shape_manual(values=c(20,18))+
  scale_color_gradient(low="white", high="black")+
  #scale_x_discrete(labels=xlbs)+
  scale_size(range=c(1,5), limits=c(0,1), breaks=c(0,.5,1))+
  guides(size = guide_legend(override.aes = list(shape = 20)),
         shape = guide_legend(overrisde.aes = list(size=5)))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="common marker genes")+
  theme_minimal()+theme(axis.title.x=element_blank(), axis.text.x=element_text(size=7),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"))
#                        axis.title.y=element_text(margin=margin(0,20,0,40,"pt")))

ggsave(file="plots/publication/Figure1/dlpfc-genes_smoothed_dotplot.pdf", p2, height=6, width=3.2)

##### now unsmoothed version #####

#load in sce for heatmap and dotplots
load("processed-data/06_pseudobulk/PRECAST/spe_n119_pseudo-dotplot_precast-n1663-k9.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
precast_levels= c("Vasc","L1","L2","L3.4","GABA","L5","L6","WM","low UMI")

sm.df = dotplotDF(spe_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
         clusters=factor(clusters, levels=precast_levels))

#dlpfc marker dotplot
c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
sm.df2 = mutate(sm.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))

p2 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f))+
  geom_tile(data=sm.df2, aes(fill=fill_color), alpha=.5)+
  scale_fill_manual(values=cpList$low.res.light, guide="none")+
  geom_count(aes(color=mean_expr_scaled, size=prop_spots))+
  #scale_shape_manual(values=c(20,18))+
  scale_color_gradient(low="white", high="black")+
  #scale_x_discrete(labels=xlbs)+
  scale_size(range=c(1,5), limits=c(0,1), breaks=c(0,.5,1))+
  guides(size = guide_legend(override.aes = list(shape = 20)),
         shape = guide_legend(overrisde.aes = list(size=5)))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots",
       y="common marker genes")+
  theme_minimal()+theme(axis.title.x=element_blank(), axis.text.x=element_text(size=7),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(15,"pt"))
#                        axis.title.y=element_text(margin=margin(0,20,0,40,"pt")))

ggsave(file="plots/publication/Figure1/supp_dlpfc-genes_original_dotplot.pdf", p2, height=6, width=3.5)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()


