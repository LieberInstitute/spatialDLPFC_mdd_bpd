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
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)),
	labels=as.character(outer(cond_sex, c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"), paste))
)
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

ggsave(file="plots/publication/Figure3/seurat-label_dlpfc-genes_dx-sex-dotplot.pdf", p1, height=7, width=7)




# now for version summarised to cluster (no dx-sex)
sm.df = dotplotDF(spe_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=T, cluster_labels="seurat_label", row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
         clusters=factor(clusters, levels=seurat_levels,
		labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg")))

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

ggsave(file="plots/publication/Figure3/seurat-label_dlpfc-genes_dotplot.pdf", p2, height=6, width=3.5)


# also plot for snRNAseq CTR
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo-dotplot_seurat-low-res.Rdata")
colnames(sce_summ) <- sce_summ$seurat_low.res
seurat_levels2 = c("Micro.Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo")

sm.df = dotplotDF(sce_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
         clusters=factor(clusters, levels=seurat_levels2,
                labels=c("M.V","Ast","L2","L3","L4","Inb","L5","L6","Olg")))

c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
sm.df2 = mutate(sm.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))

p3 <- ggplot(sm.df, aes(x=clusters, y=gene_name_f))+
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

ggsave(file="plots/publication/Figure3/label-transfer_dlpfc-genes_SZBDMultiseq-control_dotplot.pdf", p3, height=6, width=3.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()


