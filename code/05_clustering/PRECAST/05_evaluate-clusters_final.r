setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(Seurat)
	library(scater)
	library(bluster)
	library(ggspavis)
	library(escheR)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})

source("code/05_clustering/PRECAST/PRECAST_colorLists.r")
source("code/05_clustering/PRECAST/05-supp_plot-functions.r")

#load in spe
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")

spe$precast_k9_1663_f=factor(spe$precast_k9_1663, levels=precast.colorList[["n1663_k9"]][["clusters"]], labels=precast.colorList[["n1663_k9"]][["annotation"]])
spe$precast_k9_1079_f=factor(spe$precast_k9_1079, levels=precast.colorList[["n1079_k9"]][["clusters"]], labels=precast.colorList[["n1079_k9"]][["annotation"]])

#load in precast object with reduced dimension embeddings
#PRECAST 1663 reduced dims
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1663.Rdata")
reddim = seuInt@reductions$PRECAST@cell.embeddings
reducedDim(spe, "PRECAST_1663", withDimnames=F) = reddim[spe$seurat_key,]

# PRECAST 1079 reduced dims
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1079.Rdata")
reddim = seuInt@reductions$PRECAST@cell.embeddings
reducedDim(spe, "PRECAST_1079", withDimnames=F) = reddim[spe$seurat_key,]

#reduced dimension embeddings plots
out1 = reducedDimPlots(spe, "PRECAST_1663")
ggsave("plots/05_clustering/PRECAST/PRECAST_n1663-k9_reduced-dim-plots.png",
       grid.arrange(out1[[1]], out1[[2]], out1[[3]], layout_matrix=matrix(c(1,2,3,3))), bg="white", width=8, height=12, units="in")
cat("\nPRECAST n=1663 embeddings plots saved to: plots/05_clustering/PRECAST/PRECAST_n1663-k9_reduced-dim-plots.png\n")

out2 = reducedDimPlots(spe, "PRECAST_1079")
ggsave("plots/05_clustering/PRECAST/PRECAST_n1079-k9_reduced-dim-plots.png",
       grid.arrange(out2[[1]], out2[[2]], out2[[3]], layout_matrix=matrix(c(1,2,3,3))), bg="white", width=8, height=12, units="in")
cat("\nPRECAST n=1079 embeddings plots saved to: plots/05_clustering/PRECAST/PRECAST_n1079-k9_reduced-dim-plots.png\n")

#cluster example and proportion plots
(n16.samples = paste0("V13B23-30",c(1,2,8,9)))
spe_sub2 = spe[,spe$slide %in% n16.samples]

# PRECAST 1663 example spot plots
p4 <- plotSpots(spe_sub2, sample_id="sample_id", annotate="precast_k9_1663_f", point_size=.1)+
  scale_color_manual("PRECAST\nn=1663\nk=9", values=precast.colorList[["n1663_k9"]][["colors"]])+
  facet_wrap(vars(sample_id), ncol=4)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))

gabaList = locationSpotPlots(spe_sub2, "PRECAST_1663")
ggsave("plots/05_clustering/PRECAST/PRECAST_n1663-k9_example-GABA-location.png",
	do.call(grid.arrange, c(gabaList, ncol=4)), bg="white", height=9, width=9, units="in")
cat("\nPRECAST n=1663 example plots and proportion plots saved to: plots/05_clustering/PRECAST/PRECAST_n1663-k9_example-GABA-location.png\n")

# PRECAST 1663 proportion of cluster per sample (tissue composition)
prop.df = mutate(as.data.frame(colData(spe)), sample_id=as.factor(sample_id)) %>%
  group_by(sample_id, precast_k9_1663_f, .drop=F) %>% tally() %>% 
  group_by(sample_id) %>% mutate(total_n=sum(n), prop_n=n/total_n) %>%
  #have to left join condition in later otherwise it will keep 119 values for condition with ".drop=F" and that won't do
  left_join(distinct(as.data.frame(colData(spe)[,c("sample_id","condition")]))) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")))


p5 <- ggplot(prop.df, aes(x=sample_id, y=prop_n, fill=precast_k9_1663_f))+
  geom_bar(stat="identity", width=.8)+scale_fill_manual(values=precast.colorList[["n1663_k9"]][["colors"]])+
  labs(y="proportion of spots")+
  facet_grid(precast_k9_1663_f ~ condition, scales="free")+
  theme_minimal()+theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_text(size=6),
                        strip.text.y.right=element_text(angle=0), text=element_text(size=10),
                        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank())

# PRECAST 1663 UMI per cluster
p6 <- ggplot(as.data.frame(colData(spe)), aes(x=factor(condition, levels=c("NTC","MDD","BPD")), y=sum_umi, fill=precast_k9_1663_f))+
  geom_boxplot(outliers=F)+scale_fill_manual(values=precast.colorList[["n1663_k9"]][["colors"]])+
  facet_grid(cols=vars(precast_k9_1663_f))+
  scale_y_log10()+labs(x="condition", y="sum UMI (log10 scale)")+
  theme_minimal()+theme(legend.position="none", panel.grid.major.x=element_blank(),
                        text=element_text(size=10), strip.text=element_text(size=10),
                        axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        plot.margin=unit(c(6,12,6,6),"pt"))

ggsave("plots/05_clustering/PRECAST/PRECAST_n1663-k9_cluster-example-plots.png",
       grid.arrange(p4, p5, p6, layout_matrix=matrix(c(1,1,2,2,3))), bg="white", width=6, height=12, units="in")
cat("\nPRECAST n=1663 example plots and proportion plots saved to: plots/05_clustering/PRECAST/PRECAST_n1663-k9_cluster-example-plots.png\n")

#PRECAST 1079 example spot plots
p7 <- plotSpots(spe_sub2, sample_id="sample_id", annotate="precast_k9_1079_f", point_size=.1)+
  scale_color_manual("PRECAST\nn=1079\nk=9", values=precast.colorList[["n1079_k9"]][["colors"]])+
  facet_wrap(vars(sample_id), ncol=4)+
  theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
        legend.title=element_text(margin=margin(0,0,4,0,"pt")),
        legend.box.spacing = unit(2,"pt"),
        legend.margin=margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,2,0,2,"pt"))


gabaList2 = locationSpotPlots(spe_sub2, "PRECAST_1079")
ggsave("plots/05_clustering/PRECAST/PRECAST_n1079-k9_example-GABA-location.png",
        do.call(grid.arrange, c(gabaList2, ncol=4)), bg="white", height=9, width=9, units="in")
cat("\nPRECAST n=1663 example plots and proportion plots saved to: plots/05_clustering/PRECAST/PRECAST_n1079-k9_example-GABA-location.png\n")

# PRECAST 1079 as proportion of PRECAST 1663 cluster per sample
prop.df2 = mutate(as.data.frame(colData(spe)), sample_id=as.factor(sample_id)) %>%
  group_by(sample_id, precast_k9_1663_f, precast_k9_1079_f, .drop=F) %>% tally() %>% 
  group_by(sample_id) %>% mutate(total_n=sum(n), prop_total=n/total_n) %>%
  #have to left join condition in later otherwise it will keep 119 values for condition with ".drop=F" and that won't do
  left_join(distinct(as.data.frame(colData(spe)[,c("sample_id","condition")]))) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")))


p8 <- ggplot(prop.df2, aes(x=sample_id, y=prop_total, fill=precast_k9_1079_f))+
  geom_bar(stat="identity", width=.8, position="stack")+
  scale_fill_manual("PRECAST n=1079", values=precast.colorList[["n1079_k9"]][["colors"]])+
  facet_grid(precast_k9_1663_f ~ condition, scales="free")+
  labs(y="proportion of spots")+
  theme_minimal()+theme(axis.text.x=element_blank(), axis.text.y=element_text(size=6),
                        strip.text.y.right=element_text(angle=0), text=element_text(size=10),
                        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
			legend.position="bottom", legend.key.size = unit(6, 'pt'),
                        legend.title=element_text(size=8),# margin=margin(0,0,4,0,"pt")),
                        legend.box.spacing = unit(2,"pt"),
                        legend.margin=margin(0,0,0,0,"pt"))#,
                        #legend.box.margin = margin(0,2,0,2,"pt"))

# PRECAST 1079 UMI per cluster
p9 <- ggplot(as.data.frame(colData(spe)), aes(x=factor(condition, levels=c("NTC","MDD","BPD")), y=sum_umi, fill=precast_k9_1079_f))+
  geom_boxplot(outliers=F)+scale_fill_manual(values=precast.colorList[["n1079_k9"]][["colors"]])+
  facet_grid(cols=vars(precast_k9_1079_f))+
  scale_y_log10()+labs(x="condition", y="sum UMI (log10 scale)")+
  theme_minimal()+theme(legend.position="none", panel.grid.major.x=element_blank(),
                        text=element_text(size=10), strip.text=element_text(size=10),
                        axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        plot.margin=unit(c(6,12,6,6),"pt"))


ggsave("plots/05_clustering/PRECAST/PRECAST_n1079-k9_cluster-example-plots.png",
       grid.arrange(p7, p8, p9, layout_matrix=matrix(c(1,1,2,2,3))), bg="white", width=6, height=12, units="in")
cat("\nPRECAST n=1079 example plots and	proportion plots saved to: plots/05_clustering/PRECAST/PRECAST_n1079-k9_cluster-example-plots.png\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
