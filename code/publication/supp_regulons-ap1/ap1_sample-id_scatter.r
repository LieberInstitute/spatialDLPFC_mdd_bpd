setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

set.seed(123)
setAutoBlockSize(1e9)

# only regulons with >10 components
regulons <- read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons1 = filter(regulons, set_size>=10)$TF

#load aucell
aucell= read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))
aucell$sample_id = substr(rownames(aucell), start=20, stop=50)  

# plot avg aucell for each regulon for each sample
df1 = group_by(aucell, condition, sex, sample_id) %>%
  summarise_at(regulons1, mean) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")),
         sex=factor(sex, levels=c("F","M")))
  
df1 = tidyr::pivot_longer(df1, all_of(regulons1), names_to="regulon", values_to="avg_aucell")


ap1.order = filter(ungroup(df1), regulon %in% c("FOS","JUNB","FOSB")) %>% 
  group_by(sample_id) %>% summarise(avg_aucell=mean(avg_aucell)) %>%
  ungroup() %>% mutate(rank_ap1= rank(-avg_aucell)) %>%
  arrange(rank_ap1) %>% pull(sample_id)

p3 <- ggplot(filter(df1, regulon %in% c("FOS","JUNB","FOSB")) %>%
         mutate(y_labs= factor(sample_id, levels=ap1.order)), 
       aes(x=y_labs, y=avg_aucell, color=regulon))+
  geom_point(size=.3)+geom_vline(aes(xintercept=10.5))+
  labs(x="sample_id (ordered by AP1 regulons)", y="avg. AUCell", color="Initial regulon",
       title="AP1 donor-enriched regulons")+
  theme_bw()+theme(axis.text.x=element_blank(), legend.position="bottom", 
                   panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
		text=element_text(size=6), axis.ticks = element_line(linewidth=.3))

load("processed-data/06_pseudobulk/spe_n119_pseudo-dotplot-no-lowUMI_sample-id.Rdata")
cpList <- readRDS("plots/colorPalettes.rds")

npas4.df = data.frame(sample_id=spe_summ$sample_id, 
                      condition=spe_summ$condition, sex=spe_summ$sex,
                  npas4_prop.detected=assay(spe_summ, "logcounts.prop.detected")[rowData(spe_summ)$gene_name=="NPAS4",]) %>%
  mutate(y_labs=factor(sample_id, levels=ap1.order)) 


p4 <- ggplot(npas4.df, aes(x=y_labs, y=npas4_prop.detected, color=condition, shape=sex))+
  geom_point(size=.3)+scale_color_manual(values=cpList$dx.pal)+
  scale_shape_manual(values=c("F"=19, "M"=17))+
  ylim(0,1)+geom_vline(aes(xintercept=10.5))+
  labs(x="sample_id (ordered by AP1 regulons)", y="prop. spots with NPAS4 detected", 
       title="NPAS4 outliers in donor-enriched AP1 regulons")+
  theme_bw()+theme(axis.text.x=element_blank(), legend.position="bottom", 
                   panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
		text=element_text(size=6), axis.ticks = element_line(linewidth=.3))

ggsave(file="plots/publication/supp_regulons-ap1/ap1-sample-id_scatter.pdf", 
	grid.arrange(p3, p4, ncol=1), width=3, height=4)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
