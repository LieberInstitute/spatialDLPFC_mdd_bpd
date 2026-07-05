setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})


# select only regulons greater than 10 components
regulons <- read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
regulons1 = filter(regulons, set_size>=10)$TF

# load aucell
aucell = read.csv("processed-data/10_SCENIC/spe-n119_13162-no-lowUMI_regulons-top20_AUCell.csv", row.names=1)
colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)] = sapply(colnames(aucell)[1:(grep("seurat_label", colnames(aucell))-1)], function(x) substr(x, start=0, stop=nchar(x)-3))
aucell$sample_id = substr(rownames(aucell), start=20, stop=50)  

# plot avg aucell for each regulon for each sample
df1 = group_by(aucell, condition, sex, sample_id) %>%
  summarise_at(regulons1, mean) %>%
  mutate(condition=factor(condition, levels=c("NTC","MDD","BPD")),
         sex=factor(sex, levels=c("F","M")))
  
df1 = tidyr::pivot_longer(df1, all_of(regulons1), names_to="regulon", values_to="avg_aucell")
df1.1 = group_by(df1, regulon) %>% summarise(avg_aucell=mean(avg_aucell))
reg.order = arrange(df1.1, avg_aucell) %>% pull(regulon)

p1 <- ggplot(mutate(df1, regulon = factor(regulon, levels=reg.order)), 
       aes(x=regulon, y=avg_aucell))+
  geom_line(aes(group=sample_id))+scale_y_continuous(labels=function(x) sprintf("%.2f", x))+
  geom_point(data= mutate(df1.1, regulon= factor(regulon, levels=reg.order)),
             color="red")+
  labs(x="", y="avg. AUCell (per donor)", title="Initial regulons (n=119)")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5), text=element_text(size=6),
	panel.grid.minor=element_blank())

ggsave(file="plots/publication/supp_regulons-ap1/n119_avg-aucell_line-graph.pdf", p1, height=2, width=3)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
