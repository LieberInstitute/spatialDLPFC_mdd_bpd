setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

results_set="smoothed-k9-1663"
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

#results_set="seurat-pc30"
#load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")


la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

plot.genes = c("BCL6","CEBPD","APOLD1","ELK1","SST")

## previous gene lists by themes
#plot.genes= c("BCL6","HSPA1B","JUN") #WM-inflamm,
#plot.genes= c("CLDN11","MAG","SLC44A1") #WM-myelin
#plot.genes = c("UBA52","EIF5B","RPS8") #
#plot.genes= c("CEBPD","GADD45B","ANGPTL4") #BBB angio inflamm group
#plot.genes= c("SST", "CORT", "CRH") #InhN that are L-R too
#plot.genes= c("ELK1","DUSP6","RASD1") #MAPK

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         smoothed_k9_1663="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "smoothed_k9_1663", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(smoothed_k9_1663, levels=c("L-A", names(cpList$smoothed.bright))))

names(plot.genes) <- plot.genes

all.df_filt = filter(all.df, key_genes %in% plot.genes)

#calc sd/se
all.df_filt2 = group_by(all.df_filt, condition, sex, cluster, key_genes) %>% summarise(ypos=mean(logcounts), ysd=sd(logcounts), n=n(), yse=ysd/sqrt(n)) #%>%

#remove outliers
all.df_filt3 = group_by(all.df_filt, condition, sex, cluster, key_genes) %>% 
	mutate(iqr=IQR(logcounts), q1=quantile(logcounts, probs=c(.25)), q3= quantile(logcounts, probs=c(.75)),
		out_low=logcounts< (q1-1.5*iqr), out_high= logcounts> (q3+1.5*iqr)) %>%
	filter(out_low==F, out_high==F)

#ceiling(all.df_filt3$logcounts) #10
p2 <- ggplot(all.df_filt3, aes(x=cluster, y=logcounts))+
  geom_violin(aes(fill=condition), scale="width", color="transparent", position = position_dodge(width=.8), trim=F, bounds=c(0,10))+
  scale_fill_manual(values=cpList$dx.pal, guide="none")+
  facet_grid(cols=vars(sex), rows=vars(key_genes), switch="y",
             labeller= as_labeller(c("F"="Female","M"="Male", plot.genes)))+
  geom_text(data=filter(all.df_filt2, cluster=="L3.4"), aes(x=cluster, y=1, label=key_genes), 
          color="grey50", size=2, fontface="italic", hjust=.5)+
  geom_crossbar(data=all.df_filt2, aes(color=condition, y=ypos, ymax=ypos+yse, ymin=ypos-yse), 
	position = position_dodge(width=.8), linewidth=.3)+
  scale_color_manual(values=c("black","black","black"), guide="none")+
  scale_y_continuous(position="right", breaks=c(0,2,4,6,8,10))+coord_cartesian(ylim=c(0,10))+
  theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                   text=element_text(size=8), axis.text=element_text(size=6),
                   axis.title.x=element_blank(), axis.title.y=element_blank(),
                   axis.ticks = element_line(linewidth=.2), strip.text.y.left = element_blank(),
                   panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2))


ggsave(file="plots/publication/Figure2/deg-violin_5-genes.pdf", p2, height=5.5, width=3)

#ggsave(file="plots/publication/Figure2/deg-boxplot_WM-inflamm.pdf", p3, height=3, width=3)
#ggsave(file="plots/publication/Figure2/deg-boxplot_WM-myelin.pdf", p3, height=3, width=3)
#ggsave(file="plots/publication/Figure2/deg-boxplot_WM-ribo.pdf", p3, height=3, width=3)
#ggsave(file="plots/publication/Figure2/deg-boxplot_transl-starv.pdf", p3, height=3, width=3)
#ggsave(file="plots/publication/Figure2/deg-boxplot_BBB-inflamm.pdf", p3, height=3, width=3)
#ggsave(file="plots/publication/Figure2/deg-boxplot_GABA-pep.pdf", p3, height=3, width=3)
#ggsave(file="plots/publication/Figure2/deg-boxplot_MAPK.pdf", p3, height=3, width=3)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
