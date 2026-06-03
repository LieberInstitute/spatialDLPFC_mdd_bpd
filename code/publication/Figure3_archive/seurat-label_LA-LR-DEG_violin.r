setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")

#results_set="smoothed-k9-1663"
#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")

results_set="seurat-pc30"
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")


la.degs = read.csv(paste0("processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))
lr.degs = read.csv(paste0("processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_", results_set,
                          "_dx-sex_degs-F-test-t-test.csv"))

#plot.genes = c("BCL6","CEBPD","APOLD1","ELK1","SST")
plot.genes = c("CX3CR1","C3","SELPLG","RGS1","C1QB")

for (j in plot.genes) {
  colData(spe_pseudo)[[gsub("-","\\.", j)]] = logcounts(spe_pseudo)[rowData(spe_pseudo)$gene_name==j,]
}

summ.la.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes),
         seurat_label="L-A")
summ.lr.df = as.data.frame(colData(spe_pseudo)[,c("condition", "sex", "seurat_label", plot.genes)]) %>%
  tidyr::pivot_longer(all_of(plot.genes), names_to="key_genes", values_to="logcounts") %>%
  mutate(key_genes= factor(key_genes, levels=plot.genes))
all.df = bind_rows(summ.la.df, summ.lr.df) %>% mutate(cluster=factor(seurat_label, levels=c("L-A", names(cpList$transfer.bright)[c(1:4,8,5:7)])))

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
  geom_violin(aes(fill=condition), scale="width", color="transparent", position = position_dodge(width=.8), trim=F, bounds=c(0,11))+
  scale_fill_manual(values=cpList$dx.pal, guide="none")+
  facet_grid(cols=vars(sex), rows=vars(key_genes), switch="y",
             labeller= as_labeller(c("F"="Female","M"="Male", plot.genes)))+
  geom_text(data=filter(all.df_filt2, cluster=="L4"), aes(x=cluster, y=10, label=key_genes), 
          color="grey50", size=2, fontface="italic", hjust=.5)+
  geom_crossbar(data=all.df_filt2, aes(color=condition, y=ypos, ymax=ypos+yse, ymin=ypos-yse), 
	position = position_dodge(width=.8), linewidth=.3)+
  scale_color_manual(values=c("black","black","black"), guide="none")+
  scale_y_continuous(position="right", breaks=c(0,2,4,6,8,10))+coord_cartesian(ylim=c(0,11))+
  scale_x_discrete(labels=c("LA","MV","Ast","L2.3","L4","Inb","L5","L6","Olg"))+
  theme_bw()+theme(strip.background = element_rect(fill="transparent", color="transparent"),
                   text=element_text(size=8), axis.text=element_text(size=6),
                   axis.title.x=element_blank(), axis.title.y=element_blank(),
                   axis.ticks = element_line(linewidth=.2), strip.text.y.left = element_blank(),
                   panel.grid.minor=element_blank(), panel.grid.major=element_line(linewidth=.2))


#ggsave(file="plots/publication/Figure3/LA-LR-DEG_seurat-pc30_deg-violin_5-genes.pdf", p2, height=5.5, width=3)
ggsave(file="plots/publication/Figure3/LA-LR-DEG_seurat-pc30_deg-violin_immune.pdf", p2, height=5.5, width=3)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
