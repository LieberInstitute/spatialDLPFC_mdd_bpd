setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

source("code/publication/plotting_utils.r")

mod_gene="GLUL"
plot.genes = c("APOLD1","ANGPTL4","PFKFB3","TIPARP","VEGFA","DDIT4",
	"VASN","ITPKC","PDLIM4","EDN1")

o1 = filter(refined.modules, TF==mod_gene, target %in% plot.genes) %>% arrange(desc(importance))
o2 = filter(refined.modules, TF==mod_gene) %>% slice_max(n=1, importance) %>%
  mutate(target=mod_gene)
imp.df <- bind_rows(o2, o1) %>%
  mutate(y_lab=factor(target, levels=c(rev(o1$target), mod_gene)),
         bar_format= factor(ifelse(target==mod_gene, "self", "normal"), levels=c("self","normal")))

plot.genes = imp.df$target
length(plot.genes) #28

p3 <- getDotplot(plot.genes, de.df)
p3.1 <- getMeanRatioBar(plot.genes, sce_summ)
p3.2 <- getDetectedBoxplot(plot.genes, spe_summ)
p3.3 <- ggplot(mutate(imp.df, y_lab=factor(y_lab, levels=rev(plot.genes))), 
               aes(y=y_lab, x=importance, fill=bar_format, lty=bar_format, color=bar_format))+
  geom_bar(stat="identity", linewidth=.3)+
  scale_fill_manual(values=c("white","grey"), guide="none")+
  scale_linetype_manual(values=c(2,1), guide="none")+
  scale_color_manual(values=c("black","grey"), guide="none")+
  labs(title=" ")+scale_y_discrete(position="right")+
  theme_minimal()+theme(axis.text.y=element_blank(), 
                        axis.title.y=element_blank(),
                        axis.text.x=element_text(size=8), 
                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
                        plot.margin = margin(.2,.0,1.8,0, unit="cm"))

ggsave(file=paste0("plots/publication/Figure4/",mod_gene,"_dotplot.pdf"), 
       arrangeGrob(grobs=list(p3, p3.1, p3.2, p3.3), layout_matrix=matrix(c(1,1,1,1,1,2,3,4), ncol=8)),
       height=3, width=6.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
