setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

#spot level annotation assignment for MBv label transfer and PRECAST (smoothed)
cdata = read.csv("processed-data/05_clustering/PRECAST/colData_all-precast-clusters.csv", row.names=1)
cdata$smoothed_k9_1663_f = factor(cdata$smoothed_k9_1663, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","Vasc","GABA"), 
                                labels=c("L1","L2","L3.4","L5","L6","WM","dropped","dropped","dropped"))
                                #labels=c("L1","L2","L3.4","L5","L6","WM","low UMI","dropped","dropped"))

res.pc30 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-MBv_qual-genes-kanchor-50-pc30_red-precast-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(res.pc30)))

cdata$seurat_label = factor(res.pc30$predicted.id, levels=c("Micro/Vasc","Astro","L2","L3","L4","Inhb","L5","L6","Oligo"),
                           labels=c("Micro.Vasc","Astro","L2.3","L2.3","L4","Inhb","L5","L6","Oligo"))

cdata$prediction.score.max = res.pc30$prediction.score.max

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
#  scale_x_discrete("MBv label transfer", labels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))+
  scale_y_discrete("PRECAST (smoothed)")+
  labs(title="Spot-level annotation", x="MBv label transfer")+
  theme_minimal()+theme(panel.grid= element_blank(), aspect.ratio=1,
                        legend.key.size= unit(10, "pt"), legend.title = element_text(size=8),
                        legend.text = element_text(size=7),
                        plot.title=element_text(size=10), plot.subtitle = element_text(size=8),
                        axis.text.x= element_text(angle=90, hjust=1, vjust=.5))

ggsave(file="plots/publication/Figure1/spot-level-annotation-dotplot_MBv-label-transfer-vs-PRECAST-smoothed.pdf",
	p, width=3.4, height=3)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
