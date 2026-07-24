setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(ggrastr)
})

load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_control_MBv-filtered_processed-SCT.Rdata")

mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv") %>%
  group_by(gene_name) %>% slice_max(MeanRatio, n=1) %>%
  mutate(cellType.target= ifelse(MeanRatio<1.5, "multi", cellType.target)) %>% ungroup()

plot.genes = c("ITIH5","COL5A3","PLP1","C3","GAD1","RALYL","SLC6A1")

plist <- FeaturePlot(seu_con, features=plot.genes, pt.size=.1, raster=F, cols=c("lightgrey","black"), combine=F)
names(plist) <- plot.genes

plist_format <- lapply(plot.genes, function(x) {
  prop.main = filter(mratio.sn, gene_name==x)
  plot.title= paste0(x, " (", signif(prop.main$prop.detected*100, 3), "% of ", prop.main$cellType.target, ")")
  p1 <- plist[[x]]+ggtitle(plot.title)+scale_color_gradient(low="lightgrey", high="black", labels=function(x) sprintf("%.1f", x))+
    theme(text=element_text(size=6), axis.title=element_blank(), axis.text=element_blank())
  rasterize(p1, layers="Point", dpi=300)
})

ggsave(file="plots/publication/supp_mean-ratio/umap_example-DEGs.pdf", 
       marrangeGrob(plist_format, ncol=1, nrow=1, top=NULL),
       height=6, width=6)


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
