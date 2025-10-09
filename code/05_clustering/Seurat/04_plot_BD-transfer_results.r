setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(Seurat)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

load("processed-data/05_clustering/Seurat/seurat_SZBDMulti-seq_bipolar_MBv-filtered-conservative_processed-SCT.Rdata")

#load results
res1 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-azimuth.csv", row.names=1)
stopifnot(identical(rownames(res1), colnames(seu_bd)))
res2 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered-conservative_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(res2), colnames(seu_bd)))


new.levels = c("Micro","Immune","Endo","PC","SMC","VLMC",
        "Astro",
        "L2/3 IT", "L4 IT", "L5 IT", "L5 ET", "L5/6 NP",
        "L6 IT", "L6 IT Car3","L6b","L6 CT",
        "Oligo","OPC",
        "Lamp5", "Lamp5 Lhx6", "Vip", "Sncg", "Pax6",
        "Sst","Sst Chodl","Pvalb","Chandelier")

seu_bd$azimuth = factor(as.character(seu_bd$azimuth), levels=new.levels)
seu_bd$transfer_pc30_azimuth = factor(res1$predicted.id, levels=new.levels)
seu_bd$transfer_pc30_low.res = factor(res2$predicted.id, levels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Oligo","Inhb"))

#plot annotations
low.res.pal = c("Astro"="#cfa45c","Micro/Vasc"="#911223",
	"Inhb"="#9377AC",
        "L2"="#5D9940", "L3"="#5095CD", 
        "L4"="#85A0A0",
        "L5"="#ddc94e","L6"="#E45C5F",
        "Oligo"="#D1C4B0")


cat("\nUMAP plots of annotations...\n")
p1 <- DimPlot(seu_bd, group.by=c("azimuth","transfer_pc30_azimuth","transfer_pc30_low.res"),
        pt.size=.1, raster=F, label=T, combine=F)

cdata = group_by(seu_bd@meta.data, azimuth, transfer_pc30_azimuth, .drop=F) %>% tally() %>%
  group_by(azimuth) %>% mutate(total_n=sum(n)) %>%
  ungroup() %>%
  mutate(prop_n=n/total_n,
         azimuth_rev = factor(azimuth, levels=rev(levels(seu_bd$azimuth))))

cat("\nHeatmap plots of annotations...\n")
p2 <- ggplot(cdata, aes(x=transfer_pc30_azimuth, 
                   y=azimuth_rev, #fill=n/total_n))+
                   fill=log10(n+1)))+
  geom_tile(color="grey90")+
  scale_fill_gradientn(colors=colorRampPalette(c("white","grey70","black"), bias=.5)(6))+
  geom_text(data=mutate(cdata, prop_n=round(prop_n,2)*100) %>% filter(prop_n>15), 
            aes(label=prop_n), color="tomato3", size=3)+
  scale_x_discrete("Label transfer: azimuth", expand = c(0,0))+
  scale_y_discrete("Original azimuth annotation", expand = c(0,0))+
  ggtitle("BD label transfer: PC30")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   aspect.ratio=1, legend.key.size = unit(10, "pt"))

cdata2 = group_by(seu_bd@meta.data, azimuth, transfer_pc30_low.res, .drop=F) %>% tally() %>%
  group_by(azimuth) %>% mutate(total_n=sum(n)) %>%
  ungroup() %>%
  mutate(prop_n=n/total_n,
         azimuth_rev = factor(azimuth, levels=rev(levels(seu_bd$azimuth))))

p3 <- ggplot(cdata2, aes(x=transfer_pc30_low.res, 
                  y=azimuth_rev, #fill=n/total_n))+
                  fill=log10(n+1)))+
  geom_tile(color="grey90")+
  scale_fill_gradientn(colors=colorRampPalette(c("white","grey70","black"), bias=.5)(6))+
  geom_text(data=mutate(cdata2, prop_n=round(prop_n,2)*100) %>% filter(prop_n>15),
            aes(label=prop_n), color="tomato3", size=3)+
  scale_x_discrete("Label transfer: seurat_low.res", expand = c(0,0))+
  scale_y_discrete("Original azimuth annotation", expand = c(0,0))+
  ggtitle("BD label transfer: PC30")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   aspect.ratio=1, legend.key.size = unit(10, "pt"))

pdf(file="plots/05_clustering/Seurat/seu-bd_conservative_label-transfer-pc30_plots.pdf", height=7, width=7)
gridExtra::grid.arrange(p1[[1]]+theme(text=element_text(size=10))+NoLegend(), 
                        p1[[2]]+theme(text=element_text(size=10))+NoLegend(), 
                        p1[[3]]+scale_color_manual(values=low.res.pal)+theme(text=element_text(size=10))+NoLegend(), ncol=2)
p2
p3
dev.off()
cat("\nSaved to: plots/05_clustering/Seurat/seu-bd_conservative_label-transfer-pc30_plots.pdf\n")

## Reproducibility information
print("Reproducibility information:")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
