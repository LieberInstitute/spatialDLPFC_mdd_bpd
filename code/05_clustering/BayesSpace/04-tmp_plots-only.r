setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(ggspavis)
        library(gridExtra)
        library(mclust)
        library(ggplot2)
})
Sys.time()
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")

cdata = read.csv("processed-data/05_clustering/BayesSpace/colData_all-BayesSpace-clusters.csv", row.names=1)


bsp.colorList = list()
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")
bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.10k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.10k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.5k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.5k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")
bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#A6CEE3","#6A3D9A")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans) = c("6","4","2","5","9","3","7","8","1")
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.1k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.1k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")
bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.10k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.10k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST.1663 = c("#FF7F00","#1F78B4","navy","#33A02C","#A6CEE3","#6A3D9A","#B2DF8A","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST.1663) = c("5","3","4","2","7","8","1","9","6")
bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.5k_init.kmeans = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#FDBF6F","#B2DF8A","#FB9A99","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.5k_init.kmeans) <- c("4","6","5","8","1","2","9","3","7")
bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#FDBF6F","#B2DF8A","#FB9A99","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST) <- c("3","7","5","8","2","9","1","6","4")
bsp.colorList$BayesSpace_PRECAST.1663.d15_q9_iter.1k_init.kmeans = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#6A3D9A","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C")
names(bsp.colorList$BayesSpace_PRECAST.1663.d15_q9_iter.1k_init.kmeans) = c("8","5","6","2","1","9","7","4","3")
bsp.colorList$BayesSpace_PRECAST.1079.d15_q9_iter.1k_init.kmeans = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#6A3D9A","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C")
names(bsp.colorList$BayesSpace_PRECAST.1079.d15_q9_iter.1k_init.kmeans) = c("8","6","5","2","1","9","7","4","3")
bsp.colorList$BayesSpace_PCA.2k.HVGs.d15_q9_iter.5k_init.kmeans = c("#FF7F00","#1F78B4","navy","#33A02C","#A6CEE3","#6A3D9A","#B2DF8A","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_PCA.2k.HVGs.d15_q9_iter.5k_init.kmeans) = c("6","4","7","5","1","8","2","3","9")

stopifnot(length(bsp.colorList)==17)

#four sample demonstration
(n16.samples = paste0("V13B23-30",c(1,2,8,9)))
spe_sub = spe[,spe$slide %in% n16.samples]
cdata_sub = cdata[cdata$slide %in% n16.samples,]

plist = lapply(names(bsp.colorList), function(x) {
  colData(spe_sub)[[x]] = as.factor(cdata_sub[,x])
  suppressMessages(plotSpots(spe_sub, sample_id="sample_id", annotate=x,
                             point_size=.3)+
                     scale_color_manual("",values=bsp.colorList[[x]])+
                     labs(title=paste0(x," --- ARI with PRECAST n1079= ",round(adjustedRandIndex(spe$precast_k9_1079, cdata[[x]]),3)))+
                     theme(strip.background = element_rect(fill=NA, color=NA),
                           legend.key.size=unit(6,"pt"),
                           legend.title=element_text(margin=margin(0,0,4,0,"pt")),
                           legend.box.spacing = unit(2,"pt"),
                           legend.margin=margin(0,0,0,0,"pt"),
                           legend.box.margin = margin(0,2,0,2,"pt"),
                           plot.title=element_text(size=14), #plot.subtitle = element_text(size=8),
                           plot.margin = unit(c(3,3,3,3),"pt"))
                   )
})
pdf("plots/05_clustering/BayesSpace/BayesSpace_all-results_16-example_spot-plots.pdf", width=10, height=10)
for(i in 1:length(plist)) {
        plot(plist[[i]])
}
dev.off()
cat("\nExample plots saved to: processed-data/05_clustering/BayesSpace/BayesSpace_all-results_16-example_spot-plots.pdf\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
