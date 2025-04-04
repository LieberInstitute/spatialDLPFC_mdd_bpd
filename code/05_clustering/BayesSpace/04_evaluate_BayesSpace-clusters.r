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
cdata = colData(spe)[,c("sample_id","slide","array","brnum","age","PMI","RIN")]

########## compile BayesSpace results
bsp.colorList = list()
#iter= 10K
cat("\nLoad 10k results...\n")
d1 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d20_q9_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d1)))
cdata$BayesSpace_PCA.1079.d20_q9_iter.10k_init.kmeans = d1$BayesSpace_PCA.1079.d20_q9_init.kmeans
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")
  
d1.2 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d13_q9_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d1.2)))
cdata$BayesSpace_PCA.1079.d13_q9_iter.10k_init.kmeans = d1.2$BayesSpace_PCA.1079.d13_q9_init.kmeans
bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.10k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.10k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")

d2 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d20_q9_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d2)))
cdata$BayesSpace_PCA.1079.d20_q9_iter.10k_init.PRECAST = d2$BayesSpace_PCA.1079.d20_q9_init.PRECAST
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.10k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")


#init clus= kmeans
cat("\nLoad initial cluster= kmeans results...\n")
d3 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d20_q9_iter-5k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d3)))
cdata$BayesSpace_PCA.1079.d20_q9_iter.5k_init.kmeans = d3$BayesSpace_PCA.1079.d20_q9_iter.5k_init.kmeans
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.5k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.5k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")

d3.2 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d13_q9_iter-5k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d3.2)))
cdata$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans = d3.2$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans
bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans) = c("4","1","2","5","9","3","7","6","8")

d3.3 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d12_q9_iter-5k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d3.3)))
cdata$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans = d3.3$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#A6CEE3","#6A3D9A")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans) = c("6","4","2","5","9","3","7","8","1")


#init clus= precast n1079
cat("\nLoad initial cluster= PRECAST n1079 results...\n")
d4 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d20_q9_iter-1k_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d4)))
cdata$BayesSpace_PCA.1079.d20_q9_iter.1k_init.PRECAST <- d4$BayesSpace_PCA.1079.d20_q9_iter.1k_init.PRECAST
bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.1k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d20_q9_iter.1k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")

d4.2 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d13_q9_iter-1k_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d4.2)))
cdata$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST <- d4.2$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST
bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")

d4.3 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d12_q9_iter-1k_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d4.3)))
cdata$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST <- d4.3$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")

d4.4 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d12_q9_iter-5k_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d4.4)))
cdata$BayesSpace_PCA.1079.d12_q9_iter.5k_init.PRECAST = d4.4$BayesSpace_PCA.1079.d12_q9_iter.5k_init.PRECAST
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.5k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")

d4.5 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d12_q9_iter-10k_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d4.5)))
cdata$BayesSpace_PCA.1079.d12_q9_iter.10k_init.PRECAST = d4.5$BayesSpace_PCA.1079.d12_q9_iter.10k_init.PRECAST
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.10k_init.PRECAST = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.10k_init.PRECAST) = c("7","8","9","5","1","6","4","3","2")


#init clus = precast n1663
cat("\nLoad intial cluster= PRECAST n1663 results...\n")
d7 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-1079-d12_q9_iter-1k_init-PRECAST-1663.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d7)))
cdata$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST.1663 = d7$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST.1663
bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST.1663 = c("#FF7F00","#1F78B4","navy","#33A02C","#A6CEE3","#6A3D9A","#B2DF8A","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST.1663) = c("5","3","4","2","7","8","1","9","6")


#reduced dim= MNN
cat("\nLoad reduced dims= MNN results...\n")
d5 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_MNN-1079-d13_q9_iter-5k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d5)))
cdata$BayesSpace_MNN.1079.d13_q9_iter.5k_init.kmeans = d5$BayesSpace_MNN.1079.d13_q9_iter.5k_init.kmeans
bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.5k_init.kmeans = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#FDBF6F","#B2DF8A","#FB9A99","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.5k_init.kmeans) <- c("4","6","5","8","1","2","9","3","7")

d5.2 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_MNN-1079-d13_q9_iter-1k_init-PRECAST.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d5.2)))
cdata$BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST = d5.2$BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST
bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#FDBF6F","#B2DF8A","#FB9A99","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_MNN.1079.d13_q9_iter.1k_init.PRECAST) <- c("3","7","5","8","2","9","1","6","4")


#reduced dim= PRECAST embeddings
cat("\nLoad reduced dim= PRECAST embeddings results...\n")
d8 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PRECAST-1663-d15_q9_iter-1k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d8)))
cdata$BayesSpace_PRECAST.1663.d15_q9_iter.1k_init.kmeans = d8$BayesSpace_PRECAST.1663.d15_q9_iter.1k_init.kmeans
bsp.colorList$BayesSpace_PRECAST.1663.d15_q9_iter.1k_init.kmeans = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#6A3D9A","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C")
names(bsp.colorList$BayesSpace_PRECAST.1663.d15_q9_iter.1k_init.kmeans) = c("8","5","6","2","1","9","7","4","3")

d8.2 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PRECAST-1079-d15_q9_iter-1k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d8.2)))
cdata$BayesSpace_PRECAST.1079.d15_q9_iter.1k_init.kmeans = d8.2$BayesSpace_PRECAST.1079.d15_q9_iter.1k_init.kmeans
bsp.colorList$BayesSpace_PRECAST.1079.d15_q9_iter.1k_init.kmeans = c("#FF7F00","#1F78B4","#33A02C","#A6CEE3","#6A3D9A","#CAB2D6","#B2DF8A","#FB9A99","#E31A1C")
names(bsp.colorList$BayesSpace_PRECAST.1079.d15_q9_iter.1k_init.kmeans) = c("8","6","5","2","1","9","7","4","3")


#HVGs
cat("\nLoad HVGs, BSP PCA results...\n") 
d6 = read.csv("processed-data/05_clustering/BayesSpace/colData_BayesSpace_PCA-2k-HVGs-d15_q9_iter-5k_init-kmeans.csv", row.names=1)
stopifnot(identical(rownames(cdata), rownames(d6)))
cdata$BayesSpace_PCA.2k.HVGs.d15_q9_iter.5k_init.kmeans = d6$BayesSpace_PCA.2k.HVGs.d15_q9_iter.5k_init.kmeans
bsp.colorList$BayesSpace_PCA.2k.HVGs.d15_q9_iter.5k_init.kmeans = c("#FF7F00","#1F78B4","navy","#33A02C","#A6CEE3","#6A3D9A","#B2DF8A","#E31A1C","#8B0000")
names(bsp.colorList$BayesSpace_PCA.2k.HVGs.d15_q9_iter.5k_init.kmeans) = c("6","4","7","5","1","8","2","3","9")

#save bsp results
write.csv(cdata, "processed-data/05_clustering/BayesSpace/colData_all-BayesSpace-clusters.csv", row.names=T)
cat("\nCombined results saved to: processed-data/05_clustering/BayesSpace/colData_all-BayesSpace-clusters.csv\n")

#good sample demonstration
#sub_samples = c("V13B23-309_A1","V13B23-301_A1") #309 = clean, 301 = dirty
#cdata_sub = cdata[cdata$sample_id %in% sub_samples,]
#spe_sub = spe[,rownames(cdata_sub)]

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


####################################
## saved full results (done in interactive session)
####################################
#
##"BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans"
#uniquepal = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
#names(uniquepal) = c("4","1","2","5","9","3","7","6","8")
#spe$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans = as.factor(spe$BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans)
#
##"BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans
#uniquepal = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#A6CEE3","#6A3D9A")
#names(uniquepal) = c("6","4","2","5","9","3","7","8","1")
#spe$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans = as.factor(spe$BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans)
#
##"BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST"
#spe$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST = as.factor(spe$BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST)
#uniquepal = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
#names(uniquepal) = c("7","8","9","5","1","6","4","3","2")
#
##"BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST"
#spe$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST = as.factor(spe$BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST)
#uniquepal = c("#1F78B4","navy","#B2DF8A","#33A02C","#FB9A99","#8B0000","#E31A1C","#FF7F00","#A6CEE3")
#names(uniquepal) = c("7","8","9","5","1","6","4","3","2")
#
#
#spe$dummy_slide = ifelse(spe$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe$slide)
#spe$array2 = ifelse(spe$slide=="V13B23-283", "A1", spe$array)
#spe$lg10_sum_umi = log10(spe$sum_umi)
#seed = levels(as.factor(spe$dummy_slide))
#slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])
#slideList = lapply(slideList, function(x) {
#  do.call(cbind, lapply(x, function(y) spe[,spe$dummy_slide==y]))
#})
#
#cat("\nGenerate spot plots...\n")
#clusPlot = lapply(slideList, function(x) {
#  suppressMessages(
#    plotSpots(x, #annotate="BayesSpace_PCA.1079.d13_q9_iter.5k_init.kmeans",
#              #annotate="BayesSpace_PCA.1079.d12_q9_iter.5k_init.kmeans",
#              #annotate="BayesSpace_PCA.1079.d13_q9_iter.1k_init.PRECAST",
#              annotate="BayesSpace_PCA.1079.d12_q9_iter.1k_init.PRECAST",
#              #annotate="lg10_sum_umi",
#              point_size=.5, sample_id="sample_id")+
#      scale_color_manual("BayesSpace\nPCA.1079.d12\nq9 iter.1k\ninit.PRECAST", 
#                         #"BayesSpace\nPCA.1079.d13\nq9 iter.1k\ninit.PRECAST", 
#                         #"BayesSpace\nPCA.1079.d12\nq9 iter.5k\ninit.kmeans", 
#                         #"BayesSpace\nPCA.1079.d13\nq9 iter.5k\ninit.kmeans", 
#                         values=uniquepal)+
#      #scale_color_gradient("lg10(sum_umi)",low="white", high="navy")+#, 
#      #labels=function(x) paste0(x/1000,"k"))+
#      facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
#      theme(panel.background=element_rect(fill="grey30"))
#  )
#})
#
##pdf(file="plots/05_clustering/BayesSpace_PCA-1079-d13_q9_iter-5k_init-kmean_spot-plots.pdf", width=12, height=16)
##pdf(file="plots/05_clustering/BayesSpace_PCA-1079-d12_q9_iter-5k_init-kmean_spot-plots.pdf", width=12, height=16)
##pdf(file="plots/05_clustering/BayesSpace_PCA-1079-d13_q9_iter-1k_init-PRECAST_spot-plots.pdf", width=12, height=16)
#pdf(file="plots/05_clustering/BayesSpace_PCA-1079-d12_q9_iter-1k_init-PRECAST_spot-plots.pdf", width=12, height=16)
##pdf(file="plots/05_clustering/lg10-sum-umi_spot-plots.pdf", width=12, height=16)
#clusPlot[[1]]
#clusPlot[[2]]
#clusPlot[[3]]
#clusPlot[[4]]
#clusPlot[[5]]
#dev.off()
#
####################################
####################################


cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
