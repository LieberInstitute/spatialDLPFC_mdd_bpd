setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(bluster)
	library(pheatmap)
})
set.seed(123)

# azimuth to azimuth for ground truth
bd.transfer = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-azimuth.csv", row.names=1)

# load all SZBDMulti-seq obs data
cdata = read.csv("processed-data/05_clustering/Seurat/obs_SZBDMulti-Seq_filtered.csv")
# subset to BD that were transfered (QC filtered)
cdata2 = cdata[cdata$barcodekey %in% rownames(bd.transfer),]
rownames(cdata2) <- cdata2$barcodekey
cdata2 = cdata2[rownames(bd.transfer),]

# pairwise rand
check1 = pairwiseRand(cdata2$azimuth, bd.transfer$predicted.id, mode="ratio")
rand1 = pairwiseRand(cdata2$azimuth, bd.transfer$predicted.id, mode="index")

order1 = c("OPC","Oligo","Astro","Endo","SMC","PC","VLMC","Immune","Micro",
	"L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 IT","L6 IT Car3","L6b","L6 CT",
	"Lamp5","Lamp5 Lhx6","Vip","Pax6","Sncg","Sst Chodl","Chandelier","Pvalb","Sst")

# plot prediction score max
bd.transfer$true.label = cdata2$azimuth
avg.pred = group_by(bd.transfer[,c(2:28,30)], true.label) %>% summarise_all(mean)
pred1 = as.matrix(avg.pred[,-1])
rownames(pred1) <- avg.pred$true.label

phm1 = pheatmap(pred1[order1, paste0("prediction.score.", gsub("/","\\.", gsub(" ", "\\.", order1)))], 
	color = colorRampPalette(c("white",RColorBrewer::brewer.pal(n=5, "Purples")[2:5]))(100),
	breaks= seq(from = 0, to = 1, length.out = 101),
	cluster_rows=F, cluster_col=F, angle_col=90, fontsize=7,
	main=paste("BD snRNA-seq transfer Rand index:", round(rand1,4)))
#clustering_method="ward.D2")


# now the low res version that I used for SRT
bd.transfer2 = read.csv("processed-data/05_clustering/Seurat/results_label-transfer_MBv-filtered_ref-control_query-bipolar_qual-genes-kanchor-50-pc30_red-pca-kweight-50-low-res.csv", row.names=1)
stopifnot(identical(rownames(bd.transfer2), rownames(cdata2)))

rand2 = pairwiseRand(cdata2$azimuth, bd.transfer2$predicted.id, mode="index")

bd.transfer2$true.label = cdata2$azimuth
avg.pred = group_by(bd.transfer2[,c(2:10,12)], true.label) %>% summarise_all(mean)
pred2 = as.matrix(avg.pred[,-1])
rownames(pred2)	<- avg.pred$true.label

phm2 = pheatmap(pred2[order1,paste0("prediction.score.",c("Oligo","Astro","Micro.Vasc","L2","L3","L4","L5","L6","Inhb"))], 
	color = colorRampPalette(c("white",RColorBrewer::brewer.pal(n=5, "Purples")[2:5]))(100),
        breaks= seq(from = 0, to = 1, length.out = 101),
	cluster_col=F, cluster_rows=F, angle_col=90, fontsize=7,
	main=paste("BD snRNA-seq transfer Rand index:", round(rand2,4)))

pdf(file="plots/publication/Figure3/supp_BD-transfer-validation.pdf", height=4, width=4)
plot(phm1[[4]])
plot(phm2[[4]])
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

