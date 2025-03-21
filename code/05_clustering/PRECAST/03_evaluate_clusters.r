setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(Seurat)
	library(ggspavis)
})

#loss plot
x = readLines("code/05_clustering/PRECAST/logs/precast_n1079_k9_15278541.log")
iters = x[grep("iter =", x)]
iter.list = strsplit(iters, ", ")
df = data.frame(iter=2:(length(iters)+1), 
           loglik=as.numeric(sapply(iter.list, function(x) substr(x[[2]], start=9, stop=20))),
           d.loglik=as.numeric(sapply(iter.list, function(x) substr(x[[3]], start=9, stop=16))))

png(file="plots/05_clustering/PRECAST_n1079-k9_lglk-loss-plot.png", bg="white")
	plot(df[["iter"]], df[["loglik"]], main="PRECAST n=1079 k=9", xlab="iteration", ylab="loglik")
dev.off()

#load seurat results
load("processed-data/05_clustering/PRECAST/srt_precast_k-9_n1079.Rdata")
#load spe and add seurat key
spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")
spe$seurat_key = paste(spe$sample_id, colnames(spe), sep="_")
#match obs
mdata = seuInt@meta.data[spe$seurat_key,]
stopifnot(identical(rownames(mdata), spe$seurat_key))
spe$precast_k9_1079 = as.factor(mdata$cluster)

quickResaveHDF5SummarizedExperiment(spe)
cat("\n\nPRECAST clusters (n=1079 genes, k=9) updated to spe with quickResave\n\n")

table(colData(spe)[,c("precast_k9_1079","problem_area_flag")])

#plot it out bb!!
uniquepal = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#8B0000","#FF7F00","#CAB2D6")
names(uniquepal) = c("2","7","9","5","1","4","6","3","8")

spe$dummy_slide = ifelse(spe$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "A1", spe$array)

seed = levels(as.factor(spe$dummy_slide))
slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList = lapply(slideList, function(x) {
	do.call(cbind, lapply(x, function(y) spe[,spe$dummy_slide==y]))
})

cat("\nGenerate spot plots...\n")
clusPlot = lapply(slideList, function(x) {
	suppressMessages(
		plotSpots(x, annotate="precast_k9_1079", point_size=.5, sample_id="sample_id")+
		scale_color_manual(values=uniquepal)+
		facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
		theme(panel.background=element_rect(fill="grey30"))
	)
})


pdf(file="plots/05_clustering/PRECAST_n1079-k9_spot-plots.pdf", width=12, height=16)
	clusPlot[[1]]
	clusPlot[[2]]
	clusPlot[[3]]
	clusPlot[[4]]
	clusPlot[[5]]
dev.off()
cat("\nSaved spot plots to: plots/05_clustering/PRECAST_n1079-k9_spot-plots.pdf\n")

cat("\n\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
