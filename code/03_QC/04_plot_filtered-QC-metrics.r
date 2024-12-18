setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(ggspavis)
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(gridExtra)
	library(scater)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/03_QC/", prefix="spe_n120_postQC_")
spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

spe$sex = ifelse(spe$sex=="","M",spe$sex)

#QC metrics summary
p1 <- ggcells(spe, aes(x=array, y=sum_umi, fill=condition, lty=sex))+
  geom_violin()+facet_wrap(vars(dummy_slide))+
  scale_y_log10()+theme_bw()+ggtitle("Library size")

p2 <- ggcells(spe, aes(x=array, y=sum_gene, fill=condition, lty=sex))+
  geom_violin()+facet_wrap(vars(dummy_slide))+
  theme_bw()+ggtitle("Detected genes")
  
p3 <- ggcells(spe, aes(x=array, y=expr_chrM_ratio, fill=condition, lty=sex))+
  geom_violin()+facet_wrap(vars(dummy_slide))+
  theme_bw()+ggtitle("Mitochondrial fraction")

cat("Compiling QC summary plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/QC-filtered_qc-metrics_violin-plots.pdf", width=8.5, height=11)
p1
p2
p3
dev.off()

#QC metrics spot plots
dummy = c("V13B23-283","V13B23-339")
remainder = setdiff(unique(spe$slide), dummy)
slideList = list(c(dummy, remainder[1:5]),remainder[6:11],remainder[12:17],remainder[18:23],remainder[24:29])
  
slideList = lapply(slideList, function(x)
	unlist(lapply(x, function(y)
	unique(colData(spe)[spe$slide==y,"sample_id"])
	))
)

libList = lapply(slideList, function(x) {
	cat("Plotting library size...",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$sample_id==y])
  
	lapply(1:length(l1), function(z)
		suppressMessages(plotSpots(l1[[z]], annotate="sum_umi", point_size=0.1)+
		scale_color_gradient(low="grey90", high="black", labels=function(x) paste0(x/1000,"k"))+
		labs(title=names(l1)[[z]])+
		theme(legend.text=element_text(size=8))
	))
})

cat("Compiling library size plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/QC-filtered_library-size_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(libList[[1]], ncol=4))
do.call(grid.arrange, c(libList[[2]], ncol=4))
do.call(grid.arrange, c(libList[[3]], ncol=4))
do.call(grid.arrange, c(libList[[4]], ncol=4))
do.call(grid.arrange, c(libList[[5]], ncol=4))
dev.off()


geneList = lapply(slideList, function(x) {
	cat("Plotting genes...",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$sample_id==y])
  
	lapply(1:length(l1), function(z)
		suppressMessages(plotSpots(l1[[z]], annotate="sum_gene", point_size=0.1)+
		scale_color_gradient(low="grey90", high="black", labels=function(x) paste0(x/1000,"k"))+
		labs(title=names(l1)[[z]])+
		theme(legend.text=element_text(size=8))
	))
})

cat("Compiling n genes plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/QC-filtered_n-genes_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(geneList[[1]], ncol=4))
do.call(grid.arrange, c(geneList[[2]], ncol=4))
do.call(grid.arrange, c(geneList[[3]], ncol=4))
do.call(grid.arrange, c(geneList[[4]], ncol=4))
do.call(grid.arrange, c(geneList[[5]], ncol=4))
dev.off()


mitoList = lapply(slideList, function(x) {
	cat("Plotting chrM ratio",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$sample_id==y])
  
	lapply(1:length(l1), function(z)
		suppressMessages(plotSpots(l1[[z]], annotate="expr_chrM_ratio", point_size=0.1)+
		scale_color_gradient(low="grey90", high="black")+
		labs(title=names(l1)[[z]])+
		theme(legend.text=element_text(size=8))
	))
})

cat("Compiling chrM ratio plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/QC-filtered_chrM-ratio_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(mitoList[[1]], ncol=4))
do.call(grid.arrange, c(mitoList[[2]], ncol=4))
do.call(grid.arrange, c(mitoList[[3]], ncol=4))
do.call(grid.arrange, c(mitoList[[4]], ncol=4))
do.call(grid.arrange, c(mitoList[[5]], ncol=4))
dev.off()

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
