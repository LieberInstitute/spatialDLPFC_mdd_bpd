setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(ggspavis)
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(gridExtra)
	library(scater)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Dim spe:",dim(spe),"\n")
spe <- spe[,spe$in_tissue]
cat("Dim spe (in tissue):",dim(spe),"\n")

spe$dummy_slide = spe$slide
colData(spe)[spe$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe)[spe$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

#QC metrics summary
colData(spe)$facet_violin = paste(spe$round, spe$dummy_slide)
colData(spe)$facet_violin = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", spe$facet_violin)
colData(spe)$x_lab = ifelse(spe$slide=="V13B23-283", "283 B1", spe$array)

p1 <- ggcells(spe, aes(x=x_lab, y=sum_umi, fill=condition, lty=sex))+
  geom_violin()+facet_wrap(vars(facet_violin), scales="free_x")+
  scale_y_log10()+theme_bw()+ggtitle("Library size")

p2 <- ggcells(spe, aes(x=x_lab, y=sum_gene, fill=condition, lty=sex))+
  geom_violin()+facet_wrap(vars(facet_violin), scales="free_x")+
  theme_bw()+ggtitle("Detected genes")
  
p3 <- ggcells(spe, aes(x=x_lab, y=expr_chrM_ratio, fill=condition, lty=sex))+
  geom_violin()+facet_wrap(vars(facet_violin), scales="free_x")+
  theme_bw()+ggtitle("Mitochondrial fraction")

cat("Compiling QC summary plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_qc-metrics_violin-plots.pdf", width=8.5, height=11)
p1
p2
p3
dev.off()
cat("Saved to: plots/03_QC/unfiltered_qc-metrics_violin-plots.pdf\n")

#QC metrics spot plots
colData(spe)$facet_spots = paste(spe$round, spe$sample_id)
colData(spe)$facet_spots = ifelse(spe$brnum=="Br5366", paste(spe$facet_spots, spe$brnum), spe$facet_spots)

seed = levels(as.factor(spe$facet_violin))
slideList = list(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
  
slideList = lapply(slideList, function(x)
	unlist(lapply(x, function(y)
	sort(unique(colData(spe)[spe$facet_violin==y,"facet_spots"]))
	))
)

cat("\nPlotting library size...",format(Sys.time()),"\n")
libList = lapply(slideList, function(x) {
	#cat("Plotting library size...",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
	lapply(1:length(l1), function(z)
		suppressMessages(plotSpots(l1[[z]], annotate="sum_umi", point_size=0.2)+
		scale_color_gradient(low="white", high="navy", labels=function(x) paste0(x/1000,"k"))+
		labs(title=names(l1)[[z]])+
		theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
	))
})

cat("Compiling library size plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_library-size_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(libList[[1]], ncol=4))
do.call(grid.arrange, c(libList[[2]], ncol=4))
do.call(grid.arrange, c(libList[[3]], ncol=4))
do.call(grid.arrange, c(libList[[4]], ncol=4))
do.call(grid.arrange, c(libList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/unfiltered_library-size_spot-plots.pdf\n")

cat("\nPlotting # genes detected...",format(Sys.time()),"\n")
geneList = lapply(slideList, function(x) {
	#cat("Plotting genes...",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
	lapply(1:length(l1), function(z)
		suppressMessages(plotSpots(l1[[z]], annotate="sum_gene", point_size=0.2)+
		scale_color_gradient(low="white", high="navy", labels=function(x) paste0(x/1000,"k"))+
		labs(title=names(l1)[[z]])+
		theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
	))
})

cat("Compiling n genes plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_n-genes_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(geneList[[1]], ncol=4))
do.call(grid.arrange, c(geneList[[2]], ncol=4))
do.call(grid.arrange, c(geneList[[3]], ncol=4))
do.call(grid.arrange, c(geneList[[4]], ncol=4))
do.call(grid.arrange, c(geneList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/unfiltered_n-genes_spot-plots.pdf\n")

cat("\nPlotting chrM ratio",format(Sys.time()),"\n")
cat("*** Max color limit set to second highest expr_chrM_ratio per sample to help with 100% chrM ratio spots\n")
mitoList = lapply(slideList, function(x) {
	#cat("Plotting chrM ratio",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])
  
	lapply(1:length(l1), function(z) {
		max_value = sort(unique(colData(l1[[z]])$expr_chrM_ratio), decreasing=T)[2]
		suppressMessages(plotSpots(l1[[z]], annotate="expr_chrM_ratio", point_size=0.2)+
		scale_color_gradient(low="white", high="navy", limits=c(0,max_value))+
		labs(title=names(l1)[[z]])+
		theme(legend.text=element_text(size=8), panel.background=element_rect(fill="grey30"))
	)})
})

cat("Compiling chrM ratio plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_chrM-ratio_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(mitoList[[1]], ncol=4))
do.call(grid.arrange, c(mitoList[[2]], ncol=4))
do.call(grid.arrange, c(mitoList[[3]], ncol=4))
do.call(grid.arrange, c(mitoList[[4]], ncol=4))
do.call(grid.arrange, c(mitoList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/unfiltered_chrM-ratio_spot-plots.pdf\n")

#plot marker genes
cat("\nPlotting MBP...",format(Sys.time()),"\n")
markerList = lapply(slideList, function(x) {
        #cat("Plotting MBP...",format(Sys.time()),"\n")
        l1 = x; names(l1) = x
        l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

        lapply(1:length(l1), function(z)
                suppressMessages(plotSpots(l1[[z]], annotate="MBP", point_size=0.2, feature_names="gene_name", assay_name="counts")+
                        scale_color_gradient(low="white", high="navy")+
                        labs(title=names(l1)[[z]], color="MBP")+
                        theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                                panel.background=element_rect(fill="grey30"))
        ))
})

cat("Compiling MBP plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_MBP-raw-counts_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(markerList[[1]], ncol=4))
do.call(grid.arrange, c(markerList[[2]], ncol=4))
do.call(grid.arrange, c(markerList[[3]], ncol=4))
do.call(grid.arrange, c(markerList[[4]], ncol=4))
do.call(grid.arrange, c(markerList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/unfiltered_MBP-raw-counts_spot-plots.pdf\n")

cat("\nPlotting GAPDH...",format(Sys.time()),"\n")
markerList = lapply(slideList, function(x) {
        #cat("Plotting GAPDH...",format(Sys.time()),"\n")
        l1 = x; names(l1) = x
        l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

        lapply(1:length(l1), function(z)
                suppressMessages(plotSpots(l1[[z]], annotate="GAPDH", point_size=0.2, feature_names="gene_name", assay_name="counts")+
                        scale_color_gradient(low="white", high="navy")+
                        labs(title=names(l1)[[z]], color="GAPDH")+
                        theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                                panel.background=element_rect(fill="grey30"))
        ))
})

cat("Compiling GAPDH plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_GAPDH-raw-counts_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(markerList[[1]], ncol=4))
do.call(grid.arrange, c(markerList[[2]], ncol=4))
do.call(grid.arrange, c(markerList[[3]], ncol=4))
do.call(grid.arrange, c(markerList[[4]], ncol=4))
do.call(grid.arrange, c(markerList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/unfiltered_GAPDH-raw-counts_spot-plots.pdf\n")

cat("\nPlotting SYT1...",format(Sys.time()),"\n")
markerList = lapply(slideList, function(x) {
        #cat("Plotting SYT1...",format(Sys.time()),"\n")
        l1 = x; names(l1) = x
        l1 = lapply(l1, function(y) spe[,colData(spe)$facet_spots==y])

        lapply(1:length(l1), function(z)
                suppressMessages(plotSpots(l1[[z]], annotate="SYT1", point_size=0.2, feature_names="gene_name", assay_name="counts")+
                        scale_color_gradient(low="white", high="navy")+
                        labs(title=names(l1)[[z]], color="SYT1")+
                        theme(legend.text=element_text(size=8), legend.title=element_text(size=10),
                                panel.background=element_rect(fill="grey30"))
        ))
})

cat("Compiling SYT1 plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_SYT1-raw-counts_spot-plots.pdf", width=12, height=16)
do.call(grid.arrange, c(markerList[[1]], ncol=4))
do.call(grid.arrange, c(markerList[[2]], ncol=4))
do.call(grid.arrange, c(markerList[[3]], ncol=4))
do.call(grid.arrange, c(markerList[[4]], ncol=4))
do.call(grid.arrange, c(markerList[[5]], ncol=4))
dev.off()
cat("Saved to: plots/03_QC/unfiltered_SYT1-raw-counts_spot-plots.pdf\n")


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
