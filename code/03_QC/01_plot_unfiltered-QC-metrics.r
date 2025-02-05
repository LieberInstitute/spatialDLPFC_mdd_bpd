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

colData(spe)$facet_violin = paste(spe$round, spe$dummy_slide)
colData(spe)$facet_violin = ifelse(spe$dummy_slide=="joint-283-339", "joint-283-339", spe$facet_violin)

colData(spe)$facet_spots = paste(spe$round, spe$sample_id)
colData(spe)$facet_spots = ifelse(spe$brnum=="Br5366", paste(spe$facet_spots, spe$brnum), spe$facet_spots)

x_ordered = sort(unique(spe$facet_spots))
colData(spe)$x_facets = ""
colData(spe)[spe$facet_spots %in% x_ordered[1:30],"x_facets"] = "g1"
colData(spe)[spe$facet_spots %in% x_ordered[31:60],"x_facets"] = "g2"
colData(spe)[spe$facet_spots %in% x_ordered[61:90],"x_facets"] = "g3"
colData(spe)[spe$facet_spots %in% x_ordered[91:120],"x_facets"] = "g4"

#QC metrics summary
cat("Plotting QC boxplots...",format(Sys.time()),"\n")
p1 <- ggplot(as.data.frame(colData(spe)), aes(x=facet_spots, y=sum_umi, fill=condition, lty=sex))+
  geom_boxplot(outlier.size=.5, linewidth=.5)+geom_hline(aes(yintercept=1000), color="grey50", linewidth=1.5, alpha=.6)+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(0,50000),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  labs(x="", y="sum_umi (log10 scale)", title="Library size", fill="", lty="")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom",
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p2 <- ggplot(as.data.frame(colData(spe)), aes(x=facet_spots, y=sum_gene, fill=condition, lty=sex))+
  geom_boxplot(outlier.size=.5, linewidth=.5)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  labs(x="", y="sum_gene", title="Detected genes", fill="", lty="")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom",
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p3 <- ggplot(as.data.frame(colData(spe)), aes(x=facet_spots, y=expr_chrM_ratio, fill=condition, lty=sex))+
  geom_boxplot(outlier.size=.5, linewidth=.5)+
  facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
  labs(x="", y="expr_chrM_ratio", title="Mitochondrial fraction", fill="", lty="")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom",
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())


pdf(file="plots/03_QC/unfiltered_qc-metrics_boxplot.pdf", width=9, height=12)
p1
p2
p3
dev.off()
cat("Saved to: plots/03_QC/unfiltered_qc-metrics_boxplot.pdf\n")

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


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
