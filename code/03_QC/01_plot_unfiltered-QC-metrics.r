setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(ggspavis)
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(dplyr)
	library(ggplot2)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
cat("Dim spe:",dim(spe),"\n")
spe <- spe[,spe$in_tissue]
cat("Dim spe (in tissue):",dim(spe),"\n")

#modify colData for plotting
spe$lg10.umi = log10(spe$sum_umi)
spe$condition = factor(spe$condition, levels=c("NTC","MDD","BPD"))
spe$slide2 = ifelse(spe$slide=="V13B23-283","V13B23-339",spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "B_1", spe$array)
spe$sample_id2 = paste(spe$slide2, spe$array2)

recode.df = as.data.frame(colData(spe)) %>% select(brnum, sample_id2, round, sample_id) %>% distinct() %>% arrange(round, sample_id2)
recode.df$sample_id = ifelse(recode.df$brnum=="Br5366", paste("***",recode.df$sample_id), recode.df$sample_id)
recode.df$sample_id3 = factor(recode.df$sample_id2, levels=recode.df$sample_id2, labels=paste(recode.df$sample_id, recode.df$round))

spe$sample_id2 = factor(spe$sample_id2, levels=recode.df$sample_id2, labels=paste(recode.df$sample_id, recode.df$round))

colData(spe)$box_facets = ""
colData(spe)[spe$sample_id2 %in% recode.df$sample_id3[1:30],"box_facets"] = "g1"
colData(spe)[spe$sample_id2 %in% recode.df$sample_id3[31:60],"box_facets"] = "g2"
colData(spe)[spe$sample_id2 %in% recode.df$sample_id3[61:90],"box_facets"] = "g3"
colData(spe)[spe$sample_id2 %in% recode.df$sample_id3[91:120],"box_facets"] = "g4"

#QC metrics summary
cat("Plotting QC boxplots...",format(Sys.time()),"\n")
p1 <- ggplot(as.data.frame(colData(spe)), aes(x=sample_id2, y=sum_umi, fill=condition, lty=sex))+
  geom_boxplot(outlier.size=.5, linewidth=.5)+geom_hline(aes(yintercept=1000), color="grey50", linewidth=1.5, alpha=.6)+
  scale_fill_manual(values=c("grey50","#9e771b","#1b9e77"))+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(0,50000),
                     breaks=c(10^(0:5)), labels=c("1","10","100","1k","10k","100k"))+
  facet_wrap(vars(box_facets), ncol=1, scales="free_x")+
  labs(x="", y="sum_umi (log10 scale)", title="Library size", fill="", lty="")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom",
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p2 <- ggplot(as.data.frame(colData(spe)), aes(x=sample_id2, y=sum_gene, fill=condition, lty=sex))+
  geom_boxplot(outlier.size=.5, linewidth=.5)+
  scale_fill_manual(values=c("grey50","#9e771b","#1b9e77"))+
  facet_wrap(vars(box_facets), ncol=1, scales="free_x")+
  labs(x="", y="sum_gene", title="Detected genes", fill="", lty="")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom",
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p3 <- ggplot(as.data.frame(colData(spe)), aes(x=sample_id2, y=expr_chrM_ratio, fill=condition, lty=sex))+
  geom_boxplot(outlier.size=.5, linewidth=.5)+
  scale_fill_manual(values=c("grey50","#9e771b","#1b9e77"))+
  facet_wrap(vars(box_facets), ncol=1, scales="free_x")+
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
#split samples into 5 lists of 24 (5 pages of 6 slides each)
slideList = list(recode.df$sample_id3[1:24], recode.df$sample_id3[25:48], recode.df$sample_id3[49:72],
                 recode.df$sample_id3[73:96], recode.df$sample_id3[97:120])

slideList = lapply(slideList, function(x) {
  spe_sub = spe[,spe$sample_id2 %in% x]
  #modify spatialCoords so all capture areas start at 0
  mod_spatialCoords = spatialCoords(spe_sub)
  for (i in unique(spe_sub$sample_id)) {
    tmp = mod_spatialCoords[colData(spe_sub)$sample_id==i,]
    mod_spatialCoords[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
    mod_spatialCoords[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
  }
  spatialCoords(spe_sub) <- mod_spatialCoords
  return(spe_sub)
})

cat("\nPlotting library size...",format(Sys.time()),"\n")
libList = lapply(slideList, function(x) {
	suppressMessages(plotSpots(x, annotate="lg10.umi", point_size=0.4, sample_id="sample_id2")+
		scale_color_gradient("lg10.umi", low="white", high="navy")+#, labels=function(y) paste0(y/1000,"k"))+
		facet_wrap(vars(sample_id2), ncol=4)+
		theme(plot.title=element_blank(),
			strip.background = element_rect(fill="transparent", color="transparent"),
			panel.background=element_rect(fill="grey30"))
	)
})

cat("Compiling library size plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_library-size_spot-plots.pdf", width=12, height=16)
	libList[[1]]
	libList[[2]]
	libList[[3]]
	libList[[4]]
	libList[[5]]
dev.off()
cat("Saved to: plots/03_QC/unfiltered_library-size_spot-plots.pdf\n")

cat("\nPlotting # genes detected...",format(Sys.time()),"\n")
geneList = lapply(slideList, function(x) {
	suppressMessages(plotSpots(x, annotate="sum_gene", point_size=0.4, sample_id="sample_id2")+
		scale_color_gradient("sum_genes", low="white", high="navy", labels=function(y) paste0(y/1000,"k"))+
		facet_wrap(vars(sample_id2), ncol=4)+
		theme(plot.title=element_blank(),
                        strip.background = element_rect(fill="transparent", color="transparent"),
			panel.background=element_rect(fill="grey30"))
	)
})

cat("Compiling n genes plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_n-genes_spot-plots.pdf", width=12, height=16)
	geneList[[1]]
	geneList[[2]]
	geneList[[3]]
	geneList[[4]]
	geneList[[5]]
dev.off()
cat("Saved to: plots/03_QC/unfiltered_n-genes_spot-plots.pdf\n")

cat("\nPlotting chrM ratio",format(Sys.time()),"\n")
cat("*** Max color limit fixed across all samples to 0.65 to help with vis of samples where some spots have 100% chrM ratio\n")
mitoList = lapply(slideList, function(x) {
	suppressMessages(plotSpots(x, annotate="expr_chrM_ratio", point_size=0.4, sample_id="sample_id2")+
		scale_color_gradient("mito\nfraction", low="white", high="navy", limits=c(0,0.65))+
		facet_wrap(vars(sample_id2), ncol=4)+
		theme(plot.title=element_blank(), 
                        strip.background = element_rect(fill="transparent", color="transparent"),
			panel.background=element_rect(fill="grey30"))
	)
})

cat("Compiling chrM ratio plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/unfiltered_chrM-ratio_spot-plots.pdf", width=12, height=16)
	mitoList[[1]]
	mitoList[[2]]
	mitoList[[3]]
	mitoList[[4]]
	mitoList[[5]]
dev.off()
cat("Saved to: plots/03_QC/unfiltered_chrM-ratio_spot-plots.pdf\n")


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
