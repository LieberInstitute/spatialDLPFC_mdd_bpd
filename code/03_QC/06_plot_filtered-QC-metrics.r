setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(ggspavis)
  library(dplyr)
  library(gridExtra)
  library(scater)
})

#final cdata
cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv", row.names=1)
cdata$spotsweeper_outlier = cdata$umi_local.outlier | cdata$genes_local.outlier | cdata$chrM.ratio_local.outlier
cdata[is.na(cdata$spotsweeper_outlier),"spotsweeper_outlier"] = FALSE

cdata$remove_spots = cdata$in_tissue==FALSE | cdata$remove_problem.areas | cdata$spotsweeper_outlier
#cdata$problem_area_flag = ifelse(cdata$problem_areas_genes.size>20, TRUE, FALSE)
tmp = filter(cdata, in_tissue==TRUE, true_edges==FALSE) %>% #mutate(lowumi = sum_umi<=100) %>%
  group_by(problem_areas_genes.id) %>%
  summarise(n_lowumi=sum(lowumi), n_spots=n(), prop_lowumi=n_lowumi/n_spots) %>%
  filter(n_spots>20, !is.na(problem_areas_genes.id))
flag.areas = unique(filter(tmp, prop_lowumi<.5)$problem_areas_genes.id)
cdata$problem_area_flag = ifelse(cdata$problem_areas_genes.id %in% flag.areas, TRUE, FALSE)

cat("\n\nRemove spots (off tissue):\n")
table(cdata[,c("in_tissue","remove_spots")])
tmp = cdata[cdata$in_tissue,]

cat("\n\nRemove in_tissue spots (individual problem area criteria):\n")
table(tmp[,c("true_edges","remove_spots")])
table(tmp[,c("problem_areas_binary","remove_spots")])
table(tmp[,c("lowumi","remove_spots")])
cat("Remove in_tissue spots (summary problem area criteria):\n")
table(tmp[,c("remove_problem.areas","remove_spots")])

cat("\n\nRemove in_tissue spots (individual SpotSweeper criteria):\n")
table(tmp[,c("umi_local.outlier","remove_spots")])
table(tmp[,c("genes_local.outlier","remove_spots")])
table(tmp[,c("chrM.ratio_local.outlier","remove_spots")])
cat("Remove in_tissue spots (summary SpotSweeper criteria):\n")
table(tmp[,c("spotsweeper_outlier","remove_spots")])

cat("\n\nTotal number of spots removed for any reason:\n")
table(cdata$remove_spots)
cat("Remaining spots flagged by problem areas:\n")
table(cdata[cdata$remove_spots==FALSE,"problem_area_flag"])

write.csv(cdata, "processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv", row.names=T)
cat("\nFinal colData saved to: processed-data/03_QC/colData_edges-problem-areas_spotsweeper_FINAL.csv\n")

#load spe
system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
cat("\nDim spe:",dim(spe),"\n")

stopifnot(identical(rownames(cdata),rownames(colData(spe))))
spe$remove_spots = cdata$remove_spots

#modify colData for plotting
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

#boxplots
df = filter(as.data.frame(colData(spe)), in_tissue==TRUE) %>% group_by(sample_id2) %>% mutate(n_spots_removed=sum(remove_spots)) %>%
  ungroup() %>% filter(remove_spots==FALSE) %>% mutate(n_spots_removed_discrete=cut(n_spots_removed, breaks=c(0,50,100,300,600,1200)))

cat("Plotting QC boxplots...",format(Sys.time()),"\n")
p1 <- ggplot(df, aes(x=sample_id2, y=sum_umi, fill=n_spots_removed_discrete))+
  geom_boxplot(color="grey", outlier.size=.5)+geom_hline(aes(yintercept=1000), lty=2, color="red3")+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), limits=c(100,50000),
                     breaks=c(10^(2:5)), labels=c("100","1k","10k","100k"))+
  facet_wrap(vars(box_facets), ncol=1, scales="free_x")+
  scale_fill_viridis_d(option="F")+
  labs(x="", y="sum_umi (log10 scale)", title="Library size - kept spots only", fill="# spots\ndiscarded")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p2 <- ggplot(df, aes(x=sample_id2, y=sum_gene, fill=n_spots_removed_discrete))+
  geom_boxplot(color="grey", outlier.size=.5)+
  facet_wrap(vars(box_facets), ncol=1, scales="free_x")+
  scale_fill_viridis_d(option="F")+
  labs(x="", y="sum_gene", title="Detected genes - kept spots only", fill="# spots\ndiscarded")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())
p3 <- ggplot(df, aes(x=sample_id2, y=expr_chrM_ratio, fill=n_spots_removed_discrete))+
  geom_boxplot(color="grey", outlier.size=.5)+
  facet_wrap(vars(box_facets), ncol=1, scales="free_x")+
  scale_fill_viridis_d(option="F")+
  labs(x="", y="expr_chrM_ratio", title="Mitochondrial fraction - kept spots only", fill="# spots\ndiscarded")+theme_bw()+
  theme_bw()+theme(#legend.position="bottom", 
    axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
    strip.placement = "inside", strip.text=element_blank(),
    strip.background = element_blank())

pdf(file="plots/03_QC/filtered_qc-metrics_boxplot.pdf", width=9, height=12)
	p1
	p2
	p3
dev.off()
cat("Saved to: plots/03_QC/filtered_qc-metrics_boxplot.pdf\n")

#by slide spot plots
spe = spe[,spe$remove_spots==FALSE]
cat("\nDim spe remaining:",dim(spe),"\n")

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
	suppressMessages(plotSpots(x, annotate="sum_umi", point_size=0.3, sample_id="sample_id2")+
		scale_color_gradient(low="white", high="navy", labels=function(y) paste0(y/1000,"k"))+
		facet_wrap(vars(sample_id2), ncol=4)+
		theme(plot.title=element_blank(),
			strip.background = element_rect(fill="transparent", color="transparent"),
			panel.background=element_rect(fill="grey30"))
	)
})

cat("Compiling library size plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/filtered_library-size_spot-plots.pdf", width=12, height=16)
	libList[[1]]
	libList[[2]]
	libList[[3]]
	libList[[4]]
	libList[[5]]
dev.off()
cat("Saved to: plots/03_QC/filtered_library-size_spot-plots.pdf\n")

cat("\nPlotting # genes detected...",format(Sys.time()),"\n")
geneList = lapply(slideList, function(x) {
	suppressMessages(plotSpots(x, annotate="sum_gene", point_size=0.3, sample_id="sample_id2")+
		scale_color_gradient(low="white", high="navy", labels=function(y) paste0(y/1000,"k"))+
		facet_wrap(vars(sample_id2), ncol=4)+
		theme(plot.title=element_blank(),
                        strip.background = element_rect(fill="transparent", color="transparent"),
			panel.background=element_rect(fill="grey30"))
	)
})

cat("Compiling n genes plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/filtered_n-genes_spot-plots.pdf", width=12, height=16)
	geneList[[1]]
	geneList[[2]]
	geneList[[3]]
	geneList[[4]]
	geneList[[5]]
dev.off()
cat("Saved to: plots/03_QC/filtered_n-genes_spot-plots.pdf\n")

cat("\nPlotting chrM ratio",format(Sys.time()),"\n")
mitoList = lapply(slideList, function(x) {
	suppressMessages(plotSpots(x, annotate="expr_chrM_ratio", point_size=0.3, sample_id="sample_id2")+
		scale_color_gradient("mito\nfraction", low="white", high="navy")+
		facet_wrap(vars(sample_id2), ncol=4)+
		theme(plot.title=element_blank(), 
                        strip.background = element_rect(fill="transparent", color="transparent"),
			panel.background=element_rect(fill="grey30"))
	)
})

cat("Compiling chrM ratio plots...",format(Sys.time()),"\n")
pdf(file="plots/03_QC/filtered_chrM-ratio_spot-plots.pdf", width=12, height=16)
	mitoList[[1]]
	mitoList[[2]]
	mitoList[[3]]
	mitoList[[4]]
	mitoList[[5]]
dev.off()
cat("Saved to: plots/03_QC/filtered_chrM-ratio_spot-plots.pdf\n")


## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
