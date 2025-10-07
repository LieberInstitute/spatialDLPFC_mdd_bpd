setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(ggspavis)
	library(dplyr)
	library(ggrastr)
	library(gridExtra)
})

set.seed(123)

cpList = readRDS("plots/colorPalettes.rds")
#fill.palette = c(cpList$smoothed.light, "low UMI"="grey50", "GABA"="white", "Vasc"="white")
fill.palette = c(cpList$smoothed.light, "low UMI"="grey50")

cdata = read.csv("processed-data/05_clustering/PRECAST/colData_conservative_all-precast-clusters.csv", row.names=1)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC-conservative_norm_")
stopifnot(identical(rownames(cdata), colnames(spe)))
#spe$smoothed_k9_1663 = factor(cdata$smoothed_k9_1663_f, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI","GABA","Vasc"),
#	labels=c("L1","L2","L3.4","L5","L6","WM","low UMI","GABA","Vasc"))
spe$smoothed_k7_1626 = factor(cdata$smoothed_k7_1626, levels=c("L1","L2","L3/4","L5","L6","WM","low UMI"),
                                labels=c("L1","L2","L3.4","L5","L6","WM","low UMI"))

#create spe_sub for each page of plots to modify coordinates for prettier plots
#to not have 31 slides, move V13B23-283 with its original group
spe$slide2 = ifelse(spe$slide=="V13B23-283","V13B23-339",spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "A1", spe$array)
slides = sort(unique(spe$slide2))
slideList = list(slides[1:6], slides[7:12], slides[13:18], slides[19:24], slides[25:30])
str(slideList)

plotList = lapply(slideList, function(x) {
	spe_sub = spe[,spe$slide2 %in% x]
	spe_sub$facet_row = factor(spe_sub$slide2, levels=x)
	spe_sub$facet_col = factor(spe_sub$array2, levels= c("A1","B1","C1","D1"))

	mod_spatialCoords2 = spatialCoords(spe_sub)
	for (i in unique(spe_sub$sample_id)) {
		tmp = mod_spatialCoords2[colData(spe_sub)$sample_id==i,]
		mod_spatialCoords2[colData(spe_sub)$sample_id==i,1] = tmp[,1]-min(tmp[,1])
		mod_spatialCoords2[colData(spe_sub)$sample_id==i,2] = tmp[,2]-min(tmp[,2])
	}
	sample.name.df = distinct(as.data.frame(colData(spe_sub)[,c("facet_row","facet_col","sample_id")])) %>%
		mutate(x1=0, y1=-60)
	suppressMessages({
		p1 <- plotSpots(spe_sub, x_coord=mod_spatialCoords2[,1], y_coord=mod_spatialCoords2[,2],
			sample_id="sample_id", annotate="smoothed_k9_1663", point_size=.1)+
		scale_color_manual("PRECAST\n(smoothed)", values=fill.palette)+
		facet_grid(rows=vars(facet_row), cols=vars(facet_col))+
		#label each panel with small text of sample_id
		geom_text(data=sample.name.df, aes(x=x1, y=y1, label=sample_id),
		#geom_text(aes(x=0, y= -60, label= ifelse(.data[["array_row"]]==5 & .data[["array_col"]]==13, .data[["sample_id"]], "")), 
			hjust=0, vjust=0, color="black", size=2)+
		theme(text=element_text(size=10), strip.background = element_rect(fill=NA, color=NA), panel.border=element_rect(fill=NA, color=NA),
			strip.text.x = element_blank(), strip.text.y = element_blank(),
			legend.title=element_text(margin=margin(0,0,4,0,"pt")), legend.box.spacing = unit(2,"pt"),
			legend.margin=margin(0,0,0,0,"pt"), legend.box.margin = margin(0,2,0,2,"pt"))
	})
	return(rasterize(p1, dpi=150))
})

ggsave(file="plots/05_clustering/PRECAST/PRECAST_conservative_n1626-k7_smoothed_low-UMI-cluster.pdf", marrangeGrob(plotList, ncol=1, nrow=1, top=NULL), height=10, width=8)
cat("\nPlots saved to: plots/05_clustering/PRECAST/PRECAST_conservative_n1626-k7_smoothed_low-UMI-cluster.pdf\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
