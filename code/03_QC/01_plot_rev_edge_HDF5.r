setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(S4Vectors)
	library(ggspavis)
	library(dplyr)
	library(here)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
mdata = read.csv("processed-data/03_QC/spe_n120_edge-detection_colData.csv", row.names=1)
stopifnot(identical(spe$key, mdata$key))
colData(spe) <- as(mdata, "DFrame")

spe$plot_tissue = factor(paste(spe$keep_spots, spe$edge_outlier),
	levels=c("image perimeter FALSE","keep FALSE","off tissue FALSE","tissue perimeter FALSE",
		"image perimeter TRUE","keep TRUE","tissue perimeter TRUE"),
	labels=c("image perimeter","tissue","off tissue","tissue perimeter",
		"image perimeter: edge removed","tissue: edge removed","tissue perimeter: edge removed"))

#make list of lists so that 6 slides per pdf page
dummy = c("V13B23-283","V13B23-339")
remainder = setdiff(unique(spe$slide), dummy)
slideList = list(c(dummy, remainder[1:5]),remainder[6:11],remainder[12:17],remainder[18:23],remainder[24:29])

slideList = lapply(slideList, function(x) 
	unlist(lapply(x, function(y) 
		unique(colData(spe)[spe$slide==y,"sample_id"])
	))
)
color.palette = c("pink","palegoldenrod","powderblue","grey","red","goldenrod1","black")
names(color.palette) = c("image perimeter","tissue","off tissue","tissue perimeter","image perimeter: edge removed","tissue: edge removed","tissue perimeter: edge removed")

plotList = lapply(slideList, function(x) {
	cat("Plotting...",format(Sys.time()),"\n")
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe[,colData(spe)$sample_id==y])
	
	lapply(1:length(l1), function(z) 
		plotSpots(l1[[z]], annotate="plot_tissue", in_tissue=NULL, point_size=0.1,
			pal=color.palette)+
		geom_point(show.legend=TRUE, size=.1)+
		scale_color_manual(values=color.palette, drop=F)+
		labs(title=names(l1)[[z]])
	)
})

cat("Compiling plots...",format(Sys.time()),"\n")
pdf(here("plots", "03_QC", "spot-tissue_annotation.pdf"), width=12, height=16)
	PRECAST::drawFigs(plotList[[1]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
	PRECAST::drawFigs(plotList[[2]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
	PRECAST::drawFigs(plotList[[3]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
	PRECAST::drawFigs(plotList[[4]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
	PRECAST::drawFigs(plotList[[5]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()
#ggsave(filename = here("plots","03_QC","spot-tissue_annotation.png"), plot=f1, width=12, height=16, bg="white")
cat("\nplot destination:",here("plots", "03_QC", "spot-tissue_annotation.pdf"),"\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
proc.time()
options(width = 120)
sessionInfo()
