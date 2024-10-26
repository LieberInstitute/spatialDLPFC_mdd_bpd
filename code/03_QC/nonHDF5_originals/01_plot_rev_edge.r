setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(ggspavis)
	library(dplyr)
	library(here)
})

load(here("processed-data","03_QC","spe_demo.Rdata"))


spe$plot_tissue = factor(paste(spe$keep_spots, spe$edge_outlier),
	levels=c("image perimeter FALSE","keep FALSE","off tissue FALSE","tissue perimeter FALSE",
		"image perimeter TRUE","keep TRUE","tissue perimeter TRUE"),
	labels=c("image perimeter","tissue","off tissue","tissue perimeter",
		"image perimeter: edge removed","tissue: edge removed","tissue perimeter: edge removed"))
l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) spe[,colData(spe)$sample_id==x])

p1 <- lapply(1:length(l1), function(x) 
	plotSpots(l1[[x]], annotate="plot_tissue", in_tissue=NULL, point_size=0.1,
		pal=c("pink","palegoldenrod","powderblue","grey",
			"red","goldenrod1","black"))+
	labs(title=names(l1)[[x]])
)

#pdf(here("plots", "03_QC", "spot-tissue_annotation.pdf"), width=12, height=16)
f1 = PRECAST::drawFigs(p1, layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
#dev.off()
ggsave(filename = here("plots","03_QC","spot-tissue_annotation.png"), plot=f1, width=12, height=16, bg="white")
cat("\nplot destination:",here("plots", "03_QC", "spot-tissue_annotation.png"),"\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
