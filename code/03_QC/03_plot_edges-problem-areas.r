setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(ggspavis)
})

system.time(spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_"))
dim(spe)

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas.csv", row.names=1)
cdata$problem_areas_genes_binary = !is.na(cdata$problem_areas_genes.id)

stopifnot(identical(rownames(cdata),rownames(colData(spe))))


#plots of problem areas <10, <20, and bigger
colData(spe)$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size<20, "less 20", "big")
colData(spe)$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size<10, "less 10", colData(spe)$problem_areas_grouped)
colData(spe)$problem_areas_grouped = ifelse(cdata$problem_areas_genes.size==0, "none", colData(spe)$problem_areas_grouped)
colData(spe)[cdata$edge_outlier_genes,"problem_areas_grouped"] = "EDGE"

#plotting
spe2 <- spe[,spe$in_tissue]
dim(spe2)
spe2$problem_areas_grouped = as.factor(spe2$problem_areas_grouped)
spe2$dummy_slide = spe2$slide
colData(spe2)[spe2$slide=="V13B23-283","dummy_slide"] = "joint-283-339"
colData(spe2)[spe2$slide=="V13B23-339","dummy_slide"] = "joint-283-339"

colData(spe2)$facet_violin = paste(spe2$round, spe2$dummy_slide)
colData(spe2)$facet_violin = ifelse(spe2$dummy_slide=="joint-283-339", "joint-283-339", colData(spe2)$facet_violin)
seed = levels(as.factor(colData(spe2)$facet_violin))

colData(spe2)$facet_spots = paste(spe2$round, spe2$sample_id)
colData(spe2)$facet_spots = ifelse(spe2$brnum=="Br5366", paste(spe2$facet_spots, spe2$brnum), spe2$facet_spots)

slideList2 = list(c(seed[2:6],seed[1]),seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList2 = lapply(slideList2, function(x) {
	unlist(lapply(x, function(y)
		sort(unique(colData(spe2)[spe2$facet_violin==y,"facet_spots"]))
  ))
})

color.palette = c("navy","darkgreen","black","red","grey80")
names(color.palette) <- c("EDGE","big","less 10","less 20","none")

cat("\nGenerating spot plots...\n")
plotList = lapply(slideList2, function(x) {
	l1 = x; names(l1) = x
	l1 = lapply(l1, function(y) spe2[,colData(spe2)$facet_spots==y])

	lapply(1:length(l1), function(z) 
		suppressMessages(plotSpots(l1[[z]], annotate="problem_areas_grouped", #in_tissue=NULL, 
			point_size=0.2,
			pal=color.palette)+
			geom_point(show.legend=TRUE, size=.1)+
			scale_color_manual("",values=color.palette, drop=F)+
			labs(title=names(l1)[[z]])
		)
	)
})


pdf(file="plots/03_QC/edges_problem-areas_size-grouped_spot-plots.pdf", width=12, height=16)
PRECAST::drawFigs(plotList[[1]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[2]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[3]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[4]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
PRECAST::drawFigs(plotList[[5]], layout.dim = c(6, 4), common.legend = TRUE, legend.position = "right", align = "hv")
dev.off()
cat("\nPlot saved to: plots/03_QC/edges_problem-areas_size-grouped_spot-plots.pdf\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
