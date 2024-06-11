setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(ggspavis)
	library(here)
})

load(here("processed-data","03_QC","spe_demo.Rdata"))
spe_full = spe
spe_full$plot_tissue = factor(paste(spe_full$keep_spots, spe_full$edge_outlier),
	levels=c("image perimeter FALSE","keep FALSE","off tissue FALSE","tissue perimeter FALSE",
		"image perimeter TRUE","keep TRUE","tissue perimeter TRUE"),
	labels=c("image perimeter","tissue","off tissue","tissue perimeter",
		"image perimeter: edge removed","tissue: edge removed","tissue perimeter: edge removed"))

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

spe$sample = factor(paste(spe$umi_3MAD.outlier_sample | spe$genes_3MAD.outlier_sample, spe$chrM.ratio_3MAD.outlier_sample),
        levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE NA","TRUE TRUE"),
        labels=c("none","mito","umi/genes","both","both"))
spe$slide = factor(paste(spe$umi_3MAD.outlier_slide | spe$genes_3MAD.outlier_slide, spe$chrM.ratio_3MAD.outlier_slide),
        levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE NA","TRUE TRUE"),
        labels=c("none","mito","umi/genes","both","both"))
spe$local = factor(paste(spe$umi_local.outlier | spe$genes_local.outlier, spe$chrM.ratio_local.outlier),
        levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
        labels=c("none","mito","umi/genes","both"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) colData(spe)$sample_id==x)

for(x in seq_along(l1)) {
	cat("Generating plots for",names(l1)[x],"...\n")
	spe.tmp = spe[,l1[[x]]]
	spe_full.tmp = spe_full[,spe_full$brain==names(l1)[x]]
	# tissue/spot annotation
	p1 <- plotSpots(spe_full.tmp, annotate="plot_tissue", in_tissue=NULL,
		pal=c("red","goldenrod1","powderblue","black","pink","palegoldenrod","grey"))+
		ggtitle(names(l1)[x])+labs(color="tissue/spot annotation")+
		theme(legend.position="bottom", legend.direction="vertical", plot.title=element_text(face="bold"))
	# build list of ggplot objects for top row
	p.list1 = list(plotSpots(spe.tmp, annotate="lg10.umi", pal=c("black","white")),
		plotSpots(spe.tmp, annotate="lg10.genes", pal=c("black","white")),
		plotSpots(spe.tmp, annotate="expr_chrM_ratio", pal=c("white","black")))
	# build list of ggplot objects for bottom row
	mlist = c("sample","slide","local")
	color.key = c("none"="grey","mito"="dodgerblue","umi/genes"="black","both"="red")
	p.list2 = lapply(mlist, function(y) {
			plotSpots(spe.tmp, annotate=y, pal=color.key[levels(colData(spe.tmp)[,y])])+
			ggtitle(paste(y,"outliers"))+labs(color="")+
			theme(plot.title=element_text(hjust=.5),
				legend.box.spacing=unit(2,"pt"), legend.spacing=unit(0,"pt"), legend.margin=margin(0,0,0,0),legend.key.width=unit(1,"mm"))
	})
	# combine plots
	lay= rbind(c(1,2,3,4),c(1,5,6,7))
	finalPlot = gridExtra::grid.arrange(p1, p.list1[[1]], p.list1[[2]], p.list1[[3]], p.list2[[1]], p.list2[[2]], p.list2[[3]],
		layout_matrix=lay)
	ggsave(filename=here("plots", "03_QC", paste0(names(l1)[x],"_outlier-overview.png")), finalPlot, width=12, height=6)
	#dev.off()
}

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
