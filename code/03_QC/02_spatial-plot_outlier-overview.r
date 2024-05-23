library(SpatialExperiment)
library(ggspavis)
library(gridExtra)
library(parallel)
library(here)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

spe$sum.outlier_f = factor(paste(spe$sum_3MAD.outlier_sample, spe$sum_local.outlier),
	levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("none","local only","3MAD only","both"))
spe$genes.outlier_f = factor(paste(spe$genes_3MAD.outlier_sample, spe$genes_local.outlier),
	levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("none","local only","3MAD only","both"))
spe$chrM.outlier_f = factor(paste(spe$chrM.ratio_3MAD.outlier_sample, spe$chrM.ratio_local.outlier),
	levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("none","local only","3MAD only","both"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) spe[,colData(spe)$sample_id==x])

for(x in seq_along(l1)) {
	# convert tissue image to grob
	g = grid::rasterGrob(imgRaster(l1[[x]]))
	# calculate # outlier totals
	m1 = t(rbind(table(l1[[x]]$sum.outlier_f),
		table(l1[[x]]$genes.outlier_f),
		table(l1[[x]]$chrM.outlier_f)))
	colnames(m1) = c("umi","genes","chrM")
	m1 = as.data.frame(m1[c("both","3MAD only","local only","none"),])
	### create theme for outlier table
	color.theme <- c(rep(c("#e41a1c","black","#377eb8","lightgrey"), times = c(12)))
	tt <- ttheme_minimal(core=list(bg_params=list(fill= color.theme), fg_params=list(col="white")))
	### make table into grob
	gt = tableGrob(m1, theme=tt)
	# outlier metrics to plot
	mlist2 = c("sum.outlier_f","genes.outlier_f","chrM.outlier_f")
	# build list of ggplot objects for top row
	p.list1 = list(ggplot()+theme_minimal()+annotation_custom(g)+ggtitle(names(l1)[x])+theme(title=element_text(face="bold")),
		plotSpots(l1[[x]], annotate="lg10.sum", pal=c("black","white")),
		plotSpots(l1[[x]], annotate="lg10.genes", pal=c("black","white")),
		plotSpots(l1[[x]], annotate="expr_chrM_ratio", pal=c("white","black")))
	# build list of ggplot objects for bottom row
	p.list2 = c(list(ggplot()+theme_minimal()+annotation_custom(gt)+ggtitle("# outliers")),
		lapply(mlist2, function(y) {
			plotSpots(l1[[x]], annotate=y, pal=c("lightgrey","#377eb8","black","#e41a1c"))+ggtitle(y)+theme(legend.position="none")
		}))
	# for top and bottom row, use precast::draw figs to combine list of ggplot objects into a single gg/ggplot/ggarrange object,
	### then combine into single object with cowplot so that it will plot on a single pdf page
	#pdf(here("plots", "03_QC", paste0(names(l1)[x],"_outlier-overview.pdf")), width=12, height=6)
	finalPlot = cowplot::plot_grid(PRECAST::drawFigs(p.list1, layout.dim=c(1,4), common.legend=FALSE),
		PRECAST::drawFigs(p.list2, layout.dim=c(1,4), common.legend=TRUE, legend.position="none"), nrow=2)
	ggsave(filename=here("plots", "03_QC", paste0(names(l1)[x],"_outlier-overview.pdf")), finalPlot, width=12, height=6)
	#dev.off()
}
