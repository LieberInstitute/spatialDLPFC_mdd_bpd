library(SpatialExperiment)
library(ggspavis)
library(here)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

s1 = unique(spe$slide)
s1 = lapply(s1, function(x) spe[,colData(spe)$slide==x])

metric.list = list("3MAD-slide"=c("sum_3MAD.outlier_slide","genes_3MAD.outlier_slide","chrM.ratio_3MAD.outlier_slide"),
"3MAD-sample"=c("sum_3MAD.outlier_sample","genes_3MAD.outlier_sample","chrM.ratio_3MAD.outlier_sample"),
"local"=c("sum_local.outlier","genes_local.outlier","chrM.ratio_local.outlier"))

slide_function <- function(spe_object, mlist) {
	l1 = unique(spe_object$sample_id)
	names(l1) = lapply(l1, function(x) unique(colData(spe_object)[spe_object$sample_id==x,"brain"]))
	l1 = lapply(l1, function(x) spe_object[,colData(spe_object)$sample_id==x])

	glist = lapply(seq_along(l1), function(x) {
		g = grid::rasterGrob(imgRaster(l1[[x]]))
		ggplot()+annotation_custom(g)+ggtitle(names(l1)[x])
	})

	gridExtra::grid.arrange(glist[[1]],
		plotSpots(l1[[1]], annotate=mlist[1])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[1]),
		plotSpots(l1[[1]], annotate=mlist[2])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[2]),
		plotSpots(l1[[1]], annotate=mlist[3])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[3]),
		glist[[2]],
		plotSpots(l1[[2]], annotate=mlist[1])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[1]),
		plotSpots(l1[[2]], annotate=mlist[2])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[2]),
		plotSpots(l1[[2]], annotate=mlist[3])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[3]),
		glist[[3]],
		plotSpots(l1[[3]], annotate=mlist[1])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[1]),
		plotSpots(l1[[3]], annotate=mlist[2])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[2]),
		plotSpots(l1[[3]], annotate=mlist[3])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[3]),
		glist[[4]],
		plotSpots(l1[[4]], annotate=mlist[1])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[1]),
		plotSpots(l1[[4]], annotate=mlist[2])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[2]),
		plotSpots(l1[[4]], annotate=mlist[3])+scale_color_manual(values=c("grey","red"))+theme(legend.position="none")+ggtitle(mlist[3]),
		ncol=4,
		top=grid::textGrob(paste("Slide",unique(spe_object$slide)), gp=grid::gpar(fontsize=24))) 
}
pdf(here("plots", "03_QC", paste0("spatial_outliers_",names(metric.list)[1],".pdf")), width=12, height=12)
for(i in seq_along(s1)) {
	slide_function(s1[[i]], metric.list[[1]])
}
dev.off()

pdf(here("plots", "03_QC", paste0("spatial_outliers_",names(metric.list)[2],".pdf")), width=12, height=12)
for(i in seq_along(s1))	{                       
        slide_function(s1[[i]], metric.list[[2]])
}
dev.off()

pdf(here("plots", "03_QC", paste0("spatial_outliers_",names(metric.list)[3],".pdf")), width=12, height=12)
for(i in seq_along(s1))	{                       
        slide_function(s1[[i]], metric.list[[3]])
}
dev.off()
