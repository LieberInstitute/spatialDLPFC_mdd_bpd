library(SpatialExperiment)
library(ggspavis)
library(here)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

spe$exclude_outliers = spe$sum_local.outlier | spe$genes_local.outlier | spe$chrM.ratio_local.outlier | spe$exclude_east_edge
spe$extra_3MAD_outliers = spe$sum_3MAD.outlier_sample | spe$genes_3MAD.outlier_sample
spe$plot_outliers = factor(paste(spe$exclude_outliers, spe$extra_3MAD_outliers),
	levels=c("FALSE FALSE","FALSE TRUE","TRUE FALSE","TRUE TRUE"),
	labels=c("keep","3MAD flag","exclude","exclude"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) spe[,colData(spe)$sample_id==x])

lapply(seq_along(l1), function(x) {

filt.spe = l1[[x]][,colData(l1[[x]])$plot_outliers!="exclude"]
named.colors = c("keep"="lightgrey","3MAD flag"="black",exclude="red")
named.colors[unique(l1[[x]]$plot_outliers)]

p.list = list(plotSpots(l1[[x]], annotate = "plot_outliers", sample_id="sample_id",pal=named.colors[unique(l1[[x]]$plot_outliers)], point_size=.2)+labs(color="")+ggtitle(names(l1)[x])+
		theme(title=element_text(face="bold", size=16), legend.position="bottom", legend.text=element_text(size=10), 
		legend.box.spacing=unit(1,"pt"), legend.spacing=unit(0,"pt"), legend.margin=margin(0,0,0,0),legend.key.width=unit(1,"mm")), 
	plotSpots(filt.spe, annotate="lg10.sum", pal=c("black","white"))+ggtitle("post-filt lg10(UMI counts)"),
	plotSpots(filt.spe, annotate="lg10.genes", pal=c("black","white"))+ggtitle("post-filt lg10(n genes)"),
	plotSpots(filt.spe, annotate="expr_chrM_ratio", pal=c("white","black"))+ggtitle("post-filt chrM ratio"))

ggsave(do.call(gridExtra::grid.arrange, c(p.list, ncol=4, bottom="Exclude (red) spots are determined based on local outliers and eastern edge removal.")), filename=here("plots","03_QC",paste0(names(l1)[x],"_filter-preview.pdf")), width=12, height=3)
})
