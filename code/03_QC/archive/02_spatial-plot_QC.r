library(SpatialExperiment)
library(ggplot2)
library(ggspavis)
library(parallel)
library(here)

load(here("processed-data","03_QC","spe_demo-filt.Rdata"))

l1 = unique(spe$sample_id)
names(l1) = lapply(l1, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
l1 = lapply(l1, function(x) colData(spe)$sample_id==x)

mlist = c("lg10.sum","lg10.genes","expr_chrM_ratio")
names(mlist) = c("log10(UMI counts)","log10(n genes)","chrM ratio") 

p.list = mclapply(mlist, function(y) {
	lapply(seq_along(l1), function(x) {
		sub.spe = spe[,l1[[x]]]
		plotSpots(spe[,l1[[x]]], annotate = y, pal=c("white", "black"))+ggtitle(names(l1)[[x]])
})
}, mc.cores=3)

pdf(here("plots", "03_QC", "spatial-expr_QC-metrics.pdf"), width=12, height=16)
do.call(gridExtra::grid.arrange, c(p.list[[1]], ncol=4, top=names(p.list)[1]))
do.call(gridExtra::grid.arrange, c(p.list[[2]], ncol=4, top=names(p.list)[2]))
do.call(gridExtra::grid.arrange, c(p.list[[3]], ncol=4, top=names(p.list)[3]))
dev.off()
