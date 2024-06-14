setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(Seurat)
	library(PRECAST)
	library(ggspavis)
	library(cowplot)
	library(here)
})

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
sub.spe = spe[,spe$sample_id %in% unique(spe$sample_id)[c(1,5,9,13,17,21)]]
test.samples = c("bindev.2k","bindev.3k","svg")

for(i in 1:3) {
	cat("\nLoading",test.samples[i],"precast results...\n")
	load(here("processed-data","05_clustering",paste0("srt_precast_6samp_",test.samples[i],".Rdata")))
	sub.spe$seurat_key = paste(sub.spe$sample_id, colnames(sub.spe), sep="_")
	sub.spe$seurat_key_match = sub.spe$seurat_key %in% colnames(seuInt)
	table(sub.spe$seurat_key_match)
	if(length(table(sub.spe$seurat_key_match))>1) {
		cat("Warning! Unequal number of spots between subsampled spe and seuInt. Check PRECAST log to see if neighbors found for all spots.\n")
		cat("Filtering subsampled spe to match seuInt size...\n")
		sub.tmp = sub.spe[,sub.spe$seurat_key_match]
		stopifnot(identical(sub.tmp$seurat_key, colnames(seuInt)))
		sub.spe$precast_clusters = NA
		colData(sub.spe)[sub.spe$seurat_key_match,"precast_clusters"] = seuInt$cluster
		sub.spe$precast_clusters = as.factor(sub.spe$precast_clusters)
	} else {
		stopifnot(identical(sub.spe$seurat_key, colnames(seuInt)))
		sub.spe$precast_clusters = as.factor(seuInt$cluster)
}
	colnames(colData(sub.spe))[ncol(colData(sub.spe))] = paste0("precast_clusters_",test.samples[i])
}

cols_cluster <- chooseColors(palettes_name = "Classic 20", n_colors = 7, plot_colors = TRUE)

l2 = unique(sub.spe$sample_id)
names(l2) = lapply(l2, function(x) unique(colData(sub.spe)[sub.spe$sample_id==x,"brain"]))
l2 = lapply(l2, function(x) sub.spe[,colData(sub.spe)$sample_id==x])

p1 <- lapply(seq_along(l2), function(x) {
	plotSpots(l2[[x]], annotate="precast_clusters_bindev.2k", point_size=.1)+
	scale_color_manual(values=cols_cluster)+
	labs(title=names(l2)[x],color="clus")
})
f1 <- plot_grid(plotlist=p1, ncol=1)
p2 <- lapply(seq_along(l2), function(x) {
        plotSpots(l2[[x]], annotate="precast_clusters_bindev.3k", point_size=.1)+
        scale_color_manual(values=cols_cluster)+
        labs(title=names(l2)[x],color="clus")
})
f2 <- plot_grid(plotlist=p2, ncol=1)
p3 <- lapply(seq_along(l2), function(x) {
        plotSpots(l2[[x]], annotate="precast_clusters_svg", point_size=.1)+
        scale_color_manual(values=cols_cluster)+
        labs(title=names(l2)[x],color="clus")
})
f3 <- plot_grid(plotlist=p3, ncol=1)

ggsave(file=here("plots","05_clustering","precast_test-feature-list_6samp.png"), 
	plot=gridExtra::grid.arrange(f1, f2, f3, ncol=3, top="left: bindev.2k; middle: bindev.3k; right: svg"), width=7, height=12, bg="white")
cat("\nplot destination:",here("plots","05_clustering","precast_test-feature-list_6samp.png"),"\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
