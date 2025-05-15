setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(ggspavis)
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(gridExtra)
        library(scater)
})

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n119_postQC_norm_")
cat("Dim spe:",dim(spe),"\n")

spe$dummy_slide = ifelse(spe$slide %in% c("V13B23-339","V13B23-283"), "joint-283-339", spe$slide)
spe$array2 = ifelse(spe$slide=="V13B23-283", "A1", spe$array)

seed = levels(as.factor(spe$dummy_slide))
slideList = list(seed[1:6],seed[7:12], seed[13:18], seed[19:24], seed[25:30])
slideList = lapply(slideList, function(x) {
	do.call(cbind, lapply(x, function(y) spe[,spe$dummy_slide==y]))
})

marker1 = "TMSB10"
markerList1 = lapply(slideList, function(x) {
	suppressMessages(
		plotSpots(x, annotate=marker1, point_size=.5, feature_names="gene_name", assay_name="logcounts", sample_id="sample_id")+
		scale_color_viridis_c(marker1, option="F", direction=-1, limits=c(0,max(logcounts(spe)[rowData(spe)$gene_name==marker1,])))+
		facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
		theme(panel.background=element_rect(fill="grey30"))
	)
})

pdf(file=paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"), width=12, height=16)
	markerList1[[1]]
	markerList1[[2]]
	markerList1[[3]]
	markerList1[[4]]
	markerList1[[5]]
	#markerList1[[6]]
dev.off()
cat("\nSaved",marker1,"plots to:", paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"),"\n")

marker1	= "AL627171.2"
markerList1 = lapply(slideList,	function(x) {
        suppressMessages( 
                plotSpots(x, annotate=marker1, point_size=.5, feature_names="gene_name", assay_name="logcounts", sample_id="sample_id")+
                scale_color_viridis_c(marker1, option="F", direction=-1, limits=c(0,max(logcounts(spe)[rowData(spe)$gene_name==marker1,])))+
                facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
                theme(panel.background=element_rect(fill="grey30"))
        )
})

pdf(file=paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"), width=12, height=16)
        markerList1[[1]]
        markerList1[[2]]
        markerList1[[3]]
        markerList1[[4]]
        markerList1[[5]]
        #markerList1[[6]]
dev.off()
cat("\nSaved",marker1,"plots to:", paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"),"\n")

marker1	= "MTRNR2L12"
markerList1 = lapply(slideList,	function(x) {
        suppressMessages( 
                plotSpots(x, annotate=marker1, point_size=.5, feature_names="gene_name", assay_name="logcounts", sample_id="sample_id")+
                scale_color_viridis_c(marker1, option="F", direction=-1, limits=c(0,max(logcounts(spe)[rowData(spe)$gene_name==marker1,])))+
                facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
                theme(panel.background=element_rect(fill="grey30"))
        )
})

pdf(file=paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"), width=12, height=16)
        markerList1[[1]]
        markerList1[[2]]
        markerList1[[3]]
        markerList1[[4]]
        markerList1[[5]]
        #markerList1[[6]]
dev.off()
cat("\nSaved",marker1,"plots to:", paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"),"\n")

marker1	= "CST3"
markerList1 = lapply(slideList,	function(x) {
        suppressMessages( 
                plotSpots(x, annotate=marker1, point_size=.5, feature_names="gene_name", assay_name="logcounts", sample_id="sample_id")+
                scale_color_viridis_c(marker1, option="F", direction=-1, limits=c(0,max(logcounts(spe)[rowData(spe)$gene_name==marker1,])))+
                facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
                theme(panel.background=element_rect(fill="grey30"))
        )
})

pdf(file=paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"), width=12, height=16)
        markerList1[[1]]
        markerList1[[2]]
        markerList1[[3]]
        markerList1[[4]]
        markerList1[[5]]
        #markerList1[[6]]
dev.off()
cat("\nSaved",marker1,"plots to:", paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"),"\n")

marker1	= "MT3"
markerList1 = lapply(slideList,	function(x) {
        suppressMessages( 
                plotSpots(x, annotate=marker1, point_size=.5, feature_names="gene_name", assay_name="logcounts", sample_id="sample_id")+
                scale_color_viridis_c(marker1, option="F", direction=-1, limits=c(0,max(logcounts(spe)[rowData(spe)$gene_name==marker1,])))+
                facet_grid(rows=vars(dummy_slide), cols=vars(array2))+
                theme(panel.background=element_rect(fill="grey30"))
        )
})

pdf(file=paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"), width=12, height=16)
        markerList1[[1]]
        markerList1[[2]]
        markerList1[[3]]
        markerList1[[4]]
        markerList1[[5]]
        #markerList1[[6]]
dev.off()
cat("\nSaved",marker1,"plots to:", paste0("plots/04_feature_selection/",marker1,"-logcounts_spot-plots.pdf"),"\n")

#not run
##weird slide effect
#raw_spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe", prefix="spe_n120_")
#spe_example = raw_spe[,raw_spe$slide=="V13B23-301"] #AVP and OXT slide
#avp.spots = colnames(spe_example)[counts(spe_example)[rowData(spe_example)$gene_name=="AVP",]>3]
#spe_example$avp.spots = colnames(spe_example) %in% avp.spots
#oxt.spots = colnames(spe_example)[counts(spe_example)[rowData(spe_example)$gene_name=="OXT",]>3]
#spe_example$oxt.spots = colnames(spe_example) %in% oxt.spots

#plotVisium(spe_example[,spe_example$array=="D1"], spots=T, image=T, annotate=rownames(spe_example)[rowData(spe_example)$gene_name=="AVP"], 
#           highlight="avp.spots", 
#           facets="sample_id", assay = "counts")+
#  facet_wrap(vars(sample_id), ncol=4)+
#  scale_fill_gradient2(low=alpha("white",0), mid="white",high="black", midpoint=10)+
#  scale_color_manual(values=c(alpha("white",0), "navy"))

#plotVisium(spe_example[,spe_example$array=="A1"], spots=T, image=T, 
#           annotate=rownames(spe_example)[rowData(spe_example)$gene_name=="OXT"], 
#           highlight="oxt.spots", 
#           facets="sample_id", assay = "counts")+
#  facet_wrap(vars(sample_id), ncol=4)+
#  scale_fill_gradient2(low=alpha("white",0), mid="white",high="black", midpoint=10)+
#  scale_color_manual(values=c(alpha("white",0), "navy"))

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
