setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(DeconvoBuddies)
})
set.seed(123)

#PRECAST
cat("\nPRECAST (smoothed)...\n")
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo_sample-smoothed-n1663-k9_norm-filt.Rdata")
marker_stats_MeanRatio <- get_mean_ratio(
    sce = spe_pseudo, 
    assay_name = "logcounts", 
    cellType_col = "smoothed_k9_1663", 
    gene_ensembl = "gene_id", 
    gene_name = "gene_name"
)
write.csv(marker_stats_MeanRatio, "processed-data/06_pseudobulk/PRECAST_smoothed/smoothed-k9-1663_mean-ratio.csv", row.names=F)

#Seurat label transfer
cat("\nSeurat pc20...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc20_norm-filt.Rdata")
marker_stats_MeanRatio <- get_mean_ratio(
    sce = spe_pseudo, 
    assay_name = "logcounts", 
    cellType_col = "seurat_label", 
    gene_ensembl = "gene_id", 
    gene_name = "gene_name"
)
write.csv(marker_stats_MeanRatio, "processed-data/06_pseudobulk/Seurat/seurat-pc20_mean-ratio.csv", row.names=F)

cat("\nSeurat pc30...\n")
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo_sample-seurat-pc30_norm-filt.Rdata")
marker_stats_MeanRatio <- get_mean_ratio(
    sce = spe_pseudo,
    assay_name = "logcounts",
    cellType_col = "seurat_label",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)
write.csv(marker_stats_MeanRatio, "processed-data/06_pseudobulk/Seurat/seurat-pc30_mean-ratio.csv", row.names=F)

#SZBD
cat("\nSZBD snRNAseq...\n")
load("processed-data/06_pseudobulk/SZBDMulti-seq/sce_control_pseudo_indivID-low-res_norm-filt.Rdata")
marker_stats_MeanRatio <- get_mean_ratio(
    sce = sce_pseudo, 
    assay_name = "logcounts", 
    cellType_col = "seurat_low.res", 
    gene_ensembl = "gene_id", 
    gene_name = "gene_name"
)
write.csv(marker_stats_MeanRatio, "processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_mean-ratio.csv", row.names=F)

sce_pseudo$seurat_nrn = factor(sce_pseudo$seurat_low.res, levels=c("Micro.Vasc","Astro","L2","L3","L4","L5","L6","Inhb","Oligo"),
	labels=c("Micro.Vasc","Astro","Nrn","Nrn","Nrn","Nrn","Nrn","Nrn", "Oligo"))
marker_stats_MeanRatio <- get_mean_ratio(
    sce = sce_pseudo,
    assay_name = "logcounts",
    cellType_col = "seurat_nrn",
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)
write.csv(marker_stats_MeanRatio, "processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_nrns-merge_mean-ratio.csv", row.names=F)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
