################
### currently I am getting out of memory error/ failure to complete for BiasDetect featureSelect() which generates the binomial deviance results
###### code/04_feature_selection/test_BiasDetect_HDF5.r/.sh
### present solution: generate dataframe using code/04_feature_selection/test_bindev_HDF5.r/.sh
### those results will be loaded here
### in the future this code will contain both the binomial deviance calculations/ dataframe generation and the SD plots and results

setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(HDF5Array)
        library(DelayedArray)
        library(BiasDetect)
	library(dplyr)
	library(ggplot2)
})
set.seed(123)

df = read.csv("processed-data/04_feature_selection/test_bindev_batch-slide.csv")
batch_df = rename(df, gene="gene_id")

batch_df$d_diff <- (batch_df$dev_default-batch_df$dev_batch)/batch_df$dev_batch
mean_dev <- mean(batch_df$d_diff)
sd_dev <- sd(batch_df$d_diff)
batch_df$nSD_dev <- (batch_df$d_diff - mean_dev) / sd_dev

batch_df$r_diff <- batch_df$rank_batch-batch_df$rank_default
mean_rank <- mean(batch_df$r_diff)
sd_rank <- sd(batch_df$r_diff)
batch_df$nSD_rank <- (batch_df$r_diff - mean_rank) / sd_rank

outList = list()
outList$plots = biasDetect(batch_df, nSD_dev=5, visual=T)
outList$genes = biasDetect(batch_df, nSD_dev=5, visual=F)
cat("\nBiased genes (n=",paste0(length(outList$genes),")"),"found with nSD=5:\n")
outList$genes
saveRDS(outList, "processed-data/04_feature_selection/test_BiasDetect-biasDetect.rda")
cat("\n\nSaved plots and biased gene list to: processed-data/04_feature_selection/test_biasDetect_batch-slide.rda\n")

ggsave("plots/04_feature_selection/BiasDetect_5sd.png",
	plot= gridExtra::grid.arrange(outList$plots[[1]], outList$plots[[2]], ncol=2), bg="white",
	height=5, width=10, units="in")
cat("\nSaved plots (visualized) to: plots/04_feature_selection/BiasDetect_5sd.png\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
