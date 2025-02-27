setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(DelayedArray)
	library(scran)
	library(dplyr)
	library(BiocParallel)
})
set.seed(123)

spe <- loadHDF5SummarizedExperiment(dir="processed-data/04_feature_selection/", prefix="spe_n120_postQC_norm_")

fileList = list.files("processed-data/04_feature_selection/per-sample_svgs")
length(fileList) 
resList = lapply(fileList, function(x) {
  name1 = substr(x, start=0, stop=13)
  df = read.csv(paste0("processed-data/04_feature_selection/per-sample_svgs/",x), row.names=1)
  rownames(df) <- NULL
  df$sample_id = name1
  return(df)
})
results.df = do.call(rbind, resList)

exclude.genes = readRDS("processed-data/04_feature_selection/batch-effect-genes_dummyslide-sample-sex-condition_list.rds")
cat("\nNumber of excluded genes:\n")
length(unique(unlist(exclude.genes))) #40

results.df$remove_gene = results.df$gene_id %in% unlist(exclude.genes)

#regular top.genes filter
top.genes = filter(results.df, padj<.05 & rank<=500, n_spots_nonzero>500) 
cat("\nNumber of top genes (padj<.05, rank<=500, >500 spots):\n")
length(unique(top.genes$gene_name))
#additional filter to decile 10
avg.expr = read.csv("processed-data/04_feature_selection/nnSVG-filtered-genes_avg-logcounts.csv", row.names=1) %>%
	tibble::rownames_to_column(var="gene_id")
top.decile = left_join(top.genes, avg.expr) %>% filter(decile==10) %>% pull(gene_id) %>% unique() 
cat("\nNumber of top genes in top decile:\n")
length(top.decile)

#subset and run correlatePairs
spe_sub = spe[top.decile,]
cat("\nSubset spe...\n")
dim(spe_sub)
cat("\nCorrelate pairs...\n")
format(Sys.time(), tz="EST")
corr_top = correlatePairs(spe_sub, block=spe_sub$sample_id, equiweight=T, assay.type="logcounts",
	BPPARAM=MulticoreParam(workers=12))
format(Sys.time(), tz="EST")

tmp1 = avg.expr[,c("gene_id","gene_name","avg_expr")]
tmp2 = avg.expr[,c("gene_id","gene_name","avg_expr")]
colnames(tmp1) = c("gene1","gene1_name","gene1_avg_expr")
colnames(tmp2) = c("gene2","gene2_name","gene2_avg_expr")
corr_top = left_join(as.data.frame(corr_top), tmp1) %>% left_join(tmp2)
write.csv(corr_top, "processed-data/04_feature_selection/top-genes-padj05-rank500_decile10_correlation.csv", row.names=F)
cat("\nCorrelation results saved to: processed-data/04_feature_selection/top-genes-padj05-rank500_decile10_correlation.csv\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

