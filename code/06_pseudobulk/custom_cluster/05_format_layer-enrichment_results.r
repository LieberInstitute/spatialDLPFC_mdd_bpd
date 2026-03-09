setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(HDF5Array)
	library(edgeR)
	library(dplyr)
	library(ggplot2)
	library(scater)
})
set.seed(123)

# load in unfiltered spe object which should have whole gene universe
### the information here (https://www.synapse.org/Synapse:syn22963646) indicates that they used
### the hg38 genome for their reference which is what we used too (?)
spe <- loadHDF5SummarizedExperiment(dir="processed-data/02_build_spe/", prefix="spe_n120_")
rdata = as.data.frame(rowData(spe))
dim(rdata) #36601 7
rm(spe)

##load in sce for heatmap and dotplots
#load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
#cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
#precast_levels= c("L1","L2","L3.4","L5","L6","WM")
#spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$smoothed_k9_1663),
#                            levels=as.character(outer(cond_sex, precast_levels, paste)))
#colnames(spe_summ) <- spe_summ$sample_id

#load in enrichment results and format them
resList <- readRDS("processed-data/06_pseudobulk/custom_cluster/lmFit-list_custom-cluster_covars-condition-sex-nspots-pc3.rda")
enrichList = lapply(names(resList), function(x) {
  tmp = resList[[x]]
  tmp = eBayes(tmp)
  tmp = topTable(tmp, "res", #p.value=.05, 
                 n=Inf, sort.by="none")
  colnames(tmp) = paste(colnames(tmp), x, sep="_")
  tmp$gene_id = rownames(tmp)
  tmp$gene_name = rdata[rownames(tmp),"gene_name"]
  return(tmp)
})

t_sm <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^t_", colnames(x), value=T),"gene_id","gene_name")]))
t_sm <- t_sm[,grep("^t_", colnames(t_sm))]
colnames(t_sm) <- gsub("^t_", "", colnames(t_sm))

lf_sm <- do.call(cbind, lapply(enrichList, function(x) x[,c(grep("^logFC_", colnames(x), value=T),"gene_id","gene_name")]))
lf_sm <- lf_sm[,grep("^logFC_", colnames(lf_sm))]
colnames(lf_sm) <- gsub("^logFC_", "", colnames(lf_sm))

#compile results in long form dframe
enrich.df = do.call(rbind, lapply(enrichList, function(x) {
  c1 = strsplit(colnames(x)[1:6], "_")
  cell.type = c1[[1]][2]
  colnames(x)[1:6] = sapply(c1, function(y) y[[1]])
  x$custom_cluster = cell.type
  rownames(x) <- NULL
  return(x)
}))
enrich.df$adj.P.Val_bin = cut(enrich.df$adj.P.Val, breaks=c(0,1e-10,.05,1), include.lowest=T)
table(enrich.df$adj.P.Val_bin, useNA="ifany")
cat("\n\n")
enrich.df$adj.P.Val_bin2 = as.character(factor(enrich.df$adj.P.Val_bin, levels=levels(enrich.df$adj.P.Val_bin), labels=c("highly sig.","sig.","NS")))
enrich.df$sig_group = ifelse(enrich.df$logFC>1, "large effect", enrich.df$adj.P.Val_bin2)
table(enrich.df[,c("sig_group","adj.P.Val_bin2")], useNA="ifany")
cat("\n\n")


enrich.df = left_join(enrich.df, rdata[,c("gene_id","gene_type")], by=c("gene_id")) %>%
  mutate(gene_type_ptn = gene_type=="protein_coding")

#save results
write.csv(t_sm, "processed-data/06_pseudobulk/custom_cluster/layer-enrichment_custom-cluster_t-stat.csv", row.names=T)
write.csv(lf_sm, "processed-data/06_pseudobulk/custom_cluster/layer-enrichment_custom-cluster_logFC.csv", row.names=T)
write.csv(enrich.df, "processed-data/06_pseudobulk/custom_cluster/layer-enrichment_custom-cluster_all-results.csv", row.names=F)

cat("\n\nLayer enrichment results saved to:",
	"\n>> processed-data/06_pseudobulk/custom_cluster/layer-enrichment_custom-cluster_t-stat.csv",
	"\n>> processed-data/06_pseudobulk/custom_cluster/layer-enrichment_custom-cluster_logFC.csv",
	"\n>> processed-data/06_pseudobulk/custom_cluster/layer-enrichment_custom-cluster_all-results.csv\n\n")


cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
