setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(dplyr)
	library(here)
})
set.seed(123)

load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
#srt.sets = srt.sets[c("Br6529","Br6436","Br8221","Br8600","Br5799","Br6192")]
#srt.sets = srt.sets[c("Br5366","Br5916","Br6231","Br5639","Br5599",'Br5694')]
feature.list = readRDS(here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))
#l1 = list.files(here("processed-data","04_preprocessing"))
#l1 = l1[
#        unlist(lapply(l1, function(x) {
#                if(dir.exists(here("processed-data","04_preprocessing",x))) return(FALSE)
#                else {
#                        split_1 = unlist(strsplit(x, split="_"))[[1]]
#                        split.2 = unlist(strsplit(x, split="\\."))[[2]]
#                        if(split_1=="nnSVG" & split.2=="csv") return(TRUE)
#                        else {return(FALSE)}
#                }
#        }))
#]
#svg.df = do.call(rbind, lapply(l1, function(x) mutate(read.csv(here("processed-data","04_preprocessing",x)), file=x) %>% filter(padj<.05)))
#remove.svgs = distinct(svg.df, file, gene_id, gene_name) %>% group_by(gene_name, gene_id) %>% tally() %>% filter(n==1) %>% pull(gene_id)
bad.gene.list = readRDS(here("processed-data","04_preprocessing","biased_excluded_feature-list.rda"))
svg.nobias = setdiff(feature.list$svg, bad.gene.list[[1]])
#svg.more1 = setdiff(svg.nobias, remove.svgs)
#layer.markers = read.csv(here("processed-data","04_preprocessing","EXT_TableS8_sig_genes_FDR5perc_enrichment.csv"))
#lm = filter(layer.markers, spatial_domain_resolution=="Sp09") %>% group_by(test) %>% slice_min(n=100,fdr)

preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=svg.nobias,
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj <- PRECAST(PRECASTObj, K=8)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

#save(PRECASTObj, file=here("processed-data","05_clustering","srt_precast-list_k-8-max-it-50_24samp_bindev.3k.Rdata"))
save(seuInt,file=here("processed-data","05_clustering","srt_precast_24samp_k-8_svg-nobias.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
