setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(Seurat)
	library(PRECAST)
	library(dplyr)
	library(here)
})
source(here("code","04_preprocessing","helper_functions.r"))
set.seed(123)

if("srt-list_spe_counts.Rdata" %in% list.files(here("processed-data","05_clustering"))) {
        cat("\nLoading saved srt.sets object...\n")
        load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
} else {
        load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
        
        l2 = unique(spe$sample_id)
        names(l2) = lapply(l2, function(x) unique(colData(spe)[spe$sample_id==x,"brain"]))
        l2 = lapply(l2, function(x) spe[,colData(spe)$sample_id==x])

        srt.sets = lapply(l2, function(x) {
                rownames(colData(x)) <- paste(x$sample_id, rownames(colData(x)), sep="_")
                colnames(counts(x)) <- rownames(colData(x))
                colData(x)$col <- x$array_col
                colData(x)$row <- x$array_row
                count <- counts(x)
                a1 <- CreateAssayObject(count, assay = "RNA", min.features = 0, min.cells = 0)
                CreateSeuratObject(a1, meta.data = as.data.frame(colData(x)))
        })
        save(srt.sets, file=here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))
}
load(here("processed-data","05_clustering","srt-list_spe_counts.Rdata"))

# SVG based feature list
l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[
        unlist(lapply(l1, function(x) {
                if(dir.exists(here("processed-data","04_preprocessing",x))) return(FALSE)
                else {
                        split_1 = unlist(strsplit(x, split="_"))[[1]]
                        split.2 = unlist(strsplit(x, split="\\."))[[2]]
                        if(split_1=="nnSVG" & split.2=="csv") return(TRUE)
                        else {return(FALSE)}
                }
        }))
]
svg.df = do.call(rbind, lapply(l1, function(x) mutate(read.csv(here("processed-data","04_preprocessing",x)), file=x) %>% filter(padj<.05)))
svgs = unique(svg.df$gene_id)
cat("\nSVGs (per-slide nnSVG analysis, p adj<.05 in any slide):",length(svgs),"genes\n")

### layer marker based feature list
#layer.markers.100 = read.csv(here("processed-data","04_preprocessing","EXT_TableS8_sig_genes_FDR5perc_enrichment.csv")) %>%
#       filter(spatial_domain_resolution=="Sp09") %>% group_by(test) %>% slice_min(n=100, fdr)
#lm100 = unique(layer.markers.100$ensembl)
#cat("\nAdding layer markers (Huuki-Meyers top 100 fdr per cluster/test):",length(lm100),"genes\n")

### adding layer markers to svgs
#svgs = union(svgs, lm100)
#cat("\nSVGs with supplemental layer markers:",length(svgs),"genes\n")

### load in binomial deviance results run only on SVGs
brain.df = read.csv(here("processed-data","04_preprocessing","bindev_default-brain_svgs-only.csv"))

brain.df$r.diff = brain.df[,"rank_brain"]-brain.df[,"rank_default"]
brain.df$d.diff = (brain.df[,"dev_default"]-brain.df[,"dev_brain"])/brain.df$dev_brain

mean1 = mean(brain.df$r.diff)
sd1 = sd(brain.df$r.diff)
brain.df$nSD_rank = (brain.df$r.diff-mean1)/sd1

mean2 = mean(brain.df$d.diff)
sd2 = sd(brain.df$d.diff)
brain.df$nSD_dev = (brain.df$d.diff-mean2)/sd2

brain.df$dev_outlier = brain.df$nSD_dev>=5
brain.df$rank_outlier = brain.df$nSD_rank>=6

biased.genes = unique(filter(brain.df, rank_outlier==T | dev_outlier==T)$gene)
biased.genes.name = unique(filter(brain.df, rank_outlier==T | dev_outlier==T)$gene_name)

cat("\nBiased genes (bin. dev. calc. from SVGs only):",length(biased.genes),"genes\n")
sort(biased.genes.name)

### exclude biased genes from SVGs
svgs.filt = setdiff(svgs, biased.genes)
cat("\nSVGs remaining after biased gene removal:",length(svgs.filt),"genes\n")


#run precast
preobj <- CreatePRECASTObject(seuList = srt.sets, customGenelist=svgs.filt,
	premin.spots=0, premin.features=0, postmin.spots=0, postmin.features=0)
PRECASTObj <- AddAdjList(preobj, platform = "Visium")
# define model parameters and run model
PRECASTObj <- AddParSetting(PRECASTObj, maxIter = 20, verbose = TRUE, Sigma_equal=FALSE, coreNum=12)
PRECASTObj <- PRECAST(PRECASTObj, K=7)
# pick model (necessary but only changes things if more than 1 K) and integrate
PRECASTObj <- SelectModel(PRECASTObj, criteria="MBIC")
seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

#save(PRECASTObj, file=here("processed-data","05_clustering","srt_precast-list_k-8-max-it-50_24samp_bindev.3k.Rdata"))
save(seuInt,file=here("processed-data","05_clustering","srt_precast_24samp_k-7_svgs-not-loop-brain_sd-dev-5-sd-rank-6.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
