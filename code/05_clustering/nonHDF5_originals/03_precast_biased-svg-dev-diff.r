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

# load in binomial deviance results run only on SVGs
l3 = list.files(here("processed-data","04_preprocessing"))
l3 = l3[grep("^bindev_V.{9}_default-brain_svgs",l3)]

bindev.svg = do.call(rbind, lapply(l3, function(x) {
	tmp = read.csv(here("processed-data","04_preprocessing",x))
        y=substr(x,8,17)
        tmp = mutate(tmp, slide=y)
        return(tmp)
})
)
svg.prep <- prepBias(bindev.svg)
svg.bias <- findBiasedFeatures(svg.prep, sd.safe="[0,5)")
biased.genes1 = intersect(filter(svg.bias, nSD.outlier==T)$gene, filter(svg.bias, nSD.outlier_slide==T)$gene)
biased.genes.name1 = intersect(filter(svg.bias, nSD.outlier==T)$gene_name, filter(svg.bias, nSD.outlier_slide==T)$gene_name)
cat("\nBiased genes (bin. dev. calc. from SVGs only):",length(biased.genes1),"genes\n")
sort(biased.genes.name1)

#other approach
bindev.svg$d.diff = bindev.svg$dev_default-bindev.svg$dev_brain
bindev.svg$d.diff2 = bindev.svg$d.diff/bindev.svg$dev_brain
biased.df = group_by(bindev.svg, gene, gene_name) %>% summarise(d.diff2=max(d.diff2)) %>% filter(d.diff2>=.5)
biased.genes2 = biased.df$gene
biased.genes.name2 = biased.df$gene_name
cat("\nBiased genes (bin. dev. calc. from SVGs only)(% diff deviance):",length(biased.genes2),"genes\n")
sort(biased.genes.name2)

# exclude biased genes from SVGs
biased.genes = union(biased.genes1, biased.genes2)
biased.genes.name = union(biased.genes.name1, biased.genes.name2)
cat("\nTotal biased genes (union list 1 and list2):",length(biased.genes),"genes\n")
sort(biased.genes.name)

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
save(seuInt,file=here("processed-data","05_clustering","srt_precast_24samp_k-7_svgs-svg-bias3.Rdata"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
