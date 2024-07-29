setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scry)
	library(dplyr)
	library(parallel)
	library(here)
})
set.seed(123)

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))

# run with SVGs only: first load SVGs
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
length(svgs)

l2 = unique(spe$slide)
names(l2) = l2
#l2 = lapply(l2, function(x) spe[-grep("^MT-",rowData(spe)$gene_name), spe$slide==x])
# run with SVGs only: filter spe to SVGs
l2 = lapply(l2, function(x) spe[svgs, spe$slide==x])

mclapply(l2, function(x) {
	cat("\n",unique(x$slide),": Running default model...\n")
        default <- devianceFeatureSelection(x, fam="binomial", batch=NULL)

        df = cbind.data.frame("gene"=rownames(default),"gene_name"=rowData(default)$gene_name,
                "dev"= rowData(default)$binomial_deviance,
                "rank"=(nrow(default)+1)-rank(rowData(default)$binomial_deviance))
	
	cat("\n",unique(x$slide),": Running batch model...\n")
        batch <- devianceFeatureSelection(x, fam="binomial", batch=as.factor(x$brain))

        df = left_join(df, cbind.data.frame("gene"=rownames(rowData(batch)), "dev"=rowData(batch)$binomial_deviance,
                "rank"=(nrow(batch)+1)-rank(rowData(batch)$binomial_deviance)),
        by="gene", suffix=c("_default","_brain"))

        write.csv(df, here("processed-data","04_preprocessing",paste0("bindev_",unique(x$slide),"_default-brain_svgs-only.csv")),row.names=FALSE)
}, mc.cores=6)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
