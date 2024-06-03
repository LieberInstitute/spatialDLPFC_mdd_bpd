setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')

suppressPackageStartupMessages({
        library(SpatialExperiment)
        library(nnSVG)
        library(here)
})
set.seed(123)

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))
spe_small = filter_genes(spe[,spe$slide=="V13Y10-023"], filter_genes_ncounts = 3, filter_genes_pcspots = .5, filter_mito=T)

        cat("\n",unique(spe_small$slide),"Calculating nnSVG... ",format(Sys.time(),tz="UTC"),"\n")
        cat("\n",dim(spe_small),"\n")
        spe_nnSVG <- nnSVG(spe_small, n_threads=12)
        svg = rowData(spe_nnSVG)
        cat("\nSaving results...\n")
        write.csv(svg, here("processed-data","04_preprocessing","nnSVG_V13Y10-023.csv"), row.names=T)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
