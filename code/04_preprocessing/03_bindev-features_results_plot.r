setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(SpatialExperiment)
	library(here)
})

source(here("code","04_preprocessing","03_plot_bindev_functions.r"))
l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[grep("^bindev_V.{9}_default-brain",l1)]

# perform for top 3k features
cat("Find subject-biased genes in top 3k deviant features\n")
bindev.3k = do.call(rbind, lapply(l1, function(x) {
	tmp = read.csv(here("processed-data","04_preprocessing",x))
	y=substr(x,8,17)
	tmp = filter(tmp, rank_default<=3000 | rank_brain<=3000)
	tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
	return(tmp)
})
)

f1 = plotBinDevResults(bindev.3k, 3000)
cat(">>> Saving top 3k csv...\n")
write.csv(f1[[2]], here("processed-data","04_preprocessing","subject-biased_genes-3000.csv"), row.names=F)
cat(">>> Saving rank comparison results plots...\n")
ggsave(file=here("plots","04_preprocessing","bindev-3000_hist-scatter.png"), plot=f1[[1]], height=10, width=12, bg="white")

load(here("processed-data","04_preprocessing","spe_norm.Rdata"))

f2 = plotBiasedFeatures(f1[[2]], sd.cutoff=5, spe)
cat(">>> Saving biased features (nSD>5) dotplots...\n")
ggsave(file=here("plots","04_preprocessing","bindev-3000_sd-5_biased-dotplot.png"), plot=gridExtra::grid.arrange(f2[[1]],f2[[2]], ncol=2, top="Biased = SD cutoff >5"),
	height=12, width=20, bg="white")

# perform for top 2k features
cat("\nFind subject-biased genes in top 2k deviant features\n")
bindev.2k = do.call(rbind, lapply(l1, function(x) {
        tmp = read.csv(here("processed-data","04_preprocessing",x))
        y=substr(x,8,17)
        tmp = filter(tmp, rank_default<=2000 | rank_brain<=2000)
        tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
        return(tmp)
})
)

f3 = plotBinDevResults(bindev.2k, 2000)
cat(">>> Saving top 2k csv...\n")
write.csv(f3[[2]], here("processed-data","04_preprocessing","subject-biased_genes-2000.csv"), row.names=F)
cat(">>> Saving rank comparison results plots...\n")
ggsave(file=here("plots","04_preprocessing","bindev-2000_hist-scatter.png"), plot=f3[[1]], height=10, width=12, bg="white")

f4 = plotBiasedFeatures(f3[[2]], sd.cutoff=5, spe)
cat(">>> Saving biased features (nSD>5) dotplots...\n")
ggsave(file=here("plots","04_preprocessing","bindev-2000_sd-5_biased-dotplot.png"), plot=gridExtra::grid.arrange(f4[[1]],f4[[2]], ncol=2, top="Biased = SD cutoff >5"),
        height=12, width=20, bg="white")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
