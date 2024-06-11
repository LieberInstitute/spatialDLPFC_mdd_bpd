setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scran)
	library(scater)
	library(ggplot2)
	library(here)
})
set.seed(123)

load(here("processed-data","03_QC","spe_demo-filt-outliers.Rdata"))

spe <- computeLibraryFactors(spe)
p1 = ggplot(as.data.frame(colData(spe)), aes(sizeFactor, color=slide))+
	stat_ecdf(linewidth=1)+
	scale_color_brewer(palette="Paired")+ggtitle("sizeFactor (unscaled)")+
	theme_bw()
p2 = ggplot(as.data.frame(colData(spe)), aes(sizeFactor, color=slide))+
	stat_ecdf(linewidth=1)+
	scale_color_brewer(palette="Paired")+ggtitle("sizeFactor (log scale)")+
	scale_x_log10()+
	theme_bw()

spe <- logNormCounts(spe)
spe <- runPCA(spe)

plot.df = cbind.data.frame(colData(spe), reducedDims(spe)[["PCA"]][,1:2])
p3 <- ggplot(plot.df, aes(PC1, PC2, color=sizeFactor))+
	geom_point(size=.2)+
	scale_color_viridis_c(option="F")+
	facet_wrap(vars(slide), ncol=2)+
	theme_bw()

#pdf(here("plots", "04_preprocessing","sizeFactor_normalization.pdf"), width=12, height=8)
f1 = gridExtra::grid.arrange(p1, p2, p3, layout_matrix=rbind(c(1,3),c(2,3)))
ggsave(file=here("plots","04_preprocessing","sizeFactor_normalization.png"), plot=f1, width=12, height=8, bg="white")
#dev.off()

save(spe, file=here("processed-data","04_preprocessing","spe_norm.Rdata"))
write(c(paste("*********** Created spe_norm on",format(Sys.time(), tz="UTC"),"UTC"),
	paste("*********** Old file location:",here("processed-data","03_QC","spe_demo-filt-outliers.Rdata")),
	paste("*********** New file location:",here("processed-data","04_preprocessing","spe_norm.Rdata")),
	paste("*********** Source code:",here("code","04_preprocessing","01_normalization.r")),
	"*************","*************","*************"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
