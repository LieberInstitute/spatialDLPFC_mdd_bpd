setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(scran)
	library(scater)
	library(ggplot2)
	library(here)
})
set.seed(123)
load(here("processed-data","03_QC","spe_demo-filt.Rdata"))
table(spe$plot_outliers)

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

pdf(here("plots", "04_preprocessing","sizeFactor_normalization.pdf"), width=12, height=8)
gridExtra::grid.arrange(p1, p2, p3, layout_matrix=rbind(c(1,3),c(2,3)))
dev.off()

save(spe, file=here("processed-data","04_preprocessing","spe_norm.Rdata"))
write(c(paste("Created spe_norm on",format(Sys.time(), tz="UTC"),"UTC"),
	paste("File location:",here("processed-data","04_preprocessing","spe_norm.Rdata")),
	paste("Source code:",here("code","04_preprocessing","01_normalization.r")),
	"*","*","*"), here("spe_tracker_current.txt"), append=TRUE)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
