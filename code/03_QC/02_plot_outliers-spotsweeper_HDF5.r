setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(ggplot2)
        library(dplyr)
        library(here)
})

cdata = read.csv("processed-data/03_QC/spe_n120_edge-detection_spotsweeper_colData.csv", row.names=1)

cdata2 = group_by(cdata, slide, array, sample_id) %>% 
	summarise(n_umi=sum(umi_local.outlier, na.rm=T), n_genes=sum(genes_local.outlier, na.rm=T), n_chrM.ratio=sum(chrM.ratio_local.outlier, na.rm=T)) %>%
	tidyr::pivot_longer(c("n_umi","n_genes","n_chrM.ratio"), names_to="outlier_type", names_prefix="n_", values_to="n_spots")
ymax = group_by(cdata2, sample_id) %>% summarise(ymax=sum(n_spots, na.rm=T)) %>% pull(ymax) %>% max()

#make list of lists so that 12 slides per pdf page
slideList = list(unique(cdata$slide)[1:12], unique(cdata$slide)[13:24], unique(cdata$slide)[25:length(unique(cdata$slide))])
slideList = lapply(slideList, function(x)
        unlist(lapply(x, function(y)
                unique(cdata[cdata$slide==y,"sample_id"])
        ))
)

plotList = lapply(slideList, function(x) {
	tmp = filter(cdata2, sample_id %in% x)
	ggplot(tmp, aes(x=array, y=n_spots, fill=outlier_type))+geom_bar(stat="identity", position="stack", color="grey30")+facet_wrap(vars(slide), ncol=3)+ylim(0,ymax+3)+labs(y="# outlier spots", x="slide position")+theme_bw()
})

pdf(here("plots","03_QC","spotsweeper_outliers.pdf"), height=11, width=8.5)
        plotList[[1]]
	plotList[[2]]
	plotList[[3]]
dev.off()
cat("\nplot destination:",here("plots", "03_QC", "spotsweeper_outliers.pdf"),"\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
proc.time()
options(width = 120)
sessionInfo()
