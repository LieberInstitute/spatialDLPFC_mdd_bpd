setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(here)
})

l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[grep("^bindev_V.{9}_default-brain",l1)]

bindev.df = do.call(rbind, lapply(l1, function(x) {
	tmp = read.csv(here("processed-data","04_preprocessing",x))
	y=substr(x,8,17)
	tmp = filter(tmp, rank_default<2000 | rank_brain<2000)
	tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
	return(tmp)
})
)
outlier = sd(bindev.df$r.diff)*5
bindev.df = mutate(bindev.df, outlier=abs(r.diff)>outlier)

pdf(here("plots","04_preprocessing","subject-biased_genes.pdf"), height=8, width=6)
	ggplot(bindev.df, aes(x=rank_default, y=rank_brain, color=outlier))+
	geom_point()+scale_color_manual(values=c("darkgrey","red"))+
	geom_text(data=filter(bindev.df, outlier==TRUE), aes(label=gene_name), hjust=0, color="black")+
	xlim(0,2500)+
	facet_wrap(vars(slide), nrow=3)+
	labs(color=">=5 SD", title="Bin. dev. calculated per slide")+
	theme_bw()+theme(legend.position="none")
dev.off()

write.csv(filter(bindev.df, outlier==TRUE), here("processed-data","04_preprocessing","subject-biased_genes.csv"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
