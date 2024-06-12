setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(gridExtra)
	library(ggplot2)
	library(here)
})

bindev.3k <- read.csv(here("processed-data","04_preprocessing","subject-biased_genes-3000.csv"), row.names = 1)
list.3k = unique(filter(bindev.3k, all.slide.outlier>5)$gene_name)

colpal = RColorBrewer::brewer.pal(n=length(unique(bindev.3k$all.slide.outlier.color)), name="Blues")

t1 = filter(bindev.3k, gene_name %in% list.3k) %>% group_by(gene_name) %>% mutate(best.rank = min(rank_brain)) %>%
	distinct(gene_name, slide, best.rank, all.slide.outlier.color) %>%
	tidyr::pivot_wider(names_from="slide", values_from="all.slide.outlier.color", values_fill=" ") %>% arrange(best.rank)

t2 = filter(bindev.3k, gene_name %in% list.3k) %>% group_by(gene_name) %>%  mutate(best.rank = min(rank_brain)) %>%
	distinct(gene_name, slide, all.slide.outlier.color, best.rank) %>%
	mutate(all.slide.outlier.color=as.character(factor(all.slide.outlier.color, 
		levels=c("[0,5]","(5,10]","(10,15]","(15,20]","(20,25]","(25,30]","(30,35]"), 
		labels=colpal))) %>% 
	arrange(best.rank) %>%
	tidyr::pivot_wider(names_from="slide", values_from="all.slide.outlier.color", values_fill="white") %>% mutate(best.rank="white") %>% as.data.frame()

gt1 = as.data.frame(ungroup(t1)[,2:8])
rownames(gt1) = t1$gene_name

color.theme = unlist(lapply(2:8, function(x) t2[,x]))
tt <- ttheme_minimal(core=list(bg_params=list(fill=color.theme)))

g1 = ggplot()+theme_minimal()+annotation_custom(tableGrob(gt1, theme=tt))

ggsave(file=here("plots","04_preprocessing","biased_features_table.png"), plot=g1, width=10, height=20, bg="white")
cat("plot destination:",here("plots","04_preprocessing","biased_features_table.png"),"\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
