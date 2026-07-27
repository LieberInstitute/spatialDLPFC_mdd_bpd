setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})

cpList <- readRDS("plots/colorPalettes.rds")

# load DEGs
source("code/09_DEG_GRN/load_DEGs.r")

plot.df = filter(sig.df, source=="sm") %>%
	mutate(facet_group=factor(cluster, levels=c("L-A", "L1", "L2", "L3.4", "L5", "L6", "WM"),
	labels=c("all","L1","L2","L3/4","L5","L6","WM")))

# load mean ratio results
mratio.sn = read.csv("processed-data/06_pseudobulk/SZBDMulti-seq/SZBD-control_azimuth-super-broad_mean-ratio.csv") %>%
	group_by(gene_name) %>% slice_max(MeanRatio, n=1) %>%
	mutate(cellType.target= ifelse(MeanRatio<1.5, "multi", cellType.target))

color.palette = c("Vasc"=cpList$low.res.light[["Micro.Vasc"]],
            "Micro"="#C28658",
            cpList$low.res.light[c("Astro")], cpList$low.res.bright["Oligo"],
            "InhN"=cpList$low.res.light[["Inhb"]],
            "ExcN"=cpList$low.res.light[["L2"]],
            "multi"="grey"
)

# join and keep only ones present in SZBDMulti-seq dataset
plot.df = left_join(plot.df, mratio.sn[,c("gene_name","gene","MeanRatio","cellType.target")], by=c("gene_name","gene_id"="gene")) %>%
	filter(!is.na(cellType.target)) %>%
	mutate(cellType.target=factor(cellType.target, levels=c("multi","ExcN","InhN","Micro","Oligo","Astro","Vasc")))
colSums(is.na(plot.df))


plot.df = distinct(plot.df, facet_group, cellType.target, gene_name) %>% 
	group_by(facet_group, cellType.target) %>% tally()

# save as bar plot
p2 <- ggplot(plot.df, aes(x=facet_group, y=n, fill=cellType.target))+
       geom_bar(stat="identity", position="fill", color="white", linewidth=.3)+
        scale_fill_manual(values=color.palette)+
#        coord_polar(theta="y")+
#        facet_wrap(vars(facet_group), scales="free")+
        theme_minimal()+theme(text=element_text(size=6), axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		legend.position="none")

ggsave(file="plots/publication/supp_mean-ratio/mean-ratio_bar_domain-SP.pdf", p2,
	height=3, width=5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

