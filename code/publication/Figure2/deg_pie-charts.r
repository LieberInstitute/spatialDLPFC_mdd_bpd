setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})

cpList <- readRDS("plots/colorPalettes.rds")

# load DEGs
source("code/09_DEG_GRN/load_DEGs.r")

# facet by consensus clusters
all_clusters = c("L-A sm","L-A se","Micro.Vasc se","Astro se","L1 sm",
                 "L2 sm","L2.3 se","L3.4 sm","L4 se",
                 "Inhb se","L5 sm","L5 se","L6 sm","L6 se","WM sm","Oligo se")
plot.df = mutate(sig.df, facet_group=factor(cluster_source, levels=all_clusters, 
	labels=c("L-A","L-A","M.V","Ast","L1",
		"L2/3","L2/3","L3/4","L3/4",
		"Inb","L5","L5","L6","L6","WM/O","WM/O")))

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

# plot
pbase <- ggplot(plot.df, aes(x=facet_group, y=n, fill=cellType.target))+
#	geom_bar(stat="identity", position="fill")+
	scale_fill_manual(values=color.palette)+
	coord_polar(theta="y")+
	facet_wrap(vars(facet_group), scales="free")+
	theme_void()+theme(aspect.ratio=1, legend.position="bottom")

#pdf(file="plots/publication/Figure2/deg_pie-charts.pdf")
#pbase+geom_bar(stat="identity", position="fill")
#pbase+geom_bar(stat="identity", position="fill", color="white", linewidth=.3) 
#pbase+geom_bar(stat="identity", position="fill", color="black", linewidth=.3)
#dev.off()


# this time save as bar plot
p2 <- ggplot(plot.df, aes(x=facet_group, y=n, fill=cellType.target))+
       geom_bar(stat="identity", position="fill", color="white", linewidth=.3)+
        scale_fill_manual(values=color.palette)+
#        coord_polar(theta="y")+
#        facet_wrap(vars(facet_group), scales="free")+
        theme_minimal()+theme(text=element_text(size=6), axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
		legend.position="none")

ggsave(file="plots/publication/Figure2/deg_pie-charts_bar-style.pdf", p2,
	height=3, width=5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()

