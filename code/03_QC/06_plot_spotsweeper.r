setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(ggplot2)
        library(dplyr)
})

cdata = read.csv("processed-data/03_QC/colData_edges-problem-areas_spotsweeper.csv", row.names=1)

cdata$facet_spots = paste(cdata$round, cdata$sample_id)
cdata$facet_spots = ifelse(cdata$brnum=="Br5366", paste(cdata$facet_spots, cdata$brnum), cdata$facet_spots)

x_ordered = sort(unique(cdata$facet_spots))
cdata$x_facets = ""
cdata[cdata$facet_spots %in% x_ordered[1:30],"x_facets"] = "g1"
cdata[cdata$facet_spots %in% x_ordered[31:60],"x_facets"] = "g2"
cdata[cdata$facet_spots %in% x_ordered[61:90],"x_facets"] = "g3"
cdata[cdata$facet_spots %in% x_ordered[91:120],"x_facets"] = "g4"

cdata2 = group_by(cdata, x_facets, facet_spots) %>% 
	summarise(n_umi=sum(umi_local.outlier, na.rm=T), n_genes=sum(genes_local.outlier, na.rm=T), 
		n_chrM.ratio=sum(chrM.ratio_local.outlier, na.rm=T)) %>%
	tidyr::pivot_longer(c("n_umi","n_genes","n_chrM.ratio"), 
		names_to="outlier_type", names_prefix="n_", values_to="n_spots")

#also compute the total number of spots excluded (since some spots were flagged by multiple thresholds
cdata3 = group_by(cdata, x_facets, facet_spots) %>%
	summarise(any.outlier=sum(umi_local.outlier | genes_local.outlier | chrM.ratio_local.outlier, na.rm=T)) %>%
	tidyr::pivot_longer("any.outlier", names_to="outlier_type", values_to="n_spots")

cdata4 = bind_rows(cdata2, cdata3) %>% mutate(is_total = outlier_type=="any.outlier")

p1 <- ggplot(cdata4, aes(x=facet_spots, y=n_spots))+
	geom_bar(data=filter(cdata4, is_total==FALSE), aes(fill=outlier_type),
		stat="identity", position="stack", color="black", alpha=.5, linewidth=.5)+
	geom_bar(data=filter(cdata4, is_total==TRUE), stat="identity", fill="grey30", color="grey30", width=.5)+
	geom_text(data=filter(cdata4, is_total==TRUE), aes(y=0, label=n_spots), vjust=0, color="white", size=3)+
	#scale_y_continuous(expand=expansion(add=c(0,20)))+
	#scale_fill_manual(values=color.palette2)+
	facet_wrap(vars(x_facets), ncol=1, scales="free_x")+
	labs(x="", y="# spots", fill="QC", title="Spotsweeper outliers")+theme_bw()+
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=7),
		strip.placement = "inside", strip.text=element_blank(),
		strip.background = element_blank())

ggsave(filename="plots/03_QC/spotsweeper_outliers_barplot.png", p1, bg="white", units="in", height=12, width=9)
cat("\nPlot saved to: plots/03_QC/spotsweeper_outliers_barplot.png\n")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
