setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(ggrastr)
})
cpList <- readRDS("plots/colorPalettes.rds")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons


source("code/09_DEG_GRN/load_DEGs.r")

plot.df = filter(de.df, source=="se", cluster=="L-A") %>%
	mutate(sex.group=factor(coef, levels=comparisons))

p1 <- ggplot(plot.df, aes(x=logFC, y=-log10(adj.P.Val)))+
  rasterize(geom_point(size=.1, color="grey"), dpi=300)+
  rasterize(geom_point(data=filter(plot.df, adj.P.Val2<.05), size=.1, color="black"), dpi=300)+
  facet_wrap(vars(sex.group), ncol=2)+
  scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10))+
  coord_cartesian(xlim=c(-3,3))+
  theme_minimal()+theme(panel.grid.minor=element_blank(), panel.grid.major = element_line(linewidth=.3),
                        aspect.ratio=1, legend.position="none",
        text=element_text(size=6), strip.text.y.right=element_text(angle=0),
                        panel.border=element_rect(fill=NA, color="grey"), axis.ticks=element_line(color="grey", linewidth=.3),
        axis.title.x=element_blank())

ggsave(file="plots/publication/Figure2/main_volcano.pdf", p1,
	height=3, width=2)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
