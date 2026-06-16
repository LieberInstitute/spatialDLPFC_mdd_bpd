setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})

cpList = readRDS("plots/colorPalettes.rds")
source("code/09_DEG_GRN/load_DEGs.r")

comparisons = c("F_NTC.MDD","M_NTC.MDD",
                "F_NTC.BPD","M_NTC.BPD",
                "F_MDD.BPD","M_MDD.BPD")
names(comparisons) <- comparisons

plot.df = filter(sig.df, source=="se", cluster!="L-A") %>%
	mutate(cluster= factor(cluster, levels=rev(c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo"))),
		sex.group=factor(sex.group, levels=rev(comparisons))) %>%
	group_by(sex.group, cluster) %>% tally()

col.pal = c("#F5C3AF", cpList$dx.pal[["MDD"]],
	"#C8AFD7", cpList$dx.pal[["BPD"]],
	"#BFBFBF", "#7F7F7F")
names(col.pal) = names(comparisons)

p1 <- ggplot(plot.df, aes(y=cluster, x=n, fill=sex.group))+
	geom_bar(stat="identity", position="stack")+
	scale_fill_manual(values=col.pal)+
	scale_y_discrete(labels=rev(c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg")))+
	labs(x="# DEGs", y="domain-CT")+
	theme_minimal()+theme(text=element_text(size=6), panel.grid.minor=element_blank(),
		panel.grid.major.y=element_blank(), legend.position="bottom", legend.key.size=unit(6,"pt"))


ggsave(file="plots/publication/Figure2/deg_barplot.pdf", p1,
        height=3, width=2)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
