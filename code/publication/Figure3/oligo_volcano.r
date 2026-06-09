setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(ggrastr)
})
set.seed(123)

source("code/09_DEG_GRN/load_DEGs.r")

olg.df = filter(de.df, sex.group %in% c("F_NTC.MDD","F_NTC.BPD"), cluster=="Oligo") %>%
	mutate(sex.group= factor(sex.group, levels=c("F_NTC.MDD","F_NTC.BPD")))

# use this to check that filtered correctly
group_by(olg.df, sex.group) %>% tally()

plot.genes = c("BCL6","HSPA1B","ANP32B","UBA52","TPT1")

p1 <- ggplot(olg.df, aes(x=logFC, y=-log10(adj.P.Val)))+
	rasterize(geom_point(size=.1, color="grey50"), dpi=300)+
	geom_point(data=filter(olg.df, gene_name %in% plot.genes), size=.3, color="black")+
	ggrepel::geom_text_repel(data=filter(olg.df, gene_name %in% plot.genes), aes(label=gene_name),
		size=2, min.segment.length=0)+
	facet_wrap(vars(sex.group), ncol=2)+xlim(-2.5,2.5)+
	theme_bw()+theme(text=element_text(size=6), panel.grid.minor=element_blank())

ggsave(file="plots/publication/Figure3/main_Oligo_volcano.pdf", p1, width=3.5, height=2)

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
