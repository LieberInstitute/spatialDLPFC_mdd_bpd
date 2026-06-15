setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
})

cpList = readRDS("plots/colorPalettes.rds")
col.pal = c(cpList$transfer.bright, cpList$smoothed.bright)
# remove extra L5 and L6
col.pal = col.pal[c(1:11,14)]

spe_pseudo <- readRDS("processed-data/publication/iSEE_pseudobulk-spe_both-annotations.rds")

df = as.data.frame(colData(spe_pseudo)) %>%
	group_by(condition, sex, annotation, domain) %>% summarise(nspots=sum(nspots)) %>%
	mutate(domain= factor(domain, levels=rev(c("Micro.Vasc","Astro","L1","L2","L2.3","L3.4","L4","Inhb","L5","L6","WM","Oligo"))),
		annotation=factor(annotation, levels=c("domain-SP","domain-CT")),
		cond_sex=factor(paste(condition, sex), levels=rev(c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M"))))


p1 <- ggplot(df, aes(x=nspots, y=cond_sex, fill=domain))+
	geom_bar(stat="identity", position="stack")+
	scale_fill_manual(values=col.pal)+
	facet_grid(rows=vars(annotation))+
	scale_x_continuous(breaks=c(0,25000,50000,75000),labels=c("0","25k","50k","75k"))+
	theme_bw()+theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
		text=element_text(size=6), legend.position="none")
ggsave(file="plots/publication/Figure1/dx-sex_domain_barplot.pdf", p1, width=1.5, height=2) 
