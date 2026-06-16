setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

cpList = readRDS("plots/colorPalettes.rds")

egene.df = read.csv("processed-data/11_eQTL_coloc/seurat/tables/map_cis_summary.csv")


plot.df = mutate(egene.df, cluster= factor(context, levels=rev(c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")))) %>%
	select(cluster, n_MDD_gwas_strict, n_BD_gwas_strict, n_SCZD_gwas_strict) %>%
	tidyr::pivot_longer(c("n_MDD_gwas_strict", "n_BD_gwas_strict", "n_SCZD_gwas_strict"), 
		names_to="gwas_set", values_to="n_eGenes") %>%
	mutate(gwas_set= factor(gwas_set, levels=rev(c("n_MDD_gwas_strict", "n_BD_gwas_strict", "n_SCZD_gwas_strict")),
		labels=rev(c("MDD","BPD","SCZ"))))

col.pal = c(cpList$dx.pal[2:3], "SCZ"="#27418A")

p1 <- ggplot(plot.df, aes(y=cluster, x=n_eGenes, fill=gwas_set))+
	geom_bar(stat="identity", position="stack")+
	scale_fill_manual(values=col.pal)+
	scale_y_discrete(labels=rev(c("M/V","Ast","L2/3","L4","Inb","L5","L6","Olg")))+
	labs(x="# eGenes", y="domain-CT", fill="Gene-level\noverlap")+
	theme_minimal()+theme(text=element_text(size=6), panel.grid.minor=element_blank(),
		panel.grid.major.y=element_blank(), legend.position="bottom", legend.key.size=unit(6,"pt"))


ggsave(file="plots/publication/Figure2/eGene-GWAS_gene-level_barplot.pdf", p1,
        height=3, width=2)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
