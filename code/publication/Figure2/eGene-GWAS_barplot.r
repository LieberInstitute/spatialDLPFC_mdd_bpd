setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
})

cpList = readRDS("plots/colorPalettes.rds")

egene.df = read.csv("processed-data/11_eQTL_coloc/seurat/tables/map_significant_pairs.csv.gz")


var.df = mutate(egene.df, cluster= factor(context, levels=rev(c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")))) %>%
	select(cluster, MDD_gwasVar_strict, BD_gwasVar_strict, SCZD_gwasVar_strict) %>%
	tidyr::pivot_longer(c("MDD_gwasVar_strict", "BD_gwasVar_strict", "SCZD_gwasVar_strict"), 
		names_to="gwas_set", values_to="n_vars") %>%
	filter(n_vars>0) %>%
	mutate(gwas_set= factor(gwas_set, levels=rev(c("MDD_gwasVar_strict", "BD_gwasVar_strict", "SCZD_gwasVar_strict")),
		labels=rev(c("MDD","BPD","SCZ")))) %>%
	group_by(cluster, gwas_set, .drop=F) %>% tally()
gene.df = mutate(egene.df, cluster= factor(context, levels=rev(c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")))) %>%
        select(cluster, MDD_gwasGene, BD_gwasGene, SCZD_gwasGene) %>%
        tidyr::pivot_longer(c("MDD_gwasGene", "BD_gwasGene", "SCZD_gwasGene"),
                names_to="gwas_set", values_to="n_vars") %>%
        filter(n_vars>0) %>%
        mutate(gwas_set= factor(gwas_set, levels=rev(c("MDD_gwasGene", "BD_gwasGene", "SCZD_gwasGene")),
                labels=rev(c("MDD","BPD","SCZ")))) %>%
	group_by(cluster, gwas_set, .drop=F) %>% tally()

plot.df = bind_rows(mutate(var.df, type="variant"), mutate(gene.df, type="gene")) %>%
	mutate(fill_factor= factor(paste(gwas_set, type),
		levels=rev(c("MDD gene","MDD variant","BPD gene","BPD variant","SCZ gene","SCZ variant"))))

col.pal = c("#F5C3AF", cpList$dx.pal[["MDD"]],
	"#C8AFD7", cpList$dx.pal[["BPD"]],
	"#89B6DA", "#27418A")
names(col.pal) = c("MDD gene","MDD variant","BPD gene","BPD variant","SCZ gene","SCZ variant")


p1 <- ggplot(plot.df, aes(y=cluster, x=n, fill=fill_factor))+
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
