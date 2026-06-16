setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  #library(HDF5Array)
  #library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(scater)
})
set.seed(123)

cpList <- readRDS("plots/colorPalettes.rds")
#load in dlpfc marker genes (manually curated)
#source("code/06_pseudobulk/dlpfc_genes.r")
#dlpfc.genes = dlpfc.genes[c(1:5,9,6:8)]
dlpfc.genes = list("Vasc"=c("COL1A2","CLDN5"),
	"Micro"=c("GPR34","CX3CR1"),
	"Astro"=c("GJA1","AQP4"),
	"L2"=c("HPCAL1","PCDH8"),
	"L3"=c("COL5A2","NTNG1"),
	"L4"=c("TNNT2","RORB"),
	"Inhb"=c("SLC32A1","GAD1"),
	"L5"=c("PCP4","TRABD2A"),
	"L6"=c("SMIM32","NR4A2"),
	"Oligo"=c("MOBP","GJB1"))

col.pal = cpList$low.res.light
names(col.pal)[1]="Vasc"
col.pal = c(col.pal["Vasc"], "Micro"="#C28658", col.pal[c("Astro","L2","L3","L4","Inhb","L5","L6","Oligo")])
print(col.pal)

#load in sce for heatmap and dotplots
load("processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-dotplot_dx-sex-seurat-pc30.Rdata")
cond_sex = c("NTC F","NTC M","MDD F","MDD M","BPD F","BPD M")
seurat_levels= c("Micro.Vasc","Astro","L2.3","L4","Inhb","L5","L6","Oligo")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$seurat_label),
                            levels=as.character(outer(cond_sex, seurat_levels, paste)),
	labels=as.character(outer(cond_sex, c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"), paste))
)
colnames(spe_summ) <- spe_summ$sample_id

#load in dotplotDF function and format dataframe for dotplot
source("code/06_pseudobulk/custom_functions.r")
cat("\nMake dotplotDF for seurat labels...\n")
se.df = dotplotDF(spe_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
	clusters=factor(clusters, levels=levels(spe_summ$sample_id), labels=gsub(" ","\n", levels(spe_summ$sample_id))))
se.df$condition = factor(substr(se.df$clusters, start=0, stop=3), levels=c("NTC","MDD","BPD")) 
se.df$sex = factor(substr(se.df$clusters, start=5, stop=5), levels=c("F","M"))
se.df$domain = factor(substr(se.df$clusters, start=7, stop=20), levels=c("M.V","Ast","L2.3","L4","Inb","L5","L6","Olg"))
#dlpfc marker dotplot
c1 = unlist(lapply(names(dlpfc.genes), function(x) rep(x, length.out=length(dlpfc.genes[[x]]))))
se.df2 = mutate(se.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))


# load in smoothed
load("processed-data/06_pseudobulk/PRECAST_smoothed/spe_n119_pseudo-dotplot_dx-sex-smoothed-n1663-k9.Rdata")
prc_levels = c("L1","L2","L3.4","L5","L6","WM")
spe_summ$sample_id = factor(paste(spe_summ$condition, spe_summ$sex, spe_summ$smoothed_k9_1663), 
	levels=as.character(outer(cond_sex, prc_levels, paste)))
colnames(spe_summ) <- spe_summ$sample_id

# repeat for smoothed
cat("\nMake dotplotDF for precast...\n")
sm.df = dotplotDF(spe_summ, unlist(dlpfc.genes), swap_rownames="gene_name",
                  summarize_groups=F, row_data=NULL) %>%
  mutate(gene_name_f=factor(gene_name, levels=rev(unlist(dlpfc.genes))),
        clusters=factor(clusters, levels=levels(spe_summ$sample_id), labels=gsub(" ","\n", levels(spe_summ$sample_id))))
sm.df$condition = factor(substr(sm.df$clusters, start=0, stop=3), levels=c("NTC","MDD","BPD"))
sm.df$sex = factor(substr(sm.df$clusters, start=5, stop=5), levels=c("F","M"))
sm.df$domain = factor(substr(sm.df$clusters, start=7, stop=20), levels=prc_levels)
#dlpfc marker dotplot
sm.df2 = mutate(sm.df, fill_color = factor(gene_name, levels=unlist(dlpfc.genes), labels=c1))

# merge dframes
tmp = bind_rows(mutate(sm.df, annotation="domain-SP"), mutate(se.df, annotation="domain-CT")) %>%
	mutate(annotation=factor(annotation, levels=c("domain-SP", "domain-CT")),
		domain=factor(domain, levels=c("M.V","Ast","L1","L2","L2.3","L3.4","L4","Inb","L5","L6","WM","Olg")))

tmp2 = bind_rows(mutate(sm.df2, annotation="domain-SP"), mutate(se.df2, annotation="domain-CT")) %>%
	mutate(annotation=factor(annotation, levels=c("domain-SP", "domain-CT")), 
                domain=factor(domain, levels=c("M.V","Ast","L1","L2","L2.3","L3.4","L4","Inb","L5","L6","WM","Olg")))

#make labels
xlbs = levels(tmp$clusters)
## keep only the odd numbered dx indicators
odd1 = seq(1, length(xlbs), by=2)
xlbs[-odd1] = substr(xlbs[-odd1], start=4, stop=11)

p1 <- ggplot(tmp, aes(x=clusters, y=gene_name_f))+
  geom_tile(data=tmp2, aes(fill=fill_color), alpha=.5)+
  scale_fill_manual(values=col.pal, guide="none")+
  geom_count(aes(shape=sex, color=mean_expr_scaled, size=prop_spots))+
  scale_shape_manual(values=c(20,18))+
  scale_color_gradient(low="white", high="black")+
  scale_x_discrete(labels=xlbs)+
  facet_grid(cols=vars(annotation, domain), scales="free_x")+
  scale_size(range=c(.5,3), limits=c(0,1), breaks=c(0,.5,1))+
  guides(size = guide_legend(override.aes = list(shape = 20)),
	 shape = guide_legend(overrisde.aes = list(size=5)))+
  labs(color="Avg. expr.\n(scaled)", size="Prop. of\nspots")+
  theme_minimal()+theme(axis.title=element_blank(), panel.grid.major=element_blank(),
                        axis.text.y=element_text(face="italic"), legend.key.size=unit(6,"pt"), text=element_text(size=6))
#                        axis.title.y=element_text(margin=margin(0,20,0,40,"pt")))


pdf(file="plots/publication/supp_clustering_compare/dlpfc-genes_both-annotations_dx-sex-dotplot.pdf", height=4, width=6.5)
p1+theme(legend.position="none")
p1
dev.off()


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
