setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
})
set.seed(123)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

source("code/publication/plotting_utils.r")

mod_gene="CD74"
# mostly in importance order but modifying to be cleaner split between MDD up and BD down
plot.genes = c("CD74","CSF1R","HLA-DPA1","PLD4","LAPTM5","CX3CR1","RGS10","C3",
	"TREM2","TYROBP","SELPLG","P2RY13","GPR34","FOLR2","HLA-DRA","C3AR1","C1QA",
	"AIF1","C1QC","VSIG4","LYVE1","IL13RA1","SLC2A5",
	"RGS1","FCGR3A","C1QB","FCER1G","FCGR2A","ALOX5AP","FPR1",
	"SLC11A1","FYB1","CD14","HAMP")
length(setdiff(c(filter(refined.modules, TF==mod_gene)$target, mod_gene), plot.genes))
#plot.genes = c(filter(refined.modules, TF==mod_gene)$target, mod_gene)

#o1 = filter(refined.modules, TF==mod_gene, target %in% plot.genes) %>% arrange(desc(importance))
#o2 = filter(refined.modules, TF==mod_gene) %>% slice_max(n=1, importance) %>%
#  mutate(target=mod_gene)
#imp.df <- bind_rows(o2, o1) %>%
#  mutate(y_lab=factor(target, levels=c(rev(o1$target), mod_gene)),
#         bar_format= factor(ifelse(target==mod_gene, "self", "normal"), levels=c("self","normal")))
#
#plot.genes = imp.df$target
#length(plot.genes) 

p3 <- getDotplot(plot.genes, de.df)
p3.1 <- getMeanRatioBar(plot.genes, sce_summ)
p3.2 <- getDetectedBoxplot(rev(plot.genes), spe_summ)+coord_flip()+
	theme(axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=.5),
	axis.text.y=element_text(size=6), axis.title.x=element_blank())
#p3.3 <- ggplot(mutate(imp.df, y_lab=factor(y_lab, levels=rev(plot.genes))), 
#               aes(y=y_lab, x=importance, fill=bar_format, lty=bar_format, color=bar_format))+
#  geom_bar(stat="identity", linewidth=.3)+
#  scale_fill_manual(values=c("white","grey"), guide="none")+
#  scale_linetype_manual(values=c(2,1), guide="none")+
#  scale_color_manual(values=c("black","grey"), guide="none")+
#  labs(title=" ")+scale_y_discrete(position="right")+
#  theme_minimal()+theme(axis.text.y=element_blank(), 
#                        axis.title.y=element_blank(),
#                        axis.text.x=element_text(size=8), 
#                        panel.grid.minor=element_blank(), #panel.grid.major.y=element_blank(),
#                        plot.margin = margin(.2,.0,1.8,0, unit="cm"))

ggsave(file=paste0("plots/publication/Figure_inflamm/",mod_gene,"_dotplot.pdf"), 
       arrangeGrob(grobs=list(p3, p3.1), layout_matrix=matrix(c(1,1,1,1,1,1,1,2), ncol=8)),
       height=5.5, width=6.5)

ggsave(file=paste0("plots/publication/Figure_inflamm/",mod_gene,"_prop-detected.pdf"),
	p3.2, 
	height=2, width=6.5)

# just do median spots detected
ordered_genes = plot.genes
gids = rownames(spe_summ)[rowData(spe_summ)$gene_name %in% ordered_genes]
names(gids) = rowData(spe_summ)[gids,"gene_name"]
stopifnot(length(ordered_genes)==length(gids))

df1 = as.data.frame(assay(spe_summ, "logcounts.prop.detected")[gids,])
df2 = tidyr::pivot_longer(tibble::rownames_to_column(df1, var="gene_id"), all_of(colnames(df1)), names_to="sample", values_to="prop.spots.detected") %>%
	left_join(as.data.frame(rowData(spe_summ)[gids,c("gene_id","gene_name")])) %>%
        mutate(gene_name= factor(gene_name, levels=ordered_genes))
tmp = group_by(df2, gene_name) %>% summarise(med_prop.spots = median(prop.spots.detected))
p3.3 <- ggplot(tmp, aes(x=gene_name, y=med_prop.spots))+
	geom_bar(stat="identity", width=.8)+ylim(0,.15)+
	geom_hline(aes(yintercept=.01), color="red")+
	labs(title="", y="med. prop. spots")+
	theme_minimal()+theme(axis.text.y=element_text(size=6), axis.title.x=element_blank(),
                        panel.grid.minor=element_blank(),
                        axis.text.x=element_text(angle=90, hjust=1, vjust=.5, size=6),
                        plot.margin = margin(.2,.2,1.5,.2,"cm"))

ggsave(file=paste0("plots/publication/Figure_inflamm/",mod_gene,"_median-prop-detected.pdf"),
        p3.3,
        height=2, width=6.5)

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
