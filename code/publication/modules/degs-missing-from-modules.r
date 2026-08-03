setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(SpatialExperiment)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	library(clusterProfiler)
})

source("code/publication/plotting_utils.r")
mbv.degs = unique(sig.df$gene_name)

refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")

missing.degs = setdiff(mbv.degs, union(refined.modules$TF, refined.modules$target))
length(missing.degs)

# use sce_summ to order degs
gids = rownames(sce_summ)[rowData(sce_summ)$gene_name %in% missing.degs]
names(gids) = rowData(sce_summ)[gids,"gene_name"]

df2 = tibble::rownames_to_column(as.data.frame(assay(sce_summ, "logcounts.mean")[gids,]))
colnames(df2) = c("gene_id", as.character(colData(sce_summ)$azimuth_super.broad))

bar.df = tidyr::pivot_longer(df2, all_of(levels(colData(sce_summ)$azimuth_super.broad)), 
                             names_to="cellType", values_to="mean.expr") %>%
  left_join(as.data.frame(rowData(sce_summ)[gids,c("gene_id","gene_name")])) %>%
  mutate(cellType=factor(cellType, levels=c("Vasc","Astro","Oligo","Micro","InhN","ExcN")))

order1 = mutate(bar.df, is_nrn = cellType %in% c("InhN","ExcN")) %>% group_by(gene_name, is_nrn) %>%
  summarise(sum_mean= sum(mean.expr)) %>%
  tidyr::pivot_wider(names_from="is_nrn", values_from="sum_mean", names_prefix = "nrn_") %>%
  mutate(nrn.ratio= nrn_TRUE/nrn_FALSE) %>% arrange(nrn.ratio) %>%
  pull(gene_name)
order2 = order1[c(grep("HB",order1), setdiff(1:length(order1), grep("HB",order1)))]
order3 = c(setdiff(order2, c("SST","CORT","CRH","VGF")), c("VGF","CORT","SST","CRH"))
order4 = rev(c(setdiff(missing.degs, order3), order3))


# plot missing DEGs
p0 <- getDotplot(order4, de.df)+theme(axis.text.y=element_text(size=6))
p1 <- getMeanRatioBar(order4, sce_summ)
#p2 <- getDetectedBoxplot(order4, spe_summ)

ggsave(file="plots/publication/modules/supp_missing-degs_dotplot.pdf",
       arrangeGrob(grobs=list(p0,p1), layout_matrix=matrix(c(1,1,1,1,1,1,1,2), ncol=8), top=NULL),
       height=7, width=6.5)


stop("Early stopping: Just need dotplot")

# now GO ORA dotplot
ora_modules = readRDS("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_modules-DEG-subset-refined_ORA-GO-results.rda")

select.terms = list("A2M"=c("GO:0150104","GO:0003018","GO:2001214"), #A2M
                    "IFITM3"=c("GO:0062023","GO:0002526","GO:0019955"), #IFITM3
                    "CD74"=c("GO:0098883","GO:0001774","GO:0001818"), #CD74
                    "HSPA1A"=c("GO:0048545","GO:0034599","GO:0034612"), #HSPA1A
                    "MT1M"=c("GO:0010038","GO:0098754","GO:0045089"), #MT1M
                    "SNHG14"=c("GO:0001217","GO:0042393","GO:0141108"), #SNHG14
                    "GLUL"=c("GO:0036293","GO:0070374","GO:0051384"), #GLUL
                    "GAD1"=c("GO:0098982","GO:0033555","GO:0160041"), #GAD1
                    "CAMK2N1"=c("GO:0008331"), #CAMK2N1
                    "GRIN1"=c("GO:0048167","GO:0042752","GO:0014069"), #GRIN1
                    "PRKAR1A"=c("GO:0004860","GO:0031625"), #PRKAR1A
                    "COX4I1"=c("GO:0046034","GO:0004129","GO:0032543"),#COX4I1
                    "EEF1A1"=c("GO:0002181","GO:0043484"),#EEF1A1
                    "GFAP"=c("GO:0072331","GO:0010506","GO:1904018"), #GFAP
                    "PLP1"=c("GO:0008366","GO:0043209")#, #PLP1
)

any.term = unlist(select.terms)
#length(any.term) #40
#length(unique(any.term)) #40

# plot any term highlighted in main panel
any.term.missing = intersect(any.term, rownames(ora_modules[["missing"]]@result))
# additional terms from missing ORA that are helpful
plot.terms = union(rev(any.term.missing), c("GO:0005179","GO:0001664","GO:0033549","GO:0005833","GO:0101002"))

mod_subset = c("SNHG14","PRKAR1A","COX4I1","EEF1A1",
               "PLP1",
               "GFAP","GLUL","HSPA1A","MT1M",
               "IFITM3","CD74","A2M",
               "GRIN1","CAMK2N1","GAD1",
               "missing")
names(mod_subset) = mod_subset

plot.df = do.call(rbind, lapply(mod_subset, function(x) {
  tmp1 = ora_modules[[x]]@result
  tmp1 = tmp1[intersect(plot.terms, rownames(tmp1)),] %>% mutate(module=x)
  return(tmp1)
}))

plot.df = mutate(plot.df, ID= factor(ID, levels=rev(plot.terms)), 
                 module= factor(module, levels=mod_subset)) %>%
  arrange(ID)
plot.df$Description2 = factor(plot.df$Description, levels=unique(plot.df$Description))

p2 <- ggplot(plot.df, aes(y=Description2, x=module, size=Count, color=log2(FoldEnrichment)))+
  geom_count()+
  scale_color_gradient("log2\nFold\nEnrich.", low="white", high="black", limits=c(0,8))+
  scale_size("# DEGs", range=c(2,6), breaks=c(3,9,15), limits=c(2,19))+
  labs(y="")+
  theme_bw()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                   strip.text.y=element_text(angle=0))

ggsave(file="plots/publication/modules/supp_missing-degs_GO-ORA.pdf", p2, 
       width=5, height=3)


## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="EST")
proc.time()
options(width = 120)
sessionInfo()
