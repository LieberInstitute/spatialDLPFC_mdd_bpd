setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(UpSetR)
})

mod_gene = "PLP1"

# load modules for overlaps
refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")
modlist = c(mod_gene, filter(refined.modules, TF==mod_gene)$target)

# load regulons
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)

regulonList = filter(regulons, set_size>10) %>% 
	select(TF, set_str) %>% tibble::deframe()
regulonList = sapply(regulonList, strsplit, split="/")

# overlaps
reg1 = lapply(regulonList, intersect, y=modlist)
reg1[[mod_gene]] = modlist

#set_order = names(sort(sapply(reg1, length)))
#reg1 = reg1[set_order]
set_size = sort(sapply(reg1, length))

pdf(file="plots/publication/Figure3/PLP1-module_regulon_upset.pdf", height=4, width=4, onefile=F)
upset(fromList(reg1), sets=names(set_size)[set_size>0], keep.order=T, mb.ratio=c(.5,.5))
dev.off()

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()
