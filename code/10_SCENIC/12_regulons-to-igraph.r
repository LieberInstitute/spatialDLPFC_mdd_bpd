setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(igraph)
})

set.seed(123)

# load DEGs
source("code/09_DEG_GRN/load_DEGs.r")

# load modules for overlaps
refined.modules = read.csv("processed-data/09_DEG_GRN/spe-n119_13162-no-lowUMI_adj_with-logcounts-corr_DEG-modules-subset-refined.csv")
modules = unique(refined.modules$TF)
modList <- lapply(modules, function(x) c(x, filter(refined.modules, TF==x)$target))
names(modList) <- modules

# load regulons
# no regulons smaller than 10
regulons <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged.csv", row.names=1)
reg_subset = filter(regulons, set_size>=10)$TF

reg_long = filter(regulons, set_size>=10) %>% 
	tidyr::separate_rows(set_str, sep="/")


adj <- read.csv("processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_adj_with-logcounts-corr_top20percent.csv")

reg_adj = left_join(reg_long, adj, by=c("TF","set_str"="target"))


# make igraph source dfs
e.df = select(reg_adj, from=TF, to=set_str, weight=importance)
n.df = data.frame("node"=union(e.df$from, e.df$to),
        "is_TF"= union(e.df$from, e.df$to) %in% e.df$from,
	"is_DEG"= union(e.df$from, e.df$to) %in% unique(sig.df$gene_name))

## add info from modules
modList2 <- lapply(modList, intersect, y=n.df$node)
cat("\nNumber of module DEGs (including module source) in any regulon:\n\n")
print(sapply(modList2, length))

for (i in modules) {
	n.df[[paste0(i,"_module")]] = n.df$node %in% modList[[i]]
}


## add info of dx-sex group
dxsex.df = distinct(sig.df, sex.group, gene_name) %>% filter(gene_name %in% n.df$node) %>% mutate(is_sig=T)

n.df = left_join(n.df, tidyr::pivot_wider(dxsex.df, names_from="sex.group", values_from="is_sig"), by=c("node"="gene_name"))
n.df[is.na(n.df)] <- FALSE

cat("\nRegulon network by DEG status:\n\n")
colSums(n.df[,-1])


g <- graph_from_data_frame(e.df, directed=T, vertices=n.df)
saveRDS(g, file="processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged_min-size-10_igraph.rda")

cat("\nSaved igraph RDS object to: processed-data/10_SCENIC/spe-n109_13162-no-lowUMI_logcounts_regulons-top20_merged_min-size-10_igraph.rda\n")

cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()



