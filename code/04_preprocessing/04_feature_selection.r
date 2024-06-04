setwd('/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/')
suppressPackageStartupMessages({
	library(dplyr)
	library(here)
})

l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[grep("^bindev_V.{9}_default-brain",l1)]

bindev.df = do.call(rbind, lapply(l1, function(x) {
	tmp = read.csv(here("processed-data","04_preprocessing",x))
	y=substr(x,8,17)
	tmp = mutate(tmp, r.diff = rank_brain-rank_default, slide=y)
	return(tmp)
})
)

genes.2k = unique(filter(bindev.df, rank_brain<=2000)$gene_name)
genes.3k = unique(filter(bindev.df, rank_brain<=3000)$gene_name)

outlier.2k = unique(read.csv(here("processed-data","04_preprocessing","subject-biased_genes.csv"))$gene_name)
outlier.3k = unique(read.csv(here("processed-data","04_preprocessing","subject-biased_genes-3000.csv"))$gene_name)

l1 = list.files(here("processed-data","04_preprocessing"))
l1 = l1[
	unlist(lapply(l1, function(x) {
		if(dir.exists(here("processed-data","04_preprocessing",x))) return(FALSE)
		else {
			split_1 = unlist(strsplit(x, split="_"))[[1]]
			split.2 = unlist(strsplit(x, split="\\."))[[2]]
			if(split_1=="nnSVG" & split.2=="csv") return(TRUE)
			else {return(FALSE)}
		}
	}))
]

svg.df = do.call(rbind, lapply(l1, function(x) mutate(read.csv(here("processed-data","04_preprocessing",x)), file=x) %>% filter(padj<.05)))
svg.genes = unique(svg.df$gene_name)

gene.list = list("bd.2k"=genes.2k, "bd.3k"=genes.3k, "svg"=svg.genes)
cat("\nStarting feature list size:\n")
unlist(lapply(gene.list, length))

exclude.list <- list("MTRN"=grep("^MTRN",bindev.df$gene_name, value=T),
	"DNAJ"=grep("^DNAJ",bindev.df$gene_name, value=T),
	"HSP"=grep("^HSP",bindev.df$gene_name, value=T),
	"ribo"=grep("RPS|RPL",bindev.df$gene_name, value=T))
cat("\n\nExcluded genes:\n")
unlist(lapply(exclude.list, length))

gene.list2 = lapply(gene.list, setdiff, unlist(exclude.list))
cat("\n\nFiltered feature list size:\n")
unlist(lapply(gene.list2, length))

lm = read.csv(here("processed-data","04_preprocessing","EXT_layer-markers_tableS5.csv"))
layer.markers = filter(lm, rank<=3000)$gene_name

feature.list = list("bindev.2k"=union(setdiff(gene.list2[[1]], outlier.2k), layer.markers),
	"bindev.3k"=union(setdiff(gene.list2[[2]], outlier.3k), layer.markers),
	"svg"=union(setdiff(gene.list2[[3]], outlier.2k), layer.markers))
cat("\n\nFinal feature list size:\n")
unlist(lapply(feature.list, length))

saveRDS(feature.list, here("processed-data","04_preprocessing","bindev-2k-3k_svg_feature-list.rda"))

## Reproducibility information
print("Reproducibility information:")
format(Sys.time(), tz="UTC")
proc.time()
options(width = 120)
sessionInfo()
